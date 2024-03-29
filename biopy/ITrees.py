
#
# Trees.py 
#
# Copyright 2005 by Frank Kauff & Cymon J. Cox. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Tree class handles phylogenetic trees. Provides a set of methods to read and write newick-format tree
# descriptions, get information about trees (monphyly of taxon sets, congruence between trees, common ancestors,...)
# and to manipulate trees (reroot trees, split terminal nodes).
#
# Bug reports welcome: fkauff@duke.edu
#

import sys, random
#from Bio.Nexus import Nodes
import Nodes

PRECISION_BRANCHLENGTH=6
PRECISION_SUPPORT=6

class TreeError(Exception): pass

class NodeData:
    """Stores tree-relevant data associated with nodes (e.g. branches or otus)."""
    def __init__(self,taxon=None,branchlength=0.0,support=None):
        self.taxon=taxon
        self.branchlength=branchlength
        self.support=support

class Tree(Nodes.Chain):
    """Represents a tree using a chain of nodes with on predecessor (=ancestor)
    and multiple successors (=subclades).
    """ 
    # A newick tree is parsed into nested list and then converted to a node list in two stages
    # mostly due to historical reasons. This could be done in one swoop). Note: parentheses ( ) and
    # colon : are not allowed in taxon names. This is against NEXUS standard, but makes life much
    # easier when parsing trees.
    
    ## NOTE: Tree should store its data class in something like self.dataclass=data,
    ## so that nodes that are generated have easy access to the data class
    ## Some routines use automatically NodeData, this needs to be more concise

    def __init__(self,tree=None,weight=1.0,rooted=False,name='',data=NodeData,values_are_support=False,max_support=1.0):
        """Ntree(self,tree)."""
        Nodes.Chain.__init__(self)
        self.dataclass=data
        self.__values_are_support=values_are_support
        self.max_support=max_support
        self.weight=weight
        self.rooted=rooted
        self.name=name
        root=Nodes.Node(data())
        self.add(root)
        self.root=root.id
        if tree:    # use the tree we have
            # if Tree is called from outside Nexus parser, we need to get rid of linebreaks, etc
            tree=tree.strip().replace('\n','').replace('\r','')
            # there's discrepancy whether newick allows semicolons et the end
            tree=tree.rstrip(';')

            if tree.count('(') != tree.count(')'):
              raise TreeError, 'Parentheses do not match in tree: '+tree

            st = self._parse(tree)
            if isinstance(st[1][-1], dict) :
                root.data.attributes = st[1].pop(-1)

            self._add_subtree(parent_id=root.id,tree=st[0])
        
    def _parse(self,tree):
        """Parses (a,b,c...)[[[xx]:]yy] into subcomponents and travels down recursively."""

        #Remove any leading/trailing white space - want any string starting
        #with " (..." should be recognised as a leaf, "(..."
        tree = tree.strip()
        #assert not tree[0].isspace()
        if tree[0] != '(': # a leaf
            lb = tree.find('[')
            colon = tree.rfind(':')   

            if lb > -1 :
                values = self._get_values(tree[lb:])
                tree = tree[:lb if colon==-1 else min(lb,colon)]
            else :
                if colon > -1:
                    values = self._get_values(tree[colon+1:])
                    tree = tree[:colon]
                else :
                    values = [None]
            return [tree, values]
        else:
            closing=tree.rfind(')')

            val = self._get_values(tree[closing+1:])
            ## if not val:
            ##     val=[None]
            subtrees=[]
            plevel=0
            prev=1
            for p in range(1,closing):
                if tree[p] in "[(":
                    plevel+=1
                elif tree[p] in "])":
                    plevel-=1
                elif tree[p]==',' and plevel==0:
                    subtrees.append(tree[prev:p])
                    prev=p+1
            subtrees.append(tree[prev:closing])
            subclades=[self._parse(subtree) for subtree in subtrees]
            return [subclades,val]
    
    def _add_subtree(self,parent_id=None,tree=None):
        """Adds leaf or tree (in newick format) to a parent_id. (self,parent_id,tree)."""
        
        if parent_id is None:
            raise TreeError('Need node_id to connect to.')
        for st in tree:
            nd=self.dataclass()

            if len(st)>1:
                if isinstance(st[1][-1], dict) :
                    nd.attributes = st[1].pop(-1)

                if len(st[1])>=2: # if there's two values, support comes first. Is that always so?
                    nd.support=st[1][0]
                    if st[1][1] is not None:
                        nd.branchlength=st[1][1]
                elif len(st[1])==1: # otherwise it could be real branchlengths or support as branchlengths
                    if not self.__values_are_support: # default
                        if st[1][0] is not None:
                            nd.branchlength=st[1][0]
                    else:
                        nd.support=st[1][0]
            
            if type(st[0])==list: # it's a subtree
                sn=Nodes.Node(nd)
                self.add(sn,parent_id)
                self._add_subtree(sn.id,st[0])
            else: # it's a leaf
                nd.taxon=st[0].strip()
                leaf=Nodes.Node(nd)
                self.add(leaf,parent_id)
    
    def _get_values(self, text):
        """Extracts values (support/branchlength) from xx[:yyy], xx."""

        if text=='':
            return [None]
        text = text.strip()
        if text[0] == '[' :
            lcomment = _getStuff(text[1:], ']')
            end = 1+lcomment
            comment = text[1:end]
            attr = _parseAttributes(comment)
            values = [attr,]
            text = text[end+1:]
        elif text[0] == ':' and text[1:].lstrip().startswith('['):
            c0 = 1 + text[1:].find('[')
            lcomment = _getStuff(text[c0:], ']')
            end = c0+lcomment
            comment = text[c0+1:end]
            attr = _parseAttributes(comment)
            values = [attr,]
            text = ':' + text[end+1:]
        else:
            values = []

        return [float(t) for t in text.split(':') if t.strip()] + values
   
    def _walk(self,node=None):
        """Return all node_ids downwards from a node."""
        
        if node is None:
            node=self.root
        for n in self.node(node).succ:
            yield n
            for sn in self._walk(n):
                yield sn

    def node(self,node_id):
        """Return the instance of node_id.
        
        node = node(self,node_id)
        """
        if node_id not in self.chain:
            raise TreeError('Unknown node_id: %d' % node_id)
        return self.chain[node_id]

    def split(self,parent_id=None,n=2,branchlength=1.0):
        """Speciation: generates n (default two) descendants of a node.
        
        [new ids] = split(self,parent_id=None,n=2,branchlength=1.0):
        """ 
        if parent_id is None:
            raise TreeError('Missing node_id.')
        ids=[]
        parent_data=self.chain[parent_id].data
        for i in range(n):
            node=Nodes.Node()
            if parent_data:
                node.data=self.dataclass()
                # each node has taxon and branchlength attribute
                if parent_data.taxon:
                    node.data.taxon=parent_data.taxon+str(i)
                node.data.branchlength=branchlength
            ids.append(self.add(node,parent_id))
        return ids

    def search_taxon(self,taxon):
        """Returns the first matching taxon in self.data.taxon. Not restricted to terminal nodes.
        
        node_id = search_taxon(self,taxon)
        """
        for id,node in self.chain.items():
            if node.data.taxon==taxon:
                return id
        return None
   
    def prune(self,taxon):
        """Prunes a terminal taxon from the tree.
        
        id_of_previous_node = prune(self,taxon)
        If taxon is from a bifurcation, the connectiong node will be collapsed
        and its branchlength added to remaining terminal node. This might be no
        longer a meaningful value'
        """
        
        id=self.search_taxon(taxon)
        if id is None:
            raise TreeError('Taxon not found: %s' % taxon)
        elif id not in self.get_terminals():
            raise TreeError('Not a terminal taxon: %s' % taxon)
        else:
            prev=self.unlink(id)
            self.kill(id)
            if not prev==self.root and len(self.node(prev).succ)==1:
                succ=self.node(prev).succ[0]
                new_bl=self.node(prev).data.branchlength+self.node(succ).data.branchlength
                self.collapse(prev)
                self.node(succ).data.branchlength=new_bl
            return prev
        
    def get_taxa(self,node_id=None,asIDs = False):
        """Return a list of all otus downwards from a node (self, node_id).

        nodes = get_taxa(self,node_id=None)
        """

        if node_id is None:
            return [self.chain[node_id].data.taxon for node_id in self.get_terminals()]
        
        ## if node_id is None:
        ##     node_id=self.root
        ## 
        if node_id not in self.chain:
            raise TreeError('Unknown node_id: %d.' % node_id)
        
        tx = []
        waiting = [self.chain[node_id]]
        while len(waiting) :
            n = waiting.pop(0)
            if n.succ==[] :
                tx.append(n.id if asIDs else n.data.taxon)
            else :
                waiting.extend([self.chain[x] for x in n.succ])
        return tx
        
        ## if self.chain[node_id].succ==[]:
        ##     if self.chain[node_id].data:
        ##         return [self.chain[node_id].data.taxon]
        ##     else:
        ##         return None
        ## else:
        ##     list=[]
        ##     for succ in self.chain[node_id].succ:
        ##         list.extend(self.get_taxa(succ))
        ##     return list

    def get_terminals(self):
        """Return a list of all terminal nodes."""
        return [i for i in self.all_ids() if self.node(i).succ==[]]

    def sum_branchlength(self,root=None,node=None):
        """Adds up the branchlengths from root (default self.root) to node.
        
        sum = sum_branchlength(self,root=None,node=None)
        """

        if root is None:
            root=self.root
        if node is None:
            raise TreeError('Missing node id.')
        blen=0.0
        while node is not None and node is not root: 
            blen+=self.node(node).data.branchlength
            node=self.node(node).prev
        return blen

    def set_subtree(self,node):
        """Return subtree as a set of nested sets.

        sets = set_subtree(self,node)
        """
        
        if self.node(node).succ==[]:
            return self.node(node).data.taxon
        else:
            return set([self.set_subtree(n) for n in self.node(node).succ])
            
    def is_identical(self,tree2):
        """Compare tree and tree2 for identity.

        result = is_identical(self,tree2)
        """
        return self.set_subtree(self.root)==tree2.set_subtree(tree2.root)

    def is_compatible(self,tree2,threshold,strict=True):
        """Compares branches with support>threshold for compatibility.
        
        result = is_compatible(self,tree2,threshold)
        """

        # check if both trees have the same set of taxa. strict=True enforces this.
        missing2 = set(self.get_taxa()) - set(tree2.get_taxa())
        missing1 = set(tree2.get_taxa()) - set(self.get_taxa())
        if strict and (missing1 or missing2):
            if missing1: 
                print 'Taxon/taxa %s is/are missing in tree %s' % (','.join(missing1) , self.name)
            if missing2:
                print 'Taxon/taxa %s is/are missing in tree %s' % (','.join(missing2) , tree2.name)
            raise TreeError, 'Can\'t compare trees with different taxon compositions.'
        t1=[(set(self.get_taxa(n)),self.node(n).data.support) for n in self.all_ids() if \
            self.node(n).succ and\
            (self.node(n).data and self.node(n).data.support and self.node(n).data.support>=threshold)]
        t2=[(set(tree2.get_taxa(n)),tree2.node(n).data.support) for n in tree2.all_ids() if \
            tree2.node(n).succ and\
            (tree2.node(n).data and tree2.node(n).data.support and tree2.node(n).data.support>=threshold)]
        conflict=[]
        for (st1,sup1) in t1:
            for (st2,sup2) in t2:
                if not st1.issubset(st2) and not st2.issubset(st1):                     # don't hiccup on upstream nodes
                    intersect,notin1,notin2=st1 & st2, st2-st1, st1-st2                 # all three are non-empty sets
                    # if notin1==missing1 or notin2==missing2  <==> st1.issubset(st2) or st2.issubset(st1) ???
                    if intersect and not (notin1.issubset(missing1) or notin2.issubset(missing2)):         # omit conflicts due to missing taxa
                        conflict.append((st1,sup1,st2,sup2,intersect,notin1,notin2))
        return conflict
        
    def common_ancestor(self,node1,node2):
        """Return the common ancestor that connects to nodes.
        
        node_id = common_ancestor(self,node1,node2)
        """
        
        l1=[self.root]+self.trace(self.root,node1)
        l2=[self.root]+self.trace(self.root,node2)
        return [n for n in l1 if n in l2][-1]


    def distance(self,node1,node2):
        """Add and return the sum of the branchlengths between two nodes.
        dist = distance(self,node1,node2)
        """
        
        ca=self.common_ancestor(node1,node2)
        return self.sum_branchlength(ca,node1)+self.sum_branchlength(ca,node2)

    def is_monophyletic(self,taxon_list):
        """Return node_id of common ancestor if taxon_list is monophyletic, -1 otherwise.
        
        result = is_monophyletic(self,taxon_list)
        """
        if isinstance(taxon_list,str):
            taxon_set=set([taxon_list])
        else:
            taxon_set=set(taxon_list)
        node_id=self.root
        while 1:
            subclade_taxa=set(self.get_taxa(node_id))
            if subclade_taxa==taxon_set:                                        # are we there?
                return node_id
            else:                                                               # check subnodes
                for subnode in self.chain[node_id].succ:
                    if set(self.get_taxa(subnode)).issuperset(taxon_set):  # taxon_set is downstream
                        node_id=subnode
                        break   # out of for loop
                else:
                    return -1   # taxon set was not with successors, for loop exhausted

    def is_bifurcating(self,node=None):
        """Return True if tree downstream of node is strictly bifurcating."""
        if not node:
            node=self.root
        if node==self.root and len(self.node(node).succ)==3: #root can be trifurcating, because it has no ancestor
            return self.is_bifurcating(self.node(node).succ[0]) and \
                    self.is_bifurcating(self.node(node).succ[1]) and \
                    self.is_bifurcating(self.node(node).succ[2])
        if len(self.node(node).succ)==2:
            return self.is_bifurcating(self.node(node).succ[0]) and self.is_bifurcating(self.node(node).succ[1])
        elif len(self.node(node).succ)==0:
            return True
        else:
            return False



    def branchlength2support(self):
        """Move values stored in data.branchlength to data.support, and set branchlength to 0.0

        This is necessary when support has been stored as branchlength (e.g. paup), and has thus
        been read in as branchlength. 
        """

        for n in self.chain.keys():
            self.node(n).data.support=self.node(n).data.branchlength
            self.node(n).data.branchlength=0.0

    def convert_absolute_support(self,nrep):
        """Convert absolute support (clade-count) to rel. frequencies.
        
        Some software (e.g. PHYLIP consense) just calculate how often clades appear, instead of
        calculating relative frequencies."""

        for n in self._walk():
            if self.node(n).data.support:
                self.node(n).data.support/=float(nrep)

    def randomize(self,ntax=None,taxon_list=None,branchlength=1.0,branchlength_sd=None,bifurcate=True):
        """Generates a random tree with ntax taxa and/or taxa from taxlabels.
    
        new_tree = randomize(self,ntax=None,taxon_list=None,branchlength=1.0,branchlength_sd=None,bifurcate=True)
        Trees are bifurcating by default. (Polytomies not yet supported).
        """

        if not ntax and taxon_list:
            ntax=len(taxon_list)
        elif not taxon_list and ntax:
            taxon_list=['taxon'+str(i+1) for i in range(ntax)]
        elif not ntax and not taxon_list:
            raise TreeError('Either numer of taxa or list of taxa must be specified.')
        elif ntax<>len(taxon_list):
            raise TreeError('Length of taxon list must correspond to ntax.')
        # initiate self with empty root
        self.__init__()
        terminals=self.get_terminals()
        # bifurcate randomly at terminal nodes until ntax is reached
        while len(terminals)<ntax:
            newsplit=random.choice(terminals)
            new_terminals=self.split(parent_id=newsplit,branchlength=branchlength)
            # if desired, give some variation to the branch length
            if branchlength_sd:
                for nt in new_terminals:
                    bl=random.gauss(branchlength,branchlength_sd)
                    if bl<0:
                        bl=0
                    self.node(nt).data.branchlength=bl
            terminals.extend(new_terminals)
            terminals.remove(newsplit)
        # distribute taxon labels randomly
        random.shuffle(taxon_list)
        for (node,name) in zip(terminals,taxon_list):
            self.node(node).data.taxon=name

    def display(self):
        """Quick and dirty lists of all nodes."""
        table=[('#','taxon','prev','succ','brlen','blen (sum)','support')]
        for i in self.all_ids():
            n=self.node(i)
            if not n.data:
                table.append((str(i),'-',str(n.prev),str(n.succ),'-','-','-'))
            else:
                tx=n.data.taxon
                if not tx:
                    tx='-'
                blength=n.data.branchlength
                if blength is None:
                    blength='-'
                    sum_blength='-'
                else:
                    sum_blength=self.sum_branchlength(node=i)
                support=n.data.support
                if support is None:
                    support='-'
                table.append((str(i),tx,str(n.prev),str(n.succ),blength,sum_blength,support))
        print '\n'.join(['%3s %32s %15s %15s %8s %10s %8s' % l for l in table])
        print '\nRoot: ',self.root

    def to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True,plain_newick=False):
        """Return a paup compatible tree line.
       
        to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True)
        """
        # if there's a conflict in the arguments, we override plain=True
        if support_as_branchlengths or branchlengths_only:
            plain=False
        self.support_as_branchlengths=support_as_branchlengths
        self.branchlengths_only=branchlengths_only
        self.plain=plain

        def make_info_string(data,terminal=False):
            """Creates nicely formatted support/branchlengths."""
            # CHECK FORMATTING
            if self.plain: # plain tree only. That's easy.
                return ''
            elif self.support_as_branchlengths: # support as branchlengths (eg. PAUP), ignore actual branchlengths
                if terminal:    # terminal branches have 100% support
                    return ':%1.2f' % self.max_support
                else:
                    return ':%1.2f' % (data.support)
            elif self.branchlengths_only: # write only branchlengths, ignore support
                return ':%1.5f' % (data.branchlength)
            else:   # write suport and branchlengths (e.g. .con tree of mrbayes)
                if terminal:
                    return ':%1.5f' % (data.branchlength)
                else:
                    if data.branchlength is not None and data.support is not None:  # we have blen and suppport
                        return '%1.2f:%1.5f' % (data.support,data.branchlength)
                    elif data.branchlength is not None:                             # we have only blen
                        return '0.00000:%1.5f' % (data.branchlength)
                    elif data.support is not None:                                  # we have only support
                        return '%1.2f:0.00000' % (data.support)
                    else:
                        return '0.00:0.00000'

        def newickize(node):
            """Convert a node tree to a newick tree recursively."""

            if not self.node(node).succ:    #terminal
                return self.node(node).data.taxon+make_info_string(self.node(node).data,terminal=True)
            else:
                return '(%s)%s' % (','.join(map(newickize,self.node(node).succ)),make_info_string(self.node(node).data))
            return subtree
                    
        treeline=['tree']
        if self.name:
            treeline.append(self.name)
        else:
            treeline.append('a_tree')
        treeline.append('=')
        if self.weight<>1:
            treeline.append('[&W%s]' % str(round(float(self.weight),3)))
        if self.rooted:
            treeline.append('[&R]')
        treeline.append('(%s);' % ','.join(map(newickize,self.node(self.root).succ)))
        if plain_newick:
            return treeline[-2]
        else:
            return ' '.join(treeline)
        
    def __str__(self):
        """Short version of to_string(), gives plain tree"""
        return self.toNewick()
    #return self.to_string(plain=True)

    def toNewick(self, nodeId = None, topologyOnly = False, attributes = None, includeStem=False) :
      """ BioPython tree or sub-tree to unique NEWICK format.

      Child nodes are sorted (via text), so representation is unique and does not
      depend on arbitrary children ordering.
      """

      if nodeId is None:
        nodeId = self.root
      node = self.node(nodeId)
      data = node.data
      if not node.succ :
        rep = data.taxon if data.taxon is not None else ""
      else :
        reps = [self.toNewick(n, topologyOnly, attributes, True) for n in node.succ]
        reps.sort()
        rep = "(" + ",".join(reps) + ")"

      if attributes is not None :
        attrs = getattr(data, attributes, None)
        if attrs is not None and len(attrs) :
          s = '[&'
          for a in attrs:
            vl = str(attrs[a])
            if vl[0] != '{' and ',' in vl :
                vl = '{' + vl + '}'
            s += (a + '=' + vl + ",")
          s = s[:-1] + ']'
          rep += s

      if not topologyOnly and includeStem and data.branchlength is not None:
        rep = rep + (":%r" % data.branchlength)
      return rep
    
    def unroot(self):
        """Defines a unrooted Tree structure, using data of a rooted Tree."""

        # travel down the rooted tree structure and save all branches and the nodes they connect

        def _get_branches(node):
            branches=[]
            for b in self.node(node).succ:
                branches.append([node,b,self.node(b).data.branchlength,self.node(b).data.support])
                branches.extend(_get_branches(b))
            return branches
    
        self.unrooted=_get_branches(self.root)
        # if root is bifurcating, then it is eliminated
        if len(self.node(self.root).succ)==2:
            # find the two branches that connect to root
            rootbranches=[b for b in self.unrooted if self.root in b[:2]]
            b1=self.unrooted.pop(self.unrooted.index(rootbranches[0]))
            b2=self.unrooted.pop(self.unrooted.index(rootbranches[1]))
            # Connect them two each other. If both have support, it should be identical (or one set to None?).
            # If both have branchlengths, they will be added
            newbranch=[b1[1],b2[1],b1[2]+b2[2]]
            if b1[3] is None:
                newbranch.append(b2[3]) # either None (both rootbranches are unsupported) or some support
            elif b2[3] is None:
                newbranch.append(b1[3]) # dito
            elif b1[3]==b2[3]:          
                newbranch.append(b1[3]) # identical support
            elif b1[3]==0 or b2[3]==0:
                newbranch.append(b1[3]+b2[3]) # one is 0, take the other
            else:
                raise TreeError, 'Support mismatch in bifurcating root: %f, %f' % (float(b1[3]),float(b2[3]))
            self.unrooted.append(newbranch)

    def root_with_outgroup(self,outgroup=None):
        
        def _connect_subtree(parent,child):
            """Hook subtree starting with node child to parent."""
            for i,branch in enumerate(self.unrooted):
                if parent in branch[:2] and child in branch[:2]:
                    branch=self.unrooted.pop(i)
                    break 
            else:
                raise TreeError, 'Unable to connect nodes for rooting: nodes %d and %d are not connected' % (parent,child)
            self.link(parent,child)
            self.node(child).data.branchlength=branch[2]
            self.node(child).data.support=branch[3]
            #now check if there are more branches connected to the child, and if so, connect them
            child_branches=[b for b in self.unrooted if child in b[:2]]
            for b in child_branches:
                if child==b[0]:
                    succ=b[1]
                else:
                    succ=b[0]
                _connect_subtree(child,succ) 
            
        # check the outgroup we're supposed to root with
        if outgroup is None:
            return self.root
        outgroup_node=self.is_monophyletic(outgroup)
        if outgroup_node==-1:
            return -1
        # if tree is already rooted with outgroup on a bifurcating root,
        # or the outgroup includes all taxa on the tree, then we're fine
        if (len(self.node(self.root).succ)==2 and outgroup_node in self.node(self.root).succ) or outgroup_node==self.root:
            return self.root
        
        self.unroot()
        # now we find the branch that connects outgroup and ingroup
        #print self.node(outgroup_node).prev
        for i,b in enumerate(self.unrooted):
            if outgroup_node in b[:2] and self.node(outgroup_node).prev in b[:2]:
                root_branch=self.unrooted.pop(i)
                break
        else:
            raise TreeError, 'Unrooted and rooted Tree do not match'
        if outgroup_node==root_branch[1]:
            ingroup_node=root_branch[0]
        else:
            ingroup_node=root_branch[1]
        # now we destroy the old tree structure, but keep node data. Nodes will be reconnected according to new outgroup
        for n in self.all_ids():
            self.node(n).prev=None
            self.node(n).succ=[]
        # now we just add both subtrees (outgroup and ingroup) branch for branch
        root=Nodes.Node(data=NodeData())            # new root    
        self.add(root)                              # add to tree description
        self.root=root.id                           # set as root
        self.unrooted.append([root.id,ingroup_node,root_branch[2],root_branch[3]])  # add branch to ingroup to unrooted tree
        self.unrooted.append([root.id,outgroup_node,0.0,0.0])   # add branch to outgroup to unrooted tree
        _connect_subtree(root.id,ingroup_node)      # add ingroup
        _connect_subtree(root.id,outgroup_node)     # add outgroup
        # if theres still a lonely node in self.chain, then it's the old root, and we delete it
        oldroot=[i for i in self.all_ids() if self.node(i).prev is None and i!=self.root]
        if len(oldroot)>1:
            raise TreeError, 'Isolated nodes in tree description: %s' % ','.join(oldroot)
        elif len(oldroot)==1:
            self.kill(oldroot[0])
        return self.root
        
         
def consensus(trees, threshold=0.5,outgroup=None):
    """Compute a majority rule consensus tree of all clades with relative frequency>=threshold from a list of trees."""
    
    total=len(trees)
    if total==0:
        return None
    # shouldn't we make sure that it's NodeData or subclass??
    dataclass=trees[0].dataclass
    max_support=trees[0].max_support
    clades={}
    #countclades={}
    alltaxa=set(trees[0].get_taxa())
    # calculate calde frequencies
    c=0
    for t in trees:
        c+=1
        #if c%50==0:
        #    print c
        if alltaxa!=set(t.get_taxa()):
            raise TreeError, 'Trees for consensus must contain the same taxa'
        t.root_with_outgroup(outgroup=outgroup)
        for st_node in t._walk(t.root):
            subclade_taxa=t.get_taxa(st_node)
            subclade_taxa.sort()
            subclade_taxa=str(subclade_taxa) # lists are not hashable
            if subclade_taxa in clades:
                clades[subclade_taxa]+=float(t.weight)/total
            else:
                clades[subclade_taxa]=float(t.weight)/total
            #if subclade_taxa in countclades:
            #    countclades[subclade_taxa]+=t.weight
            #else:
            #    countclades[subclade_taxa]=t.weight
    # weed out clades below threshold
    for (c,p) in clades.items():
        if p<threshold:
            del clades[c]
    # create a tree with a root node
    consensus=Tree(name='consensus_%2.1f' % float(threshold),data=dataclass)
    # each clade needs a node in the new tree, add them as isolated nodes
    for (c,s) in clades.items():
        node=Nodes.Node(data=dataclass())
        node.data.support=s
        node.data.taxon=set(eval(c))
        consensus.add(node)
    # set root node data
    consensus.node(consensus.root).data.support=None
    consensus.node(consensus.root).data.taxon=alltaxa
    # we sort the nodes by no. of taxa in the clade, so root will be the last
    consensus_ids=consensus.all_ids()
    consensus_ids.sort(lambda x,y:len(consensus.node(x).data.taxon)-len(consensus.node(y).data.taxon))
    # now we just have to hook each node to the next smallest node that includes all taxa of the current 
    for i,current in enumerate(consensus_ids[:-1]): # skip the last one which is the root
        #print '----'
        #print 'current: ',consensus.node(current).data.taxon
        # search remaining nodes
        for parent in consensus_ids[i+1:]:
            #print 'parent: ',consensus.node(parent).data.taxon
            if consensus.node(parent).data.taxon.issuperset(consensus.node(current).data.taxon):
                break
        else:
            sys.exit('corrupt tree structure?')
        # internal nodes don't have taxa
        if len(consensus.node(current).data.taxon)==1:
            consensus.node(current).data.taxon=consensus.node(current).data.taxon.pop()
            # reset the support for terminal nodes to maximum
            #consensus.node(current).data.support=max_support
        else:
            consensus.node(current).data.taxon=None
        consensus.link(parent,current)
    # eliminate root taxon name
    consensus.node(consensus_ids[-1]).data.taxon=None 
    return consensus



def _getStuff(s, sep) :
  e = 0
  while s[e] != sep or  s[e-1] == '\\' :
    e += 1
  return e

def _parseAttributes(s) :
  if s[0] == '&' :
    vals = dict()
    s = s[1:]
    while len(s) :
        if s[0] == ',' :
            s = s[1:]

        nameEnd = s.find('=')
        name = s[:nameEnd]
        s = s[nameEnd+1:]
        if s[0] == '"' :
            e = _getStuff(s[1:],'"')
            v = s[1:e+1]
            s = s[e+2:]
        elif s[0] == '{' :
            e = _getStuff(s[1:], '}')
            v = s[1:e+1]
            s = s[e+2:]
        else :
            e = s.find(',')
            if e == -1 :
                e = len(s)
            v = s[:e]
            s = s[e:]

        vals[name] = v
    return vals
  return None
