## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

"""
==============
Tree Utilities
==============

Build and log trees, get clades and node heights, convert to NEWICK format and
so on. 

Unless explicitly specified, any tree is assumed to be Ultrametric
(tips are contemporaneous).

"""

from __future__ import division

import operator, sys, os.path
import random

from genericutils import fileFromName

# Bio.Nexus.Tree stuff

from ITrees import Tree, NodeData
import Nodes 

__all__ = ["TreeBuilder", "TreeLogger", "getClade", "getTreeClades",
           "getCommonAncesstor", "countNexusTrees", "toNewick", "nodeHeights",
           "nodeHeight", "treeHeight", "setLabels", "convertDemographics",
           "coalLogLike", "getPostOrder", "getPreOrder", "setSpeciesSimple",
           "resolveTree", "attributesVarName", "addAttributes",
           "rootAtMidpoint", "rootByTipVarianceOptimization", "treeDiameterInfo",
           "CAhelper"]

# Can't change, still hardwired in many places in code
attributesVarName = "attributes"

class TreeBuilder(object) :
  """ A basic helper for building BioPython trees.

  Use:
   - tb = TreeBuilder()
   - Create some terminals via leaf = tb.createLeaf(name)
   - use mergeNodes to successively merge nodes.
   - Call tb.finalize(root-node) to get the tree
  """
  
  def __init__(self, weight=1.0, rooted = True, name='') :
    self.t = Tree(weight=weight, rooted = rooted, name=name)

  def createLeaf(self, name) :
    """ Create a new leaf.

    :param name: taxon name
    :returns: node"""
    
    nd = NodeData()
    nd.taxon = name
    leaf = Nodes.Node(nd)
    self.t.add(leaf, None)
    return leaf

  def newNode(self) :
    nd = NodeData()
    node = Nodes.Node(nd)
    self.t.add(node, None)
    return node
  
  def mergeNodes(self, subtrees) :
    """ Join subtrees in subtrees into one sub-tree.

    Each element in subtrees is a (node,branch) pair, where branch is the length
    of the branch connecting the node to its parent. 

    :param subtrees: sequence of [node,branch]
    :returns: node"""
    
    nd = NodeData()
    node = Nodes.Node(nd)
    self.t.add(node, None)

    for n1,h1 in subtrees:
      n1.set_prev(node.id)
      n1.data.branchlength = h1

    node.add_succ([x.id for x,h in subtrees])
    return node

  def finalize(self, rootNode) :
    """ Convert a node to a proper tree with root rootNode.
    
    :param rootNode: node
    :returns: biopython tree
    """
    t = self.t
    rr = t.node(t.root)
    if rootNode.succ :
      rr.set_succ(rootNode.succ)
      for p in rootNode.succ :
        t.node(p).set_prev(t.root)
      if hasattr(rootNode.data,attributesVarName) :
        rr.data.attributes = rootNode.data.attributes
      t.kill(rootNode.id)
    else :
      rr.set_succ([rootNode.id])
      rootNode.set_prev(t.root)
    return t

def cleanLabel(l) :
  l = str(l)
  if l.isalnum() :
    return l
  if l[0] == l[-1] == "'":
    return l
  if l[0] == l[-1] == '"':
    return l
  
  if "'" in l :
    return '"' + l + '"'
  return "'" + l + "'"

class TreeLogger(object) :
  def __init__(self, outName = None, argv = None,
               version = None, overwrite = False,
               labels = None) :
    self.outName = outName
    if outName is not None:
      if not overwrite and os.path.isfile(outName) \
             and os.path.getsize(outName) > 0 :
        raise RuntimeError("not overwriting existing file " +  outName)
      
      self.outFile = outFile = file(outName, "w")
      self.count = 0
      if argv is not None :
        print >> outFile, "[Generated by %s%s]" % \
              (os.path.basename(argv[0]),
               (", version " + version) if version is not None else "")
        print >> outFile, "[%s]" % " ".join(argv)
        
      print >> outFile, "#NEXUS"
      print >> outFile, "begin trees;"
      if labels is not None:
        print >> outFile,"\tTranslate"
        for l,n in labels[:-1]:
          print >> outFile, "\t\t%s %s," % (str(l),cleanLabel(n))
        print >> outFile, "\t\t%s %s\n;" % (str(labels[-1][0]),cleanLabel(labels[-1][1]))
      
  def outTree(self, tree, treeAttributes = None, name = None) :
    c = ""

    # rooted attribute on by default for nexus files. If treeAttributes set by
    # called, caller should probably set it as well.
    if treeAttributes is None and self.outName :
      treeAttributes = {'R' : None}
      
    if treeAttributes is not None :
      c = " ".join(["[&%s %s]" % (k,v) if v else "[&%s]" % k
                    for k,v in treeAttributes.items()]) + " "
    if self.outName :
      if not name :
        name = "tree_%d" % self.count
      print >> self.outFile, "tree %s = %s%s ;" % (name,c,tree)
      self.count += 1
    else :
      print "%s%s" % (c,tree)
      
  def close(self) :
    if self.outName :
      print >> self.outFile, "end;"
      self.outFile.close()

# don't clobber exsiting attributes - append is much nicer
def addAttributes(nd, atrs) :
  if not hasattr(nd.data, attributesVarName) :
    nd.data.attributes = dict()
  nd.data.attributes.update(atrs)

def getClade(tree, nodeId) :
  n = tree.node(nodeId)
  if n.data.taxon :
    return [nodeId,]
  return reduce(operator.add, [getClade(tree, x) for x in n.succ])

#def getCommonAncesstor(tree, taxaIds) :
#  return reduce(tree.common_ancestor, taxaIds)

def getCommonAncesstor(tree, taxaIds) :
  ## faster for large subset of taxaIds
  if len(taxaIds) == 1 :
    return taxaIds[0]
  
  i = tree.common_ancestor(*taxaIds[:2])
  if len(taxaIds) > 2 :
    t = set(taxaIds[2:])
    t.difference_update(tree.get_taxa(i, asIDs=True))
    while len(t) :
      i = tree.common_ancestor(i, t.pop())
      t.difference_update(tree.get_taxa(i, asIDs=True))
  return i


def _getTreeClades_i(tree, nodeID) :
  node = tree.node(nodeID)
  if node.data.taxon :
    return ([node.data.taxon], [])

  cl = [_getTreeClades_i(tree, ch) for ch in node.succ]
  allt = reduce(operator.add, [x[0] for x in cl])
  clades = [(allt,node)]
  for cl1 in [x[1] for x in cl if x[1]] :
    clades.extend(cl1)
  return (allt, clades)

def getTreeClades(tree, trivialClades = False):
  """ Clades of subtree as a list of (taxa-list, tree-node).

  taxa-list is a list of strings.
  """

  c = _getTreeClades_i(tree, tree.root)[1]

  if trivialClades :
    for n in tree.get_terminals():
      nd = tree.node(n)
      c.append(([nd.data.taxon], nd))
    
  return c

def getPostOrder(tree, nodeId = None) :
  if nodeId is None:
    nodeId = tree.root
  node = tree.node(nodeId)
  p = [node]
  if node.succ :
    p = reduce(lambda x,y : x+y, [getPostOrder(tree, x) for x in node.succ] + [p])
  return p

def getPreOrder(tree, nid = None, includeTaxa = True) :
  if nid is None:
    nid = tree.root
  node = tree.node(nid)
  isLeaf = len(node.succ) == 0
  
  r = [nid]
  if isLeaf:
    if not includeTaxa :
      r = []
  else :
    r.extend(reduce(lambda x,y : x+y, [getPreOrder(tree, n, includeTaxa)
                                       for n in node.succ]))
  return r

def countNexusTrees(nexFileName) :
  """ Number of trees in a nexus file."""
  nexFile = fileFromName(nexFileName)
  c = 0
  for l in nexFile:
    if l.startswith("tree ") :
      c += 1
  return c

def toNewick(tree, nodeId = None, topologyOnly = False, attributes = None) :
  """ BioPython tree or sub-tree to unique NEWICK format.

  Child nodes are sorted (via text), so representation is unique and does not
  depend on arbitrary children ordering.
  """
  return tree.toNewick(nodeId, topologyOnly, attributes)


def _getNodeHeight(tree, n, heights, w = 0):
  if n.id in heights:
    return heights[n.id]
  
  if len(n.succ) == 0 :
    heights[n.id] = 0.0
    return 0.0

  i = n.succ[w]
  if i not in heights:
    if n.succ[1-w] in heights:
      w = 1-w
      i = n.succ[w]
      
  c = tree.node(i)
  h = _getNodeHeight(tree, c, heights, 1-w) + c.data.branchlength
  
  heights[n.id] = h
  return h

def _collectNodeHeights(tree, nid, heights) :
  node = tree.node(nid)
  
  if len(node.succ) == 0 :
    heights[node.id] = 0.0
    return (node.data.branchlength, [node.id])

  (hs0,tips0), (hs1,tips1) = [_collectNodeHeights(tree, c, heights)
                               for c in node.succ]
  if hs0 != hs1 :
    if hs0 < hs1 :
      h = hs1
      dx = hs1-hs0
      # add hs1-hs0 to all left (0) side
      for n in tips0 :
        # protect from numerical instability: insure parent is higher than son
        heights[n] = min(heights[n] + dx, h)
        assert heights[n] <= h
    else :
      # add hs0-hs1 to all right (1) side 
      h = hs0
      dx = hs0-hs1
      for n in tips1 :
        heights[n] = min(heights[n] + dx, h)
        assert heights[n] <= h
  else :
    h = hs0
  heights[node.id] = h
  return (h + node.data.branchlength, [nid] + tips0 + tips1)

def _collectNodeHeightsSimple(tree, nid, heights) :
  node = tree.node(nid)
  
  if len(node.succ) == 0 :
    heights[node.id] = 0.0
    return node.data.branchlength

  h = max([_collectNodeHeightsSimple(tree, c, heights) for c in node.succ])
  heights[node.id] = h
  return h + node.data.branchlength

def nodeHeights(tree, nids = None, allTipsZero = True) :
  """ Return a mapping from node ids to node heights.
  Without nids - for all nodes.

  The mapping may contain heights for other nodes as well.

  With !allTipsZero, handle non-ultrametric trees as well.
  """

  heights = dict()

  if allTipsZero :
    if nids == None :
      _collectNodeHeightsSimple(tree, tree.root, heights)
    else :
      w = 0
      for nid in nids :
        _getNodeHeight(tree, tree.node(nid), heights, w)
        w = 1-w
  else :
      # have to scan all tree to account for tip heights, ignore nids if
      # present
      _collectNodeHeights(tree, tree.root, heights)
     
  return heights

def nodeHeight(tree, nid) :
  """ Height of node. """
  
  node = tree.node(nid)
  if not node.succ :
    h = 0
  else :
    h = max([tree.node(c).data.branchlength + nodeHeight(tree, c)
             for c in node.succ])
    
  return h

def treeHeight(tree) :
  """ Height of tree. """
  return nodeHeight(tree, tree.root)
  # _getNodeHeight(tree, tree.node(tree.root), dict())


def setLabels(trees) :
  """ Convert labels attribute to tree node member (python list).

  Return boolean array per tree indicating if tree has complete meta-data.
"""
  hasAll = []
  for tree in trees :
    has = True 
    for i in tree.get_terminals() :
      data = tree.node(i).data
      ok = False
      if hasattr(data, attributesVarName) :
        if "labels" in data.attributes:
          l = data.attributes["labels"].split(' ')
          data.labels = l
          ok = len(l) > 0
      if not ok :
        has = False
    hasAll.append(has)

  return hasAll

def setSpeciesSimple(gtree, stree, sep = None) :
  """ Set a species for taxa in 'gtree'.

  if 'sep' is given, the gene taxon is split with respect to 'sep', where one
  and only one of the parts should match only one of the species names.
  Otherwise, the gene taxon should contain one and only one of the species
  names.
  
  data.snode is set to the taxon node in stree.
  """
  tx = stree.get_terminals()
  spNames = [stree.node(x).data.taxon for x in tx]
  for n in gtree.get_terminals() :
    sid = None
    gn = gtree.node(n)
    if sep is not None:
      s = gn.data.taxon.split(sep)[1]
      b = [x in s for x in spNames]
      if sum(b) == 1 :
        sid = b.index(True)
    else :
      b = [s in gn.data.taxon for s in spNames]
      if sum(b) == 1 :
        sid = b.index(True)
    # species tree node

    if sid is None :
      raise RuntimeError("problem with " +  gn.data.taxon)
    gn.data.snode = stree.node(tx[sid])

from demographic import LinearPiecewisePopulation, ConstantPopulation, StepFunctionPopulation

def _tod(xt, yt, b = None) :
  xt = [float(x) for x in xt.split(',') if len(x)]
  yt = [float(x) for x in yt.split(',')]
  if b is not None and len(xt) + 2 == len(yt) :
    xt.append(b)

  if len(yt) == 1 and len(xt) == 0 :
    return ConstantPopulation(yt[0])
  
  return LinearPiecewisePopulation(yt, xt)

def _toDemog(dtxt) :
  return _tod(*dtxt.split('|'))

def _toDemog1(dmt, dmv, branch) :
  return _tod(dmt if dmt is not None else "", dmv, b = branch)

def convertDemographics(tree, formatUnited = "dmf",
                        formatSeparated = ("dmv", "dmt"),
                        dattr = "demographic") :
  """ Convert demographic function stored in BEAST trees attributes to biopy
  demographic.

  Support old-style dmf attribute and new-style dmv,dmt
  """
  
  dmf = formatUnited
  dmv,dmt = formatSeparated

  missing = 0
  for i in tree.all_ids() :
    data = tree.node(i).data
    d = None
    if hasattr(data, attributesVarName) :
      if dmf in data.attributes:
        dtxt = data.attributes[dmf]
        d = _toDemog(dtxt)
      elif dmv in data.attributes:
        d = _toDemog1(data.attributes.get(dmt), data.attributes.get(dmv),
                      data.branchlength)
    if d is not None :
      setattr(data, dattr, d)
    else :
      missing += 1
  return missing


def _toAttrText(vals) :
  if len(vals) > 1 :
    return "{" + ",".join(["%f" % x for x in vals]) + "}"
  return "%f" % vals[0]
   
def revertDemographics(tree, dattr = "demographic",
                       formatSeparated = ("dmv", "dmt")) :
  dmv, dmt = formatSeparated
  
  for i in tree.all_ids() :
    data = tree.node(i).data

    d = getattr(data, dattr, None)
    if d :
      if not hasattr(data, attributesVarName) :
        data.attributes = dict()
      else :
        for l in formatSeparated :
          if l in data.attributes:
            del data.attributes[l]
            
      if isinstance(d, ConstantPopulation) :
        data.attributes[dmv] = d.population(0)
      elif isinstance(d, LinearPiecewisePopulation) :
        data.attributes[dmv] = _toAttrText(d.vals)
        if len(d.xvals) :
          data.attributes[dmt] = _toAttrText(d.xvals)
      elif isinstance(d, StepFunctionPopulation) :
        data.attributes[dmv] = _toAttrText(d.vals)
        data.attributes[dmt] = _toAttrText([0] + d.xvals)
      else :
        raise RuntimeError("unsupported")

import coalescent

def coalLogLike(tree, demog, condOnTree = False) :
  nh = nodeHeights(tree, allTipsZero = False)  
  terms = tree.get_terminals()
  
  times = sorted([(t,nid not in terms) for nid,t in nh.items()])

  return coalescent.coalLogLike(demog, times, condOnTree)

def resolveTree(tree) :
  for n in getPostOrder(tree):
    if len(n.succ) > 2 :
      cans = [tree.node(x) for x in n.succ]
      while len(cans) > 2 :
        i,j = random.sample(range(len(cans)), 2)
        if i > j :
          i,j = j,i

        nd = NodeData(branchlength=None)
        node = Nodes.Node(nd)
        node.set_succ([cans[i].id, cans[j].id])
        nid = tree.add(node, n.id)
        cans[i].set_prev(nid)
        cans[j].set_prev(nid)
        
        cans[i] = node
        del cans[j]
      n.set_succ([x.id for x in cans])

  return tree


def _maxTipPath(n) :
  d = n.data
  m = d.mtip
  return m[0] + d.branchlength, m[1]

def _maxTipPartner(tree, x) :
  """ Tip with the maximum path length to x on the tree.
  Uses pre-set mtip data. """
  #ver = False
  
  mx = (-1, None)
  h = 0
  n = tree.node(x)
  while n.id != tree.root :
    #if ver: print n.id, h, (mx[0],mx[1].data.taxon if mx[1] else None)
    p = n.prev
    h += n.data.branchlength
    np = tree.node(p)
    l = (-1,None,None)
    for c,m in zip(np.succ,np.data.mtips) :
      if c != n.id :
        if m[0] > l[0] :
          l = m
    l = l[0] + h, l[1]
    #if ver: print (l[0],l[1].data.taxon),
    if l[0] > mx[0] :
      mx = l + (np,)
    #if ver: print (mx[0],mx[1].data.taxon)
    n = np

  return mx
  
def treeDiameterInfo(tree) :
  """ Find the diameter of the tree. Return the diameter, the two tips at the
  ends and the internal common ancestor in the (rooted) 'tree'""" 
  po = getPostOrder(tree)
  for n in po:
    if not n.succ:
      n.data.mtip = (0, n)
      n.data.mtips = []
    else :
      n.data.mtips = [_maxTipPath(tree.node(c)) for c in n.succ]
      n.data.mtip = max(n.data.mtips)

  x = tree.get_terminals()[0]
  mx = _maxTipPartner(tree, x)
  o = _maxTipPartner(tree, mx[1].id)

  for n in po :
    del n.data.mtip, n.data.mtips
    
  return o[0], mx[1], o[1], o[2]

def _unrootedNewickRep(tree, n, c) :
  """ 'n' is internal. 'c' is a descendant of 'n'. Return NEWICK
  representation of the sub-tree whose root is 'c', with a single branch to 'n'.
  That is, the tree contains all non-c descendants of n, and the parent of 'n'
  "away" from 'n'. """

  if n.id == tree.root and len(n.succ) == 2 :
    nc = n.succ[0] if n.succ[0] != c.id else n.succ[1]
    snc = tree.toNewick(nc)
    b = tree.node(nc).data.branchlength + c.data.branchlength
    return snc, b
    
  p = []
  for nc in n.succ :
    if nc != c.id :
      snc = tree.toNewick(nc)
      p.append( "%s:%g" % (snc,tree.node(nc).data.branchlength) )
  if n.prev is not None:
    s,b = _unrootedNewickRep(tree, tree.node(n.prev), n)
    p.append("%s:%g" % (s,b))
  return '(' + ','.join(p) + ')', c.data.branchlength
      
def _rerootedNewickRep(tree, n, d) :
  """ Re-root inside the branch between 'n' and its parent, distance 'd' from
  'n'."""
  
  assert n.prev is not None and d <= n.data.branchlength

  s,b = _unrootedNewickRep(tree, tree.node(n.prev), n)
  lft = "%s:%g" % (s, b - d)
  s = tree.toNewick(n.id)
  rht = "%s:%g" % (s, d)
  return '(%s,%s)' % (lft,rht)
      
def rootAtMidpoint(tree) :
  """ Re-root unrooted 'tree' at "midpoint" - the middle of the maximum path
  between 2 tips.

  Return NEWICK text. Metadata not preserved.
  """
  
  diam,tip1,tip2,anode = treeDiameterInfo(tree)

  def xx(tip, d, a) :
    p = 0
    n = tip
    while n != a and p + n.data.branchlength < d :
      p += n.data.branchlength
      n = tree.node(n.prev)
    if n == a :
      return None
    return tip, n, d - p

  tip,n,e = xx(tip1, diam/2, anode) or xx(tip2, diam/2, anode)

  return _rerootedNewickRep(tree, n, e)

if 0 :
  def _populateTipDistancesFromParent(tree, n, parDists) :
    if n.id != tree.root :
      assert not n.data.dtips[-1] and parDists
      n.data.dtips[-1] = [ [a[0],a[1] + n.data.branchlength] for a in parDists]
      parDists = n.data.dtips[-1]
    else :
      assert n.data.dtips[-1] and not parDists
      parDists = []

    for i in range(len(n.succ)) :
      d = flatten([n.data.dtips[j] for j in range(len(n.succ)) if j != i] + [parDists])
      _populateTipDistancesFromParent(tree, tree.node(n.succ[i]), d)

  def _populateTreeWithNodeToTipDistances(tree) :
    for n in getPostOrder(tree) :
      if not n.succ:
        n.data.dtips = [[[n,0]],[],[]]
      else :
        ch = [tree.node(c) for c in n.succ]
        n.data.dtips = [[[a[0],a[1]+x.data.branchlength] for a in x.data.dtips[0]] +
                        [[a[0],a[1]+x.data.branchlength] for a in x.data.dtips[1]]
                        for x in ch]
        if n.id != tree.root :
           n.data.dtips.append([])

    _populateTipDistancesFromParent(tree, tree.node(tree.root), [])

import array
from itertools import imap, chain
#def farray(n, val = 0.0) :
#  return array.array('f',repeat(val,n))
if 0:
  def _populateTipDistancesFromParent(tree, n, parDists) :
    if n.id != tree.root :
      assert not n.data.dtips[-1] and parDists

      i = imap(lambda x : x + n.data.branchlength, parDists)
      parDists = n.data.dtips[-1] = array.array('f', i)
    else :
      assert n.data.dtips[-1] and not parDists
      parDists = []

    for i in range(len(n.succ)) :
      d = chain(*([n.data.dtips[j] for j in range(len(n.succ)) if j != i] + [parDists]))
      _populateTipDistancesFromParent(tree, tree.node(n.succ[i]), d)

  def _populateTreeWithNodeToTipDistances(tree) :
    for n in getPostOrder(tree) :
      if not n.succ:
        n.data.dtips = [array.array('f',[0]),[],[]]
      else :
        ch = [tree.node(c) for c in n.succ]
        n.data.dtips = [[a+x.data.branchlength for a in x.data.dtips[0]] +
                        [a+x.data.branchlength for a in x.data.dtips[1]]
                        for x in ch]
        if n.id != tree.root :
           n.data.dtips.append([])

    _populateTipDistancesFromParent(tree, tree.node(tree.root), [])

  def _cleanTreeWithNodeToTipDistances(tree) :
    for nid in tree.all_ids() :
      n = tree.node(nid)
      del n.data.dtips

if 0 :
  def _rootPointByTipVarianceOptimization(tree) :
    _populateTreeWithNodeToTipDistances(tree)
    minLoc = float('inf'),None,None

    for nid in tree.all_ids() :
      if nid == tree.root :
        continue
      n = tree.node(nid)
      ## pl,mn = [flatten([[a[1] for a in d] for d in n.data.dtips[:-1] if d])] + \
      ##         [[a[1] for a in n.data.dtips[-1]]]
      pl,mn = array.array('f',chain(*[d for d in n.data.dtips[:-1] if d])), n.data.dtips[-1]

      nl = len(mn)+len(pl)
      spl, smn = sum(pl), sum(mn)

      b,c = 2 * (spl - smn)/nl, sum([x**2 for x in pl + mn])/nl
      a1,b1 = (len(pl) - len(mn))/nl, (spl + smn)/nl

      ac,bc,cc = (1 - a1**2),  (b - (2 * a1 * b1)), (c - b1**2)

      dx = min(max(-bc / (2 * ac) , 0), n.data.branchlength)

      val = dx**2 * ac + dx * bc +  cc
      #print n.id,dx,val
      if val < minLoc[0] :
        minLoc = (val, n, dx)

    _cleanTreeWithNodeToTipDistances(tree)
    return minLoc

def _updateSums(s, x) :
  if not s :
    return []
  n = s[0]
  return [n, s[1] + x*n, s[2] + 2*s[1]*x + n*x**2]

def _addSums(sms) :
  return [sum(x) for x in zip(*sms)]

def _populateTipDistancesFromParentForm(tree, n, parDists) :
  if n.id != tree.root :
    assert not n.data.sums[-1][0] and parDists
    n.data.sums[-1] = _updateSums(parDists, n.data.branchlength)
    parDists = n.data.sums[-1]
  else :
    assert n.data.sums[-1] and not parDists
    parDists = [0,0,0]
    
  for i in range(len(n.succ)) :
    d = _addSums([n.data.sums[j] for j in range(len(n.succ)) if j != i] + [parDists])
    _populateTipDistancesFromParentForm(tree, tree.node(n.succ[i]), d)


def _populateTreeWithNodeToTipDistancesForm(tree) :
  for n in getPostOrder(tree) :
    if not n.succ:
      n.data.sums = [[1,0,0],[0,0,0],[0,0,0]]
    else :
      ch = [tree.node(c) for c in n.succ]
      n.data.sums = [_addSums([_updateSums(x.data.sums[i],x.data.branchlength)
                              for i in range(len(x.data.sums)-1)])
                     for x in ch]
      if n.id != tree.root :
         n.data.sums.append([0,0,0])
  _populateTipDistancesFromParentForm(tree, tree.node(tree.root), [])

def _cleanTreeWithNodeToTipDistancesForm(tree) :
  for nid in tree.all_ids() :
    n = tree.node(nid)
    del n.data.sums

def _rootPointByTipVarianceOptimization(tree) :
  _populateTreeWithNodeToTipDistancesForm(tree)
  minLoc = float('inf'),None,None
  
  for nid in tree.all_ids() :
    if nid == tree.root :
      continue
    n = tree.node(nid)

    npl,nmn = (sum([x[0] for x in n.data.sums[:-1]]) , n.data.sums[-1][0])
    nl = npl + nmn
    spl, smn = sum([x[1] for x in n.data.sums[:-1]]), n.data.sums[-1][1]
    
    b,c = 2 * (spl - smn) / nl, sum([x[2] for x in n.data.sums])/nl

    a1,b1 = (npl - nmn)/nl, (spl + smn)/nl

    ac,bc,cc = (1 - a1**2),  (b - 2 * a1 * b1), (c - b1**2)
    
    dx = min(max(-bc / (2 * ac) , 0), n.data.branchlength)
    
    val = dx**2 * ac + dx * bc +  cc
    #print n.id,dx,val
    if val < minLoc[0] :
      minLoc = (val, n, dx)

  _cleanTreeWithNodeToTipDistancesForm(tree)
  return minLoc

def rootByTipVarianceOptimization(tree) :
  minLoc = _rootPointByTipVarianceOptimization(tree)
  return _rerootedNewickRep(tree, minLoc[1], minLoc[2])

class CAhelper(object) :
  """ Augment tree node to allow fast search of common ancestor (use when
  performing many CA operations on the same (large) tree.

  Assumes all tip times are 0 (no checks).
  """
  
  def __init__(self, tr, th = None) :
    self.tree = tr
    self.th = th
    self.dterms = dict()

    nterms = len(tr.get_terminals())

    self.tree.size = nterms
    
    for n in getPostOrder(tr) :
      if not n.succ:
        n.data.tl = 0
        n.data.rh = 0
        n.data.level = 0
        n.data.path = []
        n.data.terms = [n]
        n.data.cladesize = 0
        self.dterms[n.data.taxon] = n
      elif not (n.id == tr.root and nterms==1) :
        ch = [tr.node(x).data for x in n.succ]
        n.data.rh = max([c.rh + c.branchlength for c in ch])
        for c in ch:
          c.branchlength = n.data.rh - c.rh
        n.data.tl = sum([c.tl + c.branchlength for c in ch])
        n.data.level = max([c.level for c in ch]) + 1
        n.data.terms = reduce(lambda x,y : x+y, [c.terms for c in ch])
        n.data.cladesize = len(n.data.terms)

        for x in n.data.terms:
          x.data.path.append(n.id)

      if th is not None:
        n.data.croot = n.data.rh <= th and n.data.rh + n.data.branchlength > th

    #del tr.node(tr.root).data.terms
    for t in tr.get_terminals() :
      d = tr.node(t).data
      d.pathset = set(d.path)

  def __call__(self, *txNodes) :
    return self.getCA(txNodes)

  def getCA(self, txNodes) :
    if len(txNodes) == 1 :
      return txNodes[0]
    v = reduce(set.intersection, [x.data.pathset for x in txNodes])
    n = [self.tree.node(x) for x in v]
    a = min(n, key = lambda x : x.data.level)
    return a

  def getCAi(self, txNodes) :
    return self.getCA([self.dterms[x] for x in txNodes])

  def taxonToNode(self, t) :
    return self.dterms[t]
  
  def clade(self, n) :
    if isinstance(n, (int,long)) :
      n = self.tree.node(n)
    return n.data.terms
