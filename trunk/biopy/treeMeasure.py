## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

"""
Tree measures
=============

Distance between trees and related utilities.
"""
from __future__ import division

from math import sqrt, log
import numpy, numpy.linalg
from collections import defaultdict
import operator

from treeutils import nodeHeights, getPostOrder

__all__ = [
  "branchScoreTreeDistance", "treeScoreDistance", "heightsScoreTreeDistance",
  "heightsScoreTreeDistanceAll",
  "allPartitions", "treeLength", "treeLengthNormed", "treeArea", "vdistance",
  "rootedAgreementScore", "conditionalCladeScore", "cladesInTreesSet",
  "setTreeClades"
  ]

def _collectCladeTaxa(tree, nodeId, taxa, partitions, withHeights) :
  """ Return a (reverse) mapping of taxa for each node in the sub-tree below
  nodeId.

  The returned mapping has the taxa as key, in the form of a set of integeres,
  each an index into taxa. The value of the mapping is the node tree id.

  If withHeights, the map value if a pair (node-id,node-height). Assumes all
  tips at 0.
  """
  
  node = tree.node(nodeId)
  if node.data.taxon :
    p = [taxa.index(node.data.taxon),]
    h = 0
  else :
    ch = [_collectCladeTaxa(tree, ch, taxa, partitions, withHeights)
          for ch in node.succ]
    p = []
    if withHeights:
      for c,h in ch:
        p.extend(c)
      # leaves h set
    else :
      for c in ch:
        p.extend(c)

  if nodeId == tree.root :
    if withHeights :
      return h
    return

  if withHeights :
    partitions[frozenset(p)] = (nodeId, h)
    return (p, h + node.data.branchlength) 
  else :
    partitions[frozenset(p)] = nodeId
    return p

def _collectCladeTaxaNonATZ(tree, taxa, partitions) :
  nhs = nodeHeights(tree, allTipsZero = False)
  for n in getPostOrder(tree):
    data = n.data
    if not n.succ:
      data.clade = [taxa.index(data.taxon),]
    else :
      p = []
      for s in n.succ :
        d = tree.node(s).data
        p.extend(d.clade)
        del d.clade
      data.clade = p
      
    partitions[frozenset(data.clade)] = (n.id, nhs[n.id])
  return nhs[tree.root]

def vdistance(a1, a2, order = 2) :
  """ Vector norm {||v||_2}"""
  return numpy.linalg.norm([x-y for x,y in zip(a1,a2)], order)

def branchScoreTreeDistance(tree1, tree2, distanceMetric = vdistance) :
  """ Branch score tree distance.

  Distance between two rooted trees based on a *minimal edit distance*, the
  total sum of branch lengths changes in a sequence of moves which transform one
  tree to the other. (`Kuhner
  1994 <http://mbe.library.arizona.edu/data/1994/1103/13kuhn.pdf>`_)

  Individual branch changes are converted to a single distance using the vector
  norm distanceMetric.
  """
  
  taxa = tree1.get_taxa()
  p1, p2 = dict(), dict()
  _collectCladeTaxa(tree1, tree1.root, taxa, p1, False)
  _collectCladeTaxa(tree2, tree2.root, taxa, p2, False)

  d1,d2 = [], []
  
  for p in p1 :
    b1 = tree1.node(p1[p]).data.branchlength
    if p in p2 :
      b2 = tree2.node(p2[p]).data.branchlength
      del p2[p]
    else :
      b2 = 0
    d1.append(b1) ; d2.append(b2)

  for n in p2.values() :
    d2.append(tree2.node(n).data.branchlength);
    d1.append(0)

  if isinstance(distanceMetric, (list,tuple)) :
    return [m(d1,d2) for m in distanceMetric]
  
  return distanceMetric(d1,d2)

def _hs_score(tree1,tree2,h1,h2,p1,p2) :
  sd = abs(h1 - h2)

  notin2 = dict.fromkeys(p2)
  for p in p1 :
    if p in p2 :
      sd += abs(p1[p][1] - p2[p][1])
      del notin2[p]
    else :
      sd += tree1.node(p1[p][0]).data.branchlength

  for c in notin2 :
    sd += tree2.node(p2[c][0]).data.branchlength
  ## for n,h in p2.itervalues() :
  ##   sd += tree2.node(n).data.branchlength

  return sd
  
def heightsScoreTreeDistance(tree1, tree2) :
  """ Hybrid Heights/Branch score tree distance.

  Distance between two rooted trees. Sum all branch lengths of clades present in
  only one tree, add diffrence in clade heights for shared clades.
  """
  
  taxa = tree1.get_taxa()
  p1, p2 = dict(), dict()
  h1 = _collectCladeTaxa(tree1, tree1.root, taxa, p1, True)
  h2 = _collectCladeTaxa(tree2, tree2.root, taxa, p2, True)

  if 0 :
    sd = _hs_score(tree1,tree2,h1,h2,p1,p2)
  else :
    sd = abs(h1 - h2)

    for p in p1 :
      if p in p2 :
        sd += abs(p1[p][1] - p2[p][1])
        del p2[p]
      else :
        sd += tree1.node(p1[p][0]).data.branchlength

    for n,h in p2.itervalues() :
      sd += tree2.node(n).data.branchlength

  return sd

def heightsScoreTreeDistanceAll(tree, trees1) :
  taxa = trees1[0].get_taxa()
  ptrees1 = [dict() for x in trees1]
  hs1 = [_collectCladeTaxa(t, t.root, taxa, d, True)
         for t,d in zip(trees1, ptrees1)]
  if tree is None :
    res = []
    for nk,(tree1,h1,p1) in enumerate(zip(trees1,hs1,ptrees1)) :
      res.append([0]*(nk+1) + [_hs_score(tree1,trees1[k],h1,hs1[k],p1,ptrees1[k])
                               for k in range(nk+1,len(trees1))])
    for k in range(len(trees1)) :
      for i in range(k) :
        res[k][i] = res[i][k]
    return res
  else:
    if isinstance(tree, (list, tuple)) :
      trees2 = tree
      ptrees2 = [dict() for x in trees2]
    else :
      trees2 = [tree]
      ptrees2 = [dict()]
    hs2 = [_collectCladeTaxa(t, t.root, taxa, d, True) for t,d in zip(trees2, ptrees2)]

    res = []
    for tree2,h2,p2 in zip(trees2,hs2,ptrees2) :
      res.append([_hs_score(tree1,tree2,h1,h2,p1,p2) for tree1,h1,p1 in
                  zip(trees1,hs1,ptrees1)])
    return res if len(res) > 1 else res[0]

def treeScoreDistance(tree1, tree2, norm = vdistance, consistencyCheck = True) :
  """ Tree distance when taking into account both branch length and population
  size on the branch.

  Similar to :py:func:`branchScoreTreeDistance`, but taking the intensity
  (integral of 1/pop-size over the branch) instead of branch length. With a
  constant population size of 1, it reduces to branch score.

  Individual scores are converted to a single distance using the vector norm
  norm.
  """
  
  taxa = tree1.get_taxa()
  p1, p2 = dict(), dict()
  _collectCladeTaxa(tree1, tree1.root, taxa, p1, False)
  _collectCladeTaxa(tree2, tree2.root, taxa, p2. False)

  d1,d2 = [], []
  
  for p in p1 :
    n1 = tree1.node(p1[p])
    demo1 = n1.data.demographic
    if consistencyCheck:
      assert demo1.naturalLimit() is None or \
             demo1.naturalLimit() == n1.data.branchlength,\
             (demo1, demo1.naturalLimit() , n1.data.branchlength)
    b1 = demo1.integrate(n1.data.branchlength)
    if p in p2 :
      n2 = tree2.node(p2[p])
      demo2 = n2.data.demographic
      if consistencyCheck:
        assert demo2.naturalLimit() is None or \
               demo2.naturalLimit() == n2.data.branchlength,\
               (demo2, demo2.naturalLimit() , n2.data.branchlength)
      b2 = demo2.integrate(n2.data.branchlength)
      del p2[p]
    else :
      b2 = 0
    d1.append(b1) ; d2.append(b2)

  for p,n in p2.items() :
    n2 = tree2.node(n)
    demo2 = n2.data.demographic
    d2.append(demo2.integrate(n2.data.branchlength))
    d1.append(0)

  return norm(d1,d2)

def _allBranches(tree) :
  return [tree.node(x).data.branchlength for x in tree.all_ids()]

def treeLength(tree) :
  """ Tree length (total sum of branch lengths)."""
  return sum(_allBranches(tree))

def treeLengthNormed(tree, norm = lambda x : numpy.linalg.norm(x, 2)) :
  """ Tree length under a norm (default is || ||_2, sqrt of sum of squares)."""
  
  return norm(_allBranches(tree))

def _nodeBranchPop(tree, nodeId) :
  data = tree.node(nodeId).data
  br = data.branchlength
  demo = data.demographic

  assert demo.naturalLimit() is None \
         or numpy.allclose(demo.naturalLimit(), br, 1e-11, 1e-11), \
         (demo, demo.naturalLimit(), br)
  
  return demo.integrate(br)
  
def treeArea(tree) :
  """ Tree area.

  Total sum of intensity over all branches.
  """
  return sum([_nodeBranchPop(tree, x) for x in tree.all_ids() if x != tree.root])

def allPartitions(referenceTree, trees, func = None, withHeights = False,
                  withRoot = False) :
  """ Clade information summary for a set of trees.

  Summerize clades from all trees in one mapping. All trees must be on the
  same taxa as the reference tree. Return a mapping whose key is the clade, and
  the value is a sequence of (tree,node) pairs. If func is given, the values are
  the results of applying func to the (tree,node) pair. The clade is a (frozen)
  set of integers, each integer is the index of the taxon in the reference tree.
  """

  if withRoot :
    rootHeights = []
  taxa = referenceTree.get_taxa()
  ultrametric = True
  if withHeights :
    nhs = nodeHeights(referenceTree, allTipsZero = False)
    nhsd = [(referenceTree.node(i).data.taxon, nhs[i])
            for i in nhs if not referenceTree.node(i).succ]
    ultrametric = all([h <= 1e-13 for t,h in nhsd])
    
  p = dict()
  for tree in trees:
    p1 = dict()
    if ultrametric :
      h = _collectCladeTaxa(tree, tree.root, taxa, p1, withHeights)
    else :
      h = _collectCladeTaxaNonATZ(tree, taxa, p1)
    for k,nd in p1.iteritems() :
      pk = p.get(k)
      v = (tree, nd)
      if func :
        v = func(*v)
      if pk is None :
        p[k] = [v]
      else :
        pk.append(v)
    if withRoot :
      rootHeights.append(h)
      
  return (p,rootHeights) if withRoot else p

def _parHeight(node) :
  br = node.data.branchlength
  return node.data.height + (0 if br is None else br)

def _setTreeClades_i(tree, nodeID) :
  # should be enhanced for tip dates
  node = tree.node(nodeID)
  if node.data.taxon :
    node.data.clade = [node.data.taxon]
    node.data.height = 0
    return [node]

  cl = [_setTreeClades_i(tree, ch) for ch in node.succ]
  allt = reduce(operator.add, [x[0].data.clade for x in cl])
  for x in cl :
    x[0].data.clade = frozenset(x[0].data.clade)
  node.data.clade = allt
  br = cl[0][0].data.branchlength
  node.data.height = max([_parHeight(n[0]) for n in cl])
  return [node] + reduce(operator.add, cl)

def setTreeClades(tree):
  """ Set clade attribute in each node data.

  clade = frozenset(taxa) and height = float
  """
  r = _setTreeClades_i(tree, tree.root)
  r[0].data.clade = frozenset(r[0].data.clade)
  return r

def cladesInTreesSet(trees, withPairs=False, tidyup=True, func=None, withTaxa = False) :
  newitem = (lambda : 0) if func is None else (lambda : [])
  clc = defaultdict(newitem)
  clc2 = defaultdict(newitem) if withPairs else None
  for xtree in trees:
    nodes = setTreeClades(xtree)
    
    for node in nodes:
      cladeSet = node.data.clade
      if not withTaxa and len(cladeSet) == 1 :
        continue
      
      if func :
        dat = func(node.data)
        clc[cladeSet].append(dat)
      else :
        clc[cladeSet] += 1
        
      if withPairs :
        c2 = (cladeSet,frozenset((frozenset(xtree.node(s).data.clade)
                                  for s in node.succ)))
        if func :
          clc2[c2].append(dat)
        else :
          clc2[c2] += 1
          
      if tidyup :
        del node.data.clade   
        del node.data.height   

  if func is None :
    clc = dict([(k,x) for k,x in clc.iteritems()])
    if withPairs :
      clc2 = dict([(k,x) for k,x in clc2.iteritems()])
      return (clc, clc2)
    return clc
  else :
    clc = dict([(k,x) for k,x in clc.iteritems()])
    if withPairs :
      clc2 = dict([(k,x) for k,x in clc2.iteritems()])
      return (clc, clc2)
    return clc

def conditionalCladeScore(tree, clc12, laplaceCorrection = False) :
  clc, clc2 = clc12

  tidyup = False
  if not hasattr(tree.node(tree.root).data, "clade") :
    setTreeClades(tree)
    tidyup = True
    
  totlog = 0
  for i in tree.all_ids():
    node = tree.node(i)
    if node.succ :
      c = node.data.clade
      if len(c) <= 2 :
        continue
      k = (c, frozenset([tree.node(s).data.clade for s in node.succ]))
      c2 = clc2.get(k)
      if c2 :
        if laplaceCorrection :
          totlog += log((c2+(1./(2.0**(len(c)-1)-1)))/(clc[c]+1))
        else :
          totlog += log(c2/clc[c])
      else :
        if laplaceCorrection :
          c1 = clc.get(c)
          totlog += -log((2.0**(len(c)-1)-1) * ((1 if c1 is None else c1+1)))
        else :
          totlog = float("-inf")
          break
  if tidyup :
    for i in tree.all_ids() :
      del tree.node(i).data.clade

  return totlog

# t1 = parseNewick('(((c,d),b),a)')
# t2 = parseNewick('(((a,c),b),d)')
# c1,c2 = biopy.treeMeasure.cladesInTreesSet([t1,t2], withPairs=1)
# 1 == sum([exp(biopy.treeMeasure.conditionalCladeScore(x, c1, c2))
#          for x in [parseNewick(toNewick(y))
#                     for y in allTrees([x for x in "abcd"])]])

def rootedAgreementScore(tree1, tree2, scaled=False) :
  t1nodes = setTreeClades(tree1)
  t2nodes = setTreeClades(tree2)
  t2map = dict([(n.data.clade,n) for n in t2nodes])

  tot,tl1,tl2 = 0.0,0,0
  for t1n in t1nodes :
    b1 = t1n.data.branchlength
    tl1 += b1
    t2n = t2map.get(t1n.data.clade)
    if t2n is None :
      #print "no",t1n.data.clade
      tot += b1
    else :
      b2 = t2n.data.branchlength
      l1,u1 = t1n.data.height, t1n.data.height + b1
      l2,u2 = t2n.data.height, t2n.data.height + b2
      # for an unknown reason, this produces 0 when tree is compared against
      # itself, while the simple form produces epsilon
      dd = min(b1 + b2 + 2*max(l1,l2) - 2*min(u1,u2), b1+b2)
      tot += dd
      #tot += b1 + b2 - 2*max(min(u1,u2) - max(l1,l2), 0)
      t2n.data.shared = 1
      #print l1,u1,l2,u2,dd, tot
      
  for t2n in t2nodes:
    b2 = t2n.data.branchlength
    tl2 += b2
    if hasattr(t2n.data, "shared") :
      del t2n.data.shared
    else :
      tot += b2
      #print tot
      
  return tot/(tl1+tl2) if scaled else tot
