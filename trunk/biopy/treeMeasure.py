## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

""" Tree measures.

  Distance between trees and related utilities.
"""

from math import sqrt
import numpy, numpy.linalg

__all__ = ["branchScoreTreeDistance", "treeScoreDistance",
           "allPartitions", "treeLength", "treeLengthNormed",
           "treeArea", "vdistance"]

def _collectCladeTaxa(tree, nodeId, taxa, partitions) :
  """ Return a (reverse) mapping of taxa for each node in the sub-tree below
  C{nodeId}.

  The returned mapping has the taxa as key, in the form of a set of integeres,
  each an index into C{taxa}. The value of the mapping is the node tree id.
  """
  
  node = tree.node(nodeId)
  if node.data.taxon :
    p = set([taxa.index(node.data.taxon)])
  else :
    p = set()
    for ch in node.succ:
      p.update(_collectCladeTaxa(tree, ch, taxa, partitions))

  if nodeId == tree.root :
    return
  
  partitions[frozenset(p)] = nodeId
  return p

def vdistance(a1, a2, order = 2) :
  """ Vector norm M{||v||_2}"""
  return numpy.linalg.norm([x-y for x,y in zip(a1,a2)], order)

def branchScoreTreeDistance(tree1, tree2, distanceMetric = vdistance) :
  """ Branch score tree distance.

  Distance between two rooted trees based on a "minimal edit distance", the
  total sum of branch lengths changes in a sequence of moves which transform one
  tree to the other. (U{Kuhner
  1994<http://mbe.library.arizona.edu/data/1994/1103/13kuhn.pdf>})

  Individual branch changes are converted to a single distance using the vector
  norm C{distanceMetric}.
  """
  
  taxa = tree1.get_taxa()
  p1, p2 = dict(), dict()
  _collectCladeTaxa(tree1, tree1.root, taxa, p1)
  _collectCladeTaxa(tree2, tree2.root, taxa, p2)

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

  return distanceMetric(d1,d2)

def treeScoreDistance(tree1, tree2, norm = vdistance, consistencyCheck = True) :
  """ Tree distance when taking into account both branch length and population
  size on the branch.

  Similar to L{branchScoreTreeDistance}, but taking the intensity (integral of
  1/pop-size over the branch) instead of branch length. (with a constant
  population size of 1, it reduces to branch score). 

  Individual scores are converted to a single distance using the vector norm
  C{norm}.
  """
  
  taxa = tree1.get_taxa()
  p1, p2 = dict(), dict()
  _collectCladeTaxa(tree1, tree1.root, taxa, p1)
  _collectCladeTaxa(tree2, tree2.root, taxa, p2)

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

def allPartitions(referenceTree, trees) :
  """ Clade information summary for a set of trees.

  Summerize clades from all C{trees} in one mapping. All trees must be on the
  same taxa as the reference tree. Return a mapping whose key is the clade, and
  the value is a sequence of C{(tree,node)} pairs.
  """
  
  taxa = referenceTree.get_taxa()
  p = dict()
  for tree in trees: 
    p1 = dict()
    _collectCladeTaxa(tree, tree.root, taxa, p1)
    for k in p1 :
      if k not in p:
        p[k] = []
      p[k].append((tree,p1[k]))

  return p
