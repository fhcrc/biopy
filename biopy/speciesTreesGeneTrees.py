## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

""" Compute probabilities of gene trees under the multispecies coalescent.

Simulate gene trees based on the multispecies coalescent.
"""

from __future__ import division

from math import exp
from combinatorics import nPairs as c2
from treeCombinatorics import toNewick

from treeCombinatorics import nLabeledHistories, numberOfLabeledForests, \
     allCompatibleLabeledHistories

__all__ = ["compatibleGeneTreesInSpeciesTree", "simulateGeneTree"]

_c_vals = dict()

def c_t(n, k, i) :
  """ Tavare coefficients """
  assert k <= i <= n

  if (n,k,i) in _c_vals :
    return _c_vals[n,k,i]
  
  if k == i and (k == 1 or n == k) :
    v = 1
  else :
    if i > k :
      def f(i) :
        return (-(2*i+1) * (k+i-1) * (n-i)) / ( (2*i-1) *(i+1-k) * (n+i) )
      v = f(i-1) * c_t(n, k, i-1)
    else :
      assert n > k 
      def f(n) :
        return ((n+1)*n)/((n-k+1)*(n+k))
      v = f(n-1) * c_t(n-1, k, i)

  _c_vals[n,k,i] = v
  return v
    
def nToKlinages(n, k, d, t) :
  """ Probability of k remaining lineages from n after time t and population
  size function d"""
  
  return sum([c_t(n,k,i) * exp(-c2(i) * d.integrate(t)) for i in range(k, n+1)])

def stripLeaves(p) :
  if len(p) == 1 :
    return p[0]
  return [stripLeaves(x) for x in p]

def stripForest(f) :
  return [stripLeaves(t) for t in f]

def standarizeTree(p) :
  """ Put tree in standard form, so that trees can be compared using '=' """
  if len(p) == 1:
    return p
  return sorted([standarizeTree(x) for x in p])

def standarizeForest(f) :
  """ Put forest in standard form, so that forests can be compared using '=' """
  return sorted([standarizeTree(t) for t in f])

# With a compat function which always returns True, this computes probabilities
# for all gene trees.

def _compatibleGeneTreesInSpeciesTree(tree, nodeId, compat) :
  node = tree.node(nodeId)
  
  if node.data.taxon :
    labels = node.data.labels
    forests = [(1.0, list([[x] for x in labels])),]
  else :
    cth0,cth1 = [_compatibleGeneTreesInSpeciesTree(tree, child, compat)
                 for child in node.succ]
    # combine them
    forests = [(p0*p1,f0+f1) for p0,f0 in cth0 for p1,f1 in cth1]

  resultForests = []
  demo = node.data.demographic
  b = node.data.branchlength
  isRoot = nodeId == tree.root

  theoDict = dict()

  for p,forest in forests:
    n = len(forest)

    forestsInds = dict([[str(fx), kx] for kx,(px,fx) in enumerate(resultForests)])

    nlf = [nLabeledHistories(n,k) for k in range(1, n+1)]

    for lf in allCompatibleLabeledHistories([[x] for x in forest], compat):
      n1 = len(lf)
      if isRoot and n1 != 1:
        continue

      ratio = numberOfLabeledForests(lf)/nlf[n1-1]
      
      th = theoDict.get( (n, n1 ) )
      if th is None :
        th = nToKlinages(n, n1, demo, b) if not isRoot else 1

        theoDict[ (n, n1) ] = th
        
      lf1 = stripForest(lf)

      p1 = p * ratio * th
      sig = standarizeForest(lf1)
      ssig = str(sig)
      ni = forestsInds.get(ssig)
      if ni is None :
        forestsInds[ssig] = len(resultForests)
        resultForests.append([p1, sig])
      else :
        resultForests[ni][0] += p1
  
  return resultForests


def compatibleGeneTreesInSpeciesTree(tree, compat = None) :
  trees = _compatibleGeneTreesInSpeciesTree(tree, tree.root, compat)
  return [(t[0], toNewick(t[1][0])) for t in trees]


from treeutils import TreeBuilder, nodeHeights
from coalescent import getArrivalTimes
import random

def _simulateGeneTreeForNode(tree, nodeId, simTree, nodeHeights) :
  """ Simulate gene tree for sub-tree of species tree C{tree} rooted at father
  of C{nodeId}.

  Return a list of pairs (s,h) where s is a sub-tree in NEWICK format and h is
  height of s.
  """
  
  node = tree.node(nodeId)
  if node.data.taxon :
    tips = node.data.geneTreeTips
    strees = [list([simTree.createLeaf(x), 0.0]) for x in tips]
  else :
    strees = [_simulateGeneTreeForNode(tree, child, simTree, nodeHeights)
              for child in node.succ]

    # Don't care where they came from, all are in one ancestral species now.
    strees = strees[0] + strees[1]

    # if verbose > 2 : print nodeId,"strees",strees
    
  # Total lineages at start
  nl = len(strees)
  
  # Coalesce upwards until end of branch 
  demo = node.data.demographic
  aTimes = getArrivalTimes(demo, nl)
  
  # if verbose > 2 : print nodeId, "demo",demo, "at",aTimes
  # time 0 relative to sub-species
  
  notRoot = nodeId != tree.root

  branch = node.data.branchlength if notRoot else float('inf')
  height = nodeHeights[node.id]

  # if verbose > 2 : print "br",branch,"ht",height
  
  for t in aTimes :
    if t > branch or len(strees) == 1 :
      break
    
    i,j = random.sample(range(len(strees)), 2)
    l1,l2 = strees[i],strees[j]

    # height in tree axis
    th = height + t
    
    # if verbose> 2 : print "th",th

    subtree = simTree.mergeNodes(l1[0], th - l1[1], l2[0], th - l2[1])
  
    strees.pop(max(i,j))
    strees[min(i,j)] = [subtree, th]

  return strees

def simulateGeneTree(sTree) :
  """ Simulate one gene tree under species (spp) tree C{sTree}.

  Return a pair of tree and root height.
"""

  nh = nodeHeights(sTree)

  simTree = TreeBuilder()
  t,rootHeight = _simulateGeneTreeForNode(sTree, sTree.root, simTree, nh)[0]
  t = simTree.finalize(t)
  return (t,rootHeight)
