## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

from __future__ import division

from treeCombinatorics import nLabeledHistories, numberOfLabeledForests, \
     allCompatibleLabeledHistories

__all__ = ["compatibleGeneTreesInSpeciesTree"]

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

def compatibleGeneTreesInSpeciesTree(tree, nodeId, compat) :
  node = tree.node(nodeId)
  
  if node.data.taxon :
    labels = node.data.labels
    forests = [(1.0, list([list(x) for x in labels])),]
  else :
    cth0,cth1 = [compatibleGeneTreesInSpeciesTree(tree, child, compat)
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
