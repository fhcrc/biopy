## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

"""
Tree posterior distance optimization
====================================

"""
__all__ = ["minPosteriorDistanceTree"]

import scipy, scipy.optimize, copy, random

from treeMeasure import allPartitions
from treeutils import treeHeight, nodeHeights

## def _treeBranchAssignmentExprs(tree, clades) :
##   allid = set(tree.all_ids())
##   terms = set(tree.get_terminals())
##   allint = allid - terms
##   # works for dated tips as well - but rest of code does not support that
##   nh = nodeHeights(tree, allTipsZero = False)
##   nhInOrder = sorted([(nh[x],x) for x in allint])
##   nInOrder = [x[1] for x in reversed(nhInOrder)]

##   # x[0] is root
  
##   sr = []
##   sr.append("h%d = x[0]" % nInOrder[0])
##   for k in nInOrder[1:]:
##     sr.append("h%d = x[%d] * h%d" % (k, nInOrder.index(k), tree.node(k).prev))
    
##   for r,k in enumerate(clades) :
##     n = clades[k][0][1] ; assert n != tree.root
##     if n in terms:
##       sr.append("b%d=(h%d) # %d" % (r, tree.node(n).prev, n))
##     else :
##       sr.append("b%d=(h%d * (1-x[%d])) # %d" %
##               (r, tree.node(n).prev, nInOrder.index(n), n))
##   return sr


def _treeBranchAssignmentExprs(tree, clades) :
  allid = set(tree.all_ids())
  terms = set(tree.get_terminals())
  allint = allid - terms

  # works for dated tips as well
  nh = nodeHeights(tree, allTipsZero = False)
  nhInOrder = sorted([(nh[x],x) for x in allint])
  nInOrder = [x[1] for x in reversed(nhInOrder)]

  # mapping for minimum node height
  mh = dict()
  for n in terms:
    mh[n] = nh[n]
  for h,n in nhInOrder:
    mh[n] = max([mh[c] for c in tree.node(n).succ])
  
  # x[0] is root
  
  sr = []
  sr.append("h%d = x[0]" % nInOrder[0])
  for i,k in enumerate(nInOrder[1:]):
    if mh[k] != 0 :
      m = "%g" % mh[k]
      sr.append("h%d = x[%d] * (h%d - %s) + %s" % (k, i+1, tree.node(k).prev, m, m))
    else :
      sr.append("h%d = x[%d] * h%d" % (k, i+1, tree.node(k).prev))
    
  for r,k in enumerate(clades) :
    n = clades[k][0][1] ; assert n != tree.root
    p = tree.node(n).prev
    if n in terms:
      if mh[n] != 0 :
        sr.append("b%d=(h%d - %g) # %d" % (r, p, mh[n], n))
      else :
        sr.append("b%d=(h%d) # %d" % (r, p, n))
    else :
      sr.append("b%d=(h%d - h%d) # %d" % (r, p, n, n))
  return sr, mh[tree.root]


def minPosteriorDistanceTree(tree, trees, limit = scipy.inf) :
  """ Find a branch length assignment for tree which minimizes the total
  distance to the set of trees.

  limit is an upper bound (presumably from a prior call here with another tree).
  If the distance is known to be larger, the optimization for this tree can be
  skipped.
  """
  
  treeParts = allPartitions(tree, [tree])

  posteriorParts = allPartitions(tree, trees)

  # Text of expression to compute the total distance. The variables are the
  # branch lengths. 
  ee = ""

  # Constant term of total distance (independent from branch lengths) 
  c0 = 0

  for r,k in enumerate(treeParts) :
    if k in posteriorParts :
      # A tree clade which appears in some posterior trees

      # All branch lengths from posterior for this clade
      br = [t.node(n).data.branchlength for t,n in posteriorParts[k]]
      
      # Number of posterior trees without the clade 
      a1 = (len(trees) - len(posteriorParts[k]))

      assert len(trees) == a1 + len(br)

      # Expanded form of the non constant part of [sum_i (b_r - b_ri)**2], where
      # b_ri is the branch length in the i'th posterior tree if the clade
      # exists, 0 otherwise.
      
      ee += "+ %.15f * b%d**2 + %.15f * b%d" % (a1 + len(br), r, -2*sum(br), r)

      # The constant term contribution 
      c0 += sum([x**2 for x in br])
    else :
      # A tree clade not appearing in posterior: contributes the full branch
      # length for each posterior tree.
      
      ee += "+(%d * b%d**2)" % (len(trees), r)

  # Total distance of branches terminating at a clade which is missing in tree.
  # This is (not necessarily a good) lower bound on the total distance.
  z0 = 0
  
  for k in posteriorParts :
    if k not in treeParts:
      a0 = sum([(t.node(n).data.branchlength-0)**2 for t,n in posteriorParts[k]])
      c0 += a0
      z0 += a0

  # Tuck in the constant
  ee = ("%.15g " % c0) + ee

  if z0 >= limit :
    return (None, 0.0)

  # Get code which transforms the heights encoding to branch lengths
  # A descendant height is specified as a fraction in [0,1] of its ancestor
  # height (but leading number is the root height).
  
  ba,minRoot = _treeBranchAssignmentExprs(tree, treeParts)

  # Define the posterior distance function on the fly.
  exec "def f(x):\n  " + "\n  ".join(ba) + "\n  return " + ee

  # Number of variables (heights)
  nx = len(tree.get_terminals())-1
  
  ## zz = scipy.optimize.fmin_tnc(f, [tree.height()] + [0.5]*(nx-1),
  ##                              approx_grad=1, bounds = [[0,None]] +
  ##                              [[0,1]]*(nx-1), xtol= 1e-14)

  maxfun = 15000
  zz = scipy.optimize.fmin_l_bfgs_b(f,
                                    [treeHeight(tree)] + [random.random() for k
                                                      in range(nx-1)],
                                    approx_grad=1,
                                    bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                    iprint=-1, maxfun=maxfun)
  if zz[2]['warnflag'] != 0 :
    print "WARNING:", zz[2]['task']

  # Function to get the branch lengths from optimized heights
  exec "def b(x):\n  " + "\n  ".join(ba) + "\n  " + \
       "return (" + ",".join(['b%d' % k for k in range(len(treeParts))]) + ")" 

  # Do not change tree passed as argument. Copy tree and set branch lengths of
  # copy.
  ss = copy.deepcopy(tree)

  brs = b(zz[0])
  for r,k in enumerate(treeParts) :
    ss.node(treeParts[k][0][1]).data.branchlength = brs[r]

  return (ss, f(zz[0]))

