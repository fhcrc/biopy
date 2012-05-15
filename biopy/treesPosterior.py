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

def _getNodeIDsDescendingHeight(tree, nid, includeTaxa = True) :
  node = tree.node(nid)
  isLeaf = len(node.succ) == 0
  
  r = [nid]
  if isLeaf:
    if not includeTaxa :
      r = []
  else :
    r.extend(reduce(lambda x,y : x+y,
                    [_getNodeIDsDescendingHeight(tree, n, includeTaxa)
                     for n in node.succ]))
  return r

def der(k, n, h) :
  return ",".join((['x[%d] * d_h%d_x[%d]' % (k,h,i) for i in range(k)] +
                   ["h%d" % h] +
                   ["0"]*(n-k-1)))


from numpy import array

def _treeBranchAssignmentExprs(tree, clades, fctr, nodesMinHeight = None,
                               withDerivative = False, paranoid = False) :
  """ nodesMinHeight: minimum height for internal nodes"""
  allid = set(tree.all_ids())
  # taxa
  terms = set(tree.get_terminals())
  # internal nodes
  allint = allid - terms

  # works for dated tips as well
  nh = nodeHeights(tree, allTipsZero = False)
  # node before its descendants
  nInOrder = _getNodeIDsDescendingHeight(tree, tree.root, includeTaxa=False)
  # reversed, child before parent
  nhInOrder = [(nh[x],x) for x in reversed(nInOrder)]

  #nhInOrder = sorted([(nh[x],x) for x in allint])
  #nInOrder = [x[1] for x in reversed(nhInOrder)]

  # mapping (per node) of minimum height of node, which is the max among all of
  # its descendants
  mh = dict()
  for n in terms:
    mh[n] = nh[n]
  for h,n in nhInOrder:
    mh[n] = max([mh[c] for c in tree.node(n).succ])

  if nodesMinHeight is not None :
    for n in nodesMinHeight:
      mh[n] = max(mh[n], nodesMinHeight[n])

  if fctr != 1 :
    for n in mh :
      mh[n] = mh[n] * fctr
    
  # x[0] is root
  
  sr = []
  if paranoid:
    # solver can send values out of range
    sr.append("x = [max(x[0],%f)] + [min(max(z,0),1.0) for z in x[1:]]" % mh[tree.root])
  sr.append("h%d = x[0]" % nInOrder[0])
  if withDerivative:
    sr.append("d_h%d_x = [1] + [0]*%d" % (nInOrder[0],len(nInOrder)-1))
    sr.append("nDer = %d" % len(nInOrder))
    
  for i,k in enumerate(nInOrder[1:]):
    h = tree.node(k).prev
    if mh[k] != 0 :
      m = "%.15g" % mh[k]
      sr.append("h%d = x[%d] * (h%d - %s) + %s" % (k, i+1, h, m, m))
    else :
      sr.append("h%d = x[%d] * h%d" % (k, i+1, h))
      
    if withDerivative:
      sr.append("d_h%d_x = [%s]" % (k,der(i+1, len(nInOrder), h)))
      if mh[k] != 0 :
        sr.append("d_h%d_x[%d] -= %s" % (k,i+1,m))
        
    if paranoid:
      sr.append("assert h%d >= 0" % k)
    
  for r,k in enumerate(clades) :
    n = clades[k][0][1] ; assert n != tree.root
    p = tree.node(n).prev
    if n in terms:
      if mh[n] != 0 :
        sr.append("b%d=(h%d - %.15g) # %d" % (r, p, mh[n], n))
      else :
        sr.append("b%d=(h%d) # %d" % (r, p, n))
        
      if withDerivative:
        sr.append("d_b%d_x = d_h%d_x" % (r,p))
    else :
      sr.append("b%d=(h%d - h%d) # %d" % (r, p, n, n))
      
      if withDerivative:
        sr.append("d_b%d_x = [u-v for u,v in zip(d_h%d_x, d_h%d_x)]" % (r,p,n))
        
    if paranoid:
      sr.append("assert b%d >= 0, (%d,b%d,h%d,x)" % (r,r,r,p))
    
  return sr, mh[tree.root]

ver = False

def minPosteriorDistanceTree(tree, trees, limit = scipy.inf, norm = True,
                             nodesMinHeight = None, withDerivative = False) :
  """ Find a branch length assignment for tree which minimizes the total
  distance to the set of trees.

  limit is an upper bound (presumably from a prior call here with another tree).
  If the distance is known to be larger, the optimization for this tree can be
  skipped.
  """
  
  treeParts = allPartitions(tree, [tree])

  posteriorParts = allPartitions(tree, trees)

  # For numerical stability sake in computing gradients, scale trees
  # so that mean root height is 1 
  fctr = len(trees)/ sum([treeHeight(x) for x in trees]) if norm else 1
  if ver: print fctr
  
  # Text of expression to compute the total distance. The variables are the
  # branch lengths. 
  ee = ""
  dee = ""
  
  # Constant term of total distance (independent from branch lengths) 
  c0 = 0

  for r,k in enumerate(treeParts) :
    if k in posteriorParts :
      # A tree clade which appears in some posterior trees

      # All branch lengths from posterior for this clade
      br = [t.node(n).data.branchlength * fctr for t,n in posteriorParts[k]]
      
      # Number of posterior trees without the clade 
      a1 = (len(trees) - len(posteriorParts[k]))

      assert len(trees) == a1 + len(br)

      # Expanded form of the non constant part of [sum_i (b_r - b_ri)**2], where
      # b_ri is the branch length in the i'th posterior tree if the clade
      # exists, 0 otherwise.
      
      ee += "+ %.15f * b%d**2 + %.15f * b%d" % (a1 + len(br), r, -2*sum(br), r)

      if withDerivative:
        dee += ("+ %.15f * 2 * b%d * d_b%d_x[k] + %.15f * d_b%d_x[k]" %
                (a1 + len(br), r, r, -2*sum(br), r))

      # The constant term contribution 
      c0 += sum([x**2 for x in br])
    else :
      # A tree clade not appearing in posterior: contributes the full branch
      # length for each posterior tree.
      
      ee += "+(%d * b%d**2)" % (len(trees), r)

      if withDerivative:
        dee += "+(%d * 2 * b%d * d_b%d_x[k])" % (len(trees), r, r)
      
  # Total distance of branches terminating at a clade which is missing in tree.
  # This is (not necessarily good) lower bound on the total distance.
  z0 = 0
  
  for k in posteriorParts :
    if k not in treeParts:
      a0 = sum([(t.node(n).data.branchlength * fctr - 0)**2 for t,n in posteriorParts[k]])
      c0 += a0
      z0 += a0

  # Tuck in the constant
  ee = ("%.15g " % c0) + ee

  if z0 >= limit :
    return (None, 0.0)

  # Get code which transforms the heights encoding to branch lengths
  # A descendant height is specified as a fraction in [0,1] of its ancestor
  # height (but leading number is the root height).
  
  ba,minRoot = _treeBranchAssignmentExprs(tree, treeParts, fctr,
                                          nodesMinHeight = nodesMinHeight,
                                          withDerivative = withDerivative)

  # Define the posterior distance function on the fly.
  cod = ("def f(x):\n  " + "\n  ".join(ba) + "\n  return " +
         (('(' + ee + ", array([(" + dee + ") for k in range(nDer)]) )")
          if withDerivative else ee))
  exec cod
  
  # Number of variables (heights)
  nx = len(tree.get_terminals())-1
  
  if ver :
    print "@@",nx, minRoot, treeHeight(tree) * fctr
    print cod
  
  ## zz = scipy.optimize.fmin_tnc(f, [tree.height()] + [0.5]*(nx-1),
  ##                              approx_grad=1, bounds = [[0,None]] +
  ##                              [[0,1]]*(nx-1), xtol= 1e-14)

  maxfun = 15000
  if 1 :
    zz = scipy.optimize.fmin_l_bfgs_b(f,
                                    [1 if norm else treeHeight(tree)] +
                                    [random.random() for k in range(nx-1)],
                                    #[.7 for k in range(nx-1)],
                                    approx_grad=0 if withDerivative else 1,
                                    bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                    iprint=-1, maxfun=maxfun)
    if zz[2]['warnflag'] != 0 :
      print "WARNING:", zz[2]['task']
    #print zz[0]
    
  if 0:
    zz = scipy.optimize.fmin_tnc(f,
                                 [treeHeight(tree)*fctr] +
                                 [random.random() for k in range(nx-1)],
                                 # [.7 for k in range(nx-1)],
                                 approx_grad=1,
                                 bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                 maxfun=maxfun,
                                 messages=0)
    assert zz[2] == 1

  # Function to get the branch lengths from optimized heights
  cod = "def b(x):\n  " + "\n  ".join(ba) + "\n  " + \
       "return (" + ",".join(['b%d' % k for k in range(len(treeParts))]) + ")" 
  exec cod
  if ver: print cod
  
  # Do not change tree passed as argument. Copy tree and set branch lengths of
  # the copy.
  ss = copy.deepcopy(tree)

  brs = b(zz[0])
  for r,k in enumerate(treeParts) :
    ss.node(treeParts[k][0][1]).data.branchlength = brs[r]/fctr

  return (ss, f(zz[0]))

