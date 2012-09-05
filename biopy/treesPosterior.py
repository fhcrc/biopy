## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

"""
Tree posterior distance optimization
====================================

"""

# really, really needs factoring
__all__ = ["minPosteriorDistanceTree", "minPosteriorRADistanceTree", "minPosteriorNorm1DistanceTree"]

import scipy, scipy.optimize, copy, random
from numpy import array

from treeMeasure import allPartitions
from treeutils import treeHeight, nodeHeights

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


def _treeBranchAssignmentExprs(tree, clades, fctr, nodesMinHeight = None,
                               withDerivative = False, withInit = False,
                               withHeights = False, paranoid = False,
                               warnings = True) :
  """ nodesMinHeight: minimum height (lower bound) for internal nodes"""
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
  if withInit :
    htox = []
  
  sr = []
  if paranoid:
    # solver can send values out of range
    sr.append("x = [max(x[0],%f)] + [min(max(z,0),1.0) for z in x[1:]]" % mh[tree.root])
  sr.append("h%d = x[0]" % nInOrder[0])
  if withInit:
    htox.append("x[0] = %.15g" % nh[nInOrder[0]])
  
  if withDerivative:
    sr.append("d_h%d_x = [1] + [0]*%d" % (nInOrder[0],len(nInOrder)-1))
    sr.append("nDer = %d" % len(nInOrder))
    
  for i,k in enumerate(nInOrder[1:]):
    h = tree.node(k).prev
    if mh[k] != 0 :
      m = "%.15g" % mh[k]
      sr.append("h%d = x[%d] * (h%d - %s) + %s" % (k, i+1, h, m, m))
      if withInit:
        htox.append("x[%d] = %.15g" % (i+1, min(max((nh[k] - mh[k])/(nh[h] - mh[k]), 0), 1)))
    else :
      sr.append("h%d = x[%d] * h%d" % (k, i+1, h))
      if withInit:
        htox.append("x[%d] = %.15g" % (i+1, min(max(nh[k] / nh[h], 0),1)))
      
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
        if withHeights: sr.append("bs%d = %.15g" % (r, mh[n]))
      else :
        sr.append("b%d=(h%d) # %d" % (r, p, n)) 
        if withHeights: sr.append("bs%d = 0" % (r))
       
      if withDerivative:
        sr.append("d_b%d_x = d_h%d_x" % (r,p))
    else :
      sr.append("b%d=(h%d - h%d) # %d" % (r, p, n, n))
      if withHeights: sr.append("bs%d = h%d" % (r,n))
      
      if withDerivative:
        sr.append("d_b%d_x = [u-v for u,v in zip(d_h%d_x, d_h%d_x)]" % (r,p,n))
        
    if paranoid:
      sr.append("assert b%d >= 0, (%d,b%d,h%d,x)" % (r,r,r,p))

  if withInit:
    return sr, mh[tree.root], htox
  return sr, mh[tree.root]

verbose = False

def minPosteriorDistanceTree(tree, trees, limit = scipy.inf, norm = True,
                             nodesMinHeight = None, withDerivative = False,
                             factr=10000000.0, warnings = True) :
  """ Find a branch length assignment for tree which minimizes the total
  distance to the set of trees.

  limit is an upper bound (presumably from a prior call here with another tree).
  If the distance is known to be larger, the optimization for this tree can be
  skipped.
  """
  
  treeParts = allPartitions(tree, [tree])

  posteriorParts = allPartitions(tree, trees, lambda t,n : t.node(n).data.branchlength)

  # For numerical stability sake in computing gradients, scale trees
  # so that mean root height is 1 
  fctr = len(trees)/ sum([treeHeight(x) for x in trees]) if norm else 1
  if verbose: print fctr
  
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
      #br = [t.node(n).data.branchlength * fctr for t,n in posteriorParts[k]]
      br = [b * fctr for b in posteriorParts[k]]
      
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
      #a0 = sum([(t.node(n).data.branchlength * fctr - 0)**2 for t,n in
      #posteriorParts[k]])
      a0 = sum([(b * fctr - 0)**2 for b in posteriorParts[k]])      
      c0 += a0
      z0 += a0

  del posteriorParts
  
  # Tuck in the constant
  ee = ("%.15g " % c0) + ee

  if z0 >= limit :
    return (None, 0.0)

  # Get code which transforms the heights encoding to branch lengths
  # A descendant height is specified as a fraction in [0,1] of its ancestor
  # height (but leading number is the root height).
  
  ba,minRoot = _treeBranchAssignmentExprs(tree, treeParts, fctr,
                                          nodesMinHeight = nodesMinHeight,
                                          withDerivative = withDerivative,
                                          withInit = False)
  # ba,minRoot,htox
  
  # Define the posterior distance function on the fly.
  cod = ("def f(x):\n  " + "\n  ".join(ba) + "\n  return " +
         (('(' + ee + ", array([(" + dee + ") for k in range(nDer)]) )")
          if withDerivative else ee))
  exec cod

  # Number of variables (heights)
  nx = len(tree.get_terminals())-1


#    xcod = "def htox():\n  x = [0]*%d\n  " % nx + "\n  ".join(htox) \
#         + "\n  return x"
#    exec xcod
  
  if verbose :
    print "@@",nx, minRoot, treeHeight(tree) * fctr
    print cod
  
  ## zz = scipy.optimize.fmin_tnc(f, [tree.height()] + [0.5]*(nx-1),
  ##                              approx_grad=1, bounds = [[0,None]] +
  ##                              [[0,1]]*(nx-1), xtol= 1e-14)

  maxfun = 15000
  if 1 :
    x0 = [1 if norm else treeHeight(tree)] + [random.random() for k in
                                              range(nx-1)]
    #[.7 for k in range(nx-1)],
    
    #print f(x0)[0], 
    if 0 :
      x0 = htox()
      if norm:
        x0[0] = 1
      #print f(x0)[0]

    assert x0[0] >= minRoot
    zz = scipy.optimize.fmin_l_bfgs_b(f, x0,
                                      approx_grad=0 if withDerivative else 1,
                                      bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                      factr = factr,
                                      iprint=-1, maxfun=maxfun)
    if warnings and zz[2]['warnflag'] != 0 :
      print "WARNING:", zz[2]['task']
    #print zz[0]
    
  if 0:
    zz = scipy.optimize.fmin_tnc(f,
                                 [treeHeight(tree)*fctr] +
                                 [random.uniform(.8,.9) for k in range(nx-1)],
                                 #[.8 for k in range(nx-1)],
                                 approx_grad=1,
                                 bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                 maxfun=maxfun,
                                 messages=0)
    assert zz[2] == 1

  # Function to get the branch lengths from optimized heights
  cod = "def b(x):\n  " + "\n  ".join(ba) + "\n  " + \
       "return (" + ",".join(['b%d' % k for k in range(len(treeParts))]) + ")" 
  exec cod
  if verbose: print cod
  
  # Do not change tree passed as argument. Copy tree and set branch lengths of
  # the copy.
  ss = copy.deepcopy(tree)

  brs = b(zz[0])
  for r,k in enumerate(treeParts) :
    ss.node(treeParts[k][0][1]).data.branchlength = brs[r]/fctr

  val = f(zz[0])
  if withDerivative :
    val = val[0]
  return (ss, val)



# def _sumNonIntersect(l1,u1,lu2) :
#   tot = 0.0
#   b1 = (u1-l1)
#   for l2,u2 in lu2 :
#     b1b2 = b1 + (u2-l2)
#     tot += min( b1b2 + 2*(max(l1,l2) - min(u1,u2)), b1b2)
#   return tot

from cchelp import sumNonIntersect as _sumNonIntersect

def minPosteriorRADistanceTree(tree, trees, limit = scipy.inf, norm = True,
                               nodesMinHeight = None, withDerivative = False,
                               withInit = True, factr=10000000.0) :
  """ Find a branch length assignment for tree which minimizes the total
  distance to the set of trees.

  limit is an upper bound (presumably from a prior call here with another tree).
  If the distance is known to be larger, the optimization for this tree can be
  skipped.
  """
  assert not withDerivative
  # not correct for tip/node lower bounds
  assert nodesMinHeight is None
  
  treeParts = allPartitions(tree, [tree])

  posteriorParts = allPartitions(tree, trees,
                                 func = lambda t,(n,h) : (h, t.node(n).data.branchlength),
                                 withHeights = True)

  # For numerical stability sake in computing gradients, scale trees
  # so that mean root height is 1 
  fctr = len(trees)/ sum([treeHeight(x) for x in trees]) if norm else 1
  if verbose: print fctr
  
  # Text of expression to compute the total distance. The variables are the
  # branch lengths. 
  ee = ""
  dee = ""
  
  for r,k in enumerate(treeParts) :
    if k in posteriorParts :
      # A tree clade which is in some posterior trees

      # Branchs from posterior for this clade
      br = [(h*fctr,b*fctr) for h,b in posteriorParts[k]]
      
      # Number of posterior trees without the clade 
      a1 = (len(trees) - len(posteriorParts[k]))

      assert len(trees) == a1 + len(br)

      dt = '[' + ','.join(["(%.15g,%.15g)" % (h,h+b) for h,b in br]) + ']'
      ee += "+ _sumNonIntersect(bs%d,bs%d+b%d,%s)" % (r,r,r,dt)
      ee += "+ %d * b%d" % (a1,r)
      
      if withDerivative:
        pass
        #dee += ("+ %.15f * 2 * b%d * d_b%d_x[k] + %.15f * d_b%d_x[k]" %
        #        (a1 + len(br), r, r, -2*sum(br), r))

      # The constant term contribution 
      #c0 += sum([x**2 for x in br])
    else :
      # A tree clade not appearing in posterior: contributes the full branch
      # length for each posterior tree.
      
      ee += "+(%d * b%d)" % (len(trees), r)

      if withDerivative:
        pass
        #dee += "+(%d * 2 * b%d * d_b%d_x[k])" % (len(trees), r, r)
      
  
  # Constant term of total distance (independent from branch lengths) 
  c0 = 0
  
  for k in posteriorParts :
    if k not in treeParts:
      c0 += sum([b * fctr for h,b in posteriorParts[k]])      

  # Total distance of branches terminating at a clade which is missing in tree.
  # This is (not necessarily good) lower bound on the total distance.
  z0 = c0

  del posteriorParts
  
  # Tuck in the constant
  ee = ("%.15g " % c0) + ee

  if z0 >= limit :
    return (None, 0.0)

  # Get code which transforms the heights encoding to branch lengths
  # A descendant height is specified as a fraction in [0,1] of its ancestor
  # height (but leading number is the root height).
  ba,minRoot,htox = _treeBranchAssignmentExprs(tree, treeParts, fctr,
                                               nodesMinHeight = nodesMinHeight,
                                               withDerivative = withDerivative,
                                               withInit = True,
                                               withHeights = True)
    
  # Define the posterior distance function on the fly.
  cod = ("def rascore(x):\n  " + "\n  ".join(ba) + "\n  return " +
         (('(' + ee + ", array([(" + dee + ") for k in range(nDer)]) )")
          if withDerivative else ee))
  exec cod
  if verbose: print cod

  # Number of variables (heights)
  nx = len(tree.get_terminals())-1

  xcod = "def htox():\n  x = [0]*%d\n  " % nx + "\n  ".join(htox) \
         + "\n  return x"
  exec xcod

  # Function to get the branch lengths from optimized heights
  cod = "def code2branches(x):\n  " + "\n  ".join(ba) + "\n  " + \
       "return (" + ",".join(['b%d' % k for k in range(len(treeParts))]) + ")" 
  exec cod
  if verbose: print cod

  if verbose :
    print "@@",nx, minRoot, treeHeight(tree) * fctr
    print cod
  
  maxfun = 15000
  if 1 :
    # small protection against disasters
    while True:
      if withInit :
        x0 = htox()
        if norm:
          x0[0] = 1
      else :
        x0 = [1 if norm else treeHeight(tree)] + \
             [random.random() for k in range(nx-1)]

      initialVal = rascore(x0)
      assert x0[0] >= minRoot
      zz = scipy.optimize.fmin_l_bfgs_b(rascore, x0,
                                        approx_grad=0 if withDerivative else 1,
                                        bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                        factr = factr,
                                        iprint=-1, maxfun=maxfun)
      if warnings and zz[2]['warnflag'] != 0 :
        print "WARNING:", zz[2]['task']
        
      finaleVal = rascore(zz[0])
      if finaleVal < initialVal :
        break
      withInit = False
      factr /= 10
      if factr < 1e6 :
        # failed, leave as is
        zz = htox()
        finaleVal = rascore(zz[0])
  else :
    zz = scipy.optimize.fmin_tnc(rascore,
                                 [treeHeight(tree)*fctr] +
                                 [random.uniform(.8,.9) for k in range(nx-1)],
                                 #[.8 for k in range(nx-1)],
                                 approx_grad=1,
                                 bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                 maxfun=maxfun,
                                 messages=0)
    assert zz[2] == 1

  
  # Do not change tree passed as argument. Copy tree and set branch lengths of
  # the copy.
  ss = copy.deepcopy(tree)

  brs = code2branches(zz[0])
  for r,k in enumerate(treeParts) :
    ss.node(treeParts[k][0][1]).data.branchlength = brs[r]/fctr

  val = finaleVal
  if withDerivative :
    val = val[0]
  return (ss, val)


# sort branches and pre compute cumulative sum, so we can compute distance using
# a binary search

def _prepare(bs) :
  sbs = sorted(bs)
  tot = 0
  cums = [0]
  for b in sbs :
    tot += b
    cums.append(tot)
  return  "((" + ','.join(["%.15g" % b for b in sbs]) + ',),(' + \
         ','.join(["%.15g" % s for s in cums]) + ',))'

import bisect
def _absDiffBranch(b1,(sbs,cums)) :
  if len(sbs) <= 16 :
    tot = 0
    for b in sbs :
      tot += abs(b1 - b)
  else :
    i = bisect.bisect(sbs, b1)
    tot = b1*(i+i-len(sbs)) + cums[-1] - 2*cums[i]
  return tot

## import numpy.linalg
## from treeMeasure import branchScoreTreeDistance
## def truebs1(tree, trees) :
##   vd = lambda a1,a2 : numpy.linalg.norm([x-y for x,y in zip(a1,a2)], 1)
##   return sum(branchScoreTreeDistance(tree, x, vd) for x in trees)
## def v1score_ck(v1score, zz0, tree, trees) :
##   global ssx, code2branches, treeParts, fctr
  
##   t1 = truebs1(tree, trees)*fctr
##   t2 = v1score(zz0)
##   #import pdb ; pdb.set_trace()
##   assert abs(t1-t2)/t1 < 1e-8, ( t1,t2, t1-t2)
##   return t2

## def totr(zz0) :
##   global ssx, code2branches, treeParts, fctr
##   brs = code2branches(zz0)
##   for r,k in enumerate(treeParts) :
##     ssx.node(treeParts[k][0][1]).data.branchlength = brs[r]/fctr
##   return ssx

def minPosteriorNorm1DistanceTree(tree, trees, limit = scipy.inf, norm = True,
                                  nodesMinHeight = None, withDerivative = False,
                                  withInit = True, factr=10000000.0,
                                  warnings = True) :
  """ Find a branch length assignment for tree which minimizes the total
  distance to the set of trees.

  limit is an upper bound (presumably from a prior call here with another tree).
  If the distance is known to be larger, the optimization for this tree can be
  skipped.
  """
  #global ssx, treeParts, fctr
  
  assert not withDerivative
  # not correct for tip/node lower bounds
  assert nodesMinHeight is None
  
  treeParts = allPartitions(tree, [tree])

  posteriorParts = allPartitions(tree, trees,
                                 func = lambda t,n : t.node(n).data.branchlength,
                                 withHeights = False)

  # For numerical stability sake in computing gradients, scale trees
  # so that mean root height is 1 
  fctr = len(trees)/ sum([treeHeight(x) for x in trees]) if norm else 1

  if verbose: print fctr
  
  # Text of expression to compute the total distance. The variables are the
  # branch lengths. 
  ee = ""
  dee = ""
  
  for r,k in enumerate(treeParts) :
    if k in posteriorParts :
      # A tree clade which is in some posterior trees

      # Branchs from posterior for this clade
      br = [b*fctr for b in posteriorParts[k]]
      
      # Number of posterior trees without the clade 
      a1 = (len(trees) - len(posteriorParts[k]))

      assert len(trees) == a1 + len(br)

      ee += "+ _absDiffBranch(b%d,%s)" % (r,_prepare(br))
      ee += "+ %d * b%d" % (a1,r)
      
      if withDerivative:
        pass
        #dee += ("+ %.15f * 2 * b%d * d_b%d_x[k] + %.15f * d_b%d_x[k]" %
        #        (a1 + len(br), r, r, -2*sum(br), r))

      # The constant term contribution 
      #c0 += sum([x**2 for x in br])
    else :
      # A tree clade not appearing in posterior: contributes the full branch
      # length for each posterior tree.
      
      ee += "+(%d * b%d)" % (len(trees), r)

      if withDerivative:
        pass
        #dee += "+(%d * 2 * b%d * d_b%d_x[k])" % (len(trees), r, r)
      
  # Constant term of total distance (independent from branch lengths) 
  c0 = 0
  
  for k in posteriorParts :
    if k not in treeParts:
      c0 += sum([b * fctr for b in posteriorParts[k]])      

  # Total distance of branches terminating at a clade which is missing in tree.
  # This is (not necessarily good) lower bound on the total distance.
  z0 = c0

  del posteriorParts
  
  # Tuck in the constant
  ee = ("%.15g " % c0) + ee

  if z0 >= limit :
    return (None, 0.0)

  # Get code which transforms the heights encoding to branch lengths
  # A descendant height is specified as a fraction in [0,1] of its ancestor
  # height (but leading number is the root height).
  ba,minRoot,htox = _treeBranchAssignmentExprs(tree, treeParts, fctr,
                                               nodesMinHeight = nodesMinHeight,
                                               withDerivative = withDerivative,
                                               withInit = True,
                                               withHeights = True)
    
  #global v1score
  # Define the posterior distance function on the fly.
  cod = ("def v1score(x):\n  " + "\n  ".join(ba) + "\n  return " +
         (('(' + ee + ", array([(" + dee + ") for k in range(nDer)]) )")
          if withDerivative else ee))
  exec cod
  #v1score = v1score1
  if verbose: print cod

  # Number of variables (heights)
  nx = len(tree.get_terminals())-1

  xcod = "def htox():\n  x = [0]*%d\n  " % nx + "\n  ".join(htox) \
         + "\n  return x"
  exec xcod

  global code2branches
  # Function to get the branch lengths from optimized heights
  codb = "def code2branchesa(x):\n  " + "\n  ".join(ba) + "\n  " + \
       "return (" + ",".join(['b%d' % k for k in range(len(treeParts))]) + ")" 
  exec codb
  if verbose: print cod
  code2branches = code2branchesa
  if verbose :
    print "@@",nx, minRoot, treeHeight(tree) * fctr
    print cod

  if 0 :
    cod8 = "def tv1score(tree):\n"
    for r,k in enumerate(treeParts) :
      cod8 += "  b%d = tree.node(%d).data.branchlength\n" % (r, treeParts[k][0][1])
    cod8 += "  return " +  ee
    exec cod8
  
  #ssx = copy.deepcopy(tree)
  #global xtrees
  #xtrees = trees
  #exec """def v1scorex(zz0) : return v1score_ck(v1score, zz0, totr(zz0), xtrees)"""
  
  maxfun = 15000
  if 1 :
    # small protection against disasters
    while True:
      if withInit :
        x0 = htox()
        if norm:
          x0[0] = 1
      else :
        x0 = [1 if norm else treeHeight(tree)] + \
             [random.random() for k in range(nx-1)]

      initialVal = v1score(x0)
      assert x0[0] >= minRoot
      zz = scipy.optimize.fmin_l_bfgs_b(v1score, x0,
                                        approx_grad=0 if withDerivative else 1,
                                        bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                        factr = factr,
                                        iprint=-1, maxfun=maxfun)
      if warnings and zz[2]['warnflag'] != 0 :
        print "WARNING:", zz[2]['task']
        
      finaleVal = v1score(zz[0])
      if finaleVal < initialVal :
        break
      withInit = False
      factr /= 10
      if factr < 1e6 :
        # failed, leave as is
        zz = htox()
        finaleVal = v1score(zz[0])
  else :
    zz = scipy.optimize.fmin_tnc(v1score,
                                 [treeHeight(tree)*fctr] +
                                 [random.uniform(.8,.9) for k in range(nx-1)],
                                 #[.8 for k in range(nx-1)],
                                 approx_grad=1,
                                 bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                 maxfun=maxfun,
                                 messages=0)
    assert zz[2] == 1

  
  # Do not change tree passed as argument. Copy tree and set branch lengths of
  # the copy.
  ss = copy.deepcopy(tree)

  brs = code2branches(zz[0])
  for r,k in enumerate(treeParts) :
    ss.node(treeParts[k][0][1]).data.branchlength = brs[r]/fctr

  val = finaleVal
  if withDerivative :
    val = val[0]
  return (ss, val/fctr if norm else val)

#  return (ss, val/fctr if norm else val,
#          (v1score, htox, code2branches, fctr, tv1score, cod8, zz[0]))

