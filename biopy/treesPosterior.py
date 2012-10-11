## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

from __future__ import division

"""
Tree posterior distance optimization
====================================

Find a centeroid tree for a set of trees, a tree which minimizes the sum of
distances to all trees in the set under some tree distance.

4 distances are available:
BRANCH_SCORE: 
BRANCH_SCORE_2:
HEIGHTS_SCORE:
HEIGHTS_ONLY:
ROOTED_AGREEMENT:

Fixed target topology.
"""


__all__ = ["minDistanceTree", "BRANCH_SCORE" , "BRANCH_SCORE_2" ,
           "HEIGHTS_SCORE" , "ROOTED_AGREEMENT", "HEIGHTS_ONLY" ]

BRANCH_SCORE , BRANCH_SCORE_2 , HEIGHTS_SCORE , ROOTED_AGREEMENT, HEIGHTS_ONLY = range(5)

import scipy, scipy.optimize, copy, random
from scipy.optimize import fminbound

from numpy import array, median

from treeMeasure import allPartitions
from treeutils import treeHeight, nodeHeights, getPreOrder, getPostOrder

from cchelp import sumNonIntersect as _sumNonIntersect, sumNonIntersectDer as _sumNonIntersectDer

verbose = False


def _setTreeHeights(tree, opts) :
  order = getPostOrder(tree)
  for node in order :
    if not node.succ:
      node.data.height = 0
    else :
      hs = [tree.node(x).data.height for x in node.succ]
      mn = sum([h + opts[x] for x,h in zip(node.succ,hs)])/len(hs)
      node.data.height = max(*(hs + [mn]))

  for n in tree.all_ids() :
    node = tree.node(n)
    if node.prev is not None:
      p = tree.node(node.prev) 
      node.data.branchlength = p.data.height - node.data.height
      assert node.data.branchlength >= 0
  return tree

def _setTreeHeightsForTargets(tree, ftargets) :
  for i,h in ftargets() :
    tree.node(i).data.height = h
  for i in tree.get_terminals() :
     tree.node(i).data.height = 0

  for n in tree.all_ids() :
    node = tree.node(n)
    if node.prev is not None:
      p = tree.node(node.prev) 
      node.data.branchlength = p.data.height - node.data.height
      assert node.data.branchlength >= 0
  for i in tree.all_ids() :
    del tree.node(i).data.height
  
# Sort branches and pre-compute cumulative sum, so distance can be computed
# in O(lg(n)), using a binary search, instead of O(n).

def _prepare(bs, stringOnly = True) :
  sbs = sorted(bs)
  tot = 0
  cums = [0]
  for b in sbs :
    tot += b
    cums.append(tot)
  asStr =  "((" + ','.join(["%.15g" % b for b in sbs]) + ',),(' + \
          ','.join(["%.15g" % s for s in cums]) + ',))'
  if stringOnly :
    return asStr
  return (asStr, (sbs,cums))

import bisect
def _absDiffBranch(b1,(sbs,cums)) :
  # Assume binary search overhead for small lists is 16 (not verified)
  if len(sbs) <= 16 :
    tot = 0
    for b in sbs :
      tot += abs(b1 - b)
  else :
    i = bisect.bisect(sbs, b1)
    tot = b1*(i+i-len(sbs)) + cums[-1] - 2*cums[i]
  return tot

# Value and derivative
def _absDiffBranchDer(b1,(sbs,cums)) :
  if len(sbs) <= 16 :
    tot,dtot = 0,0
    for b in sbs :
      if b1 >= b :
        tot += b1 - b
        dtot += 1
      else :
        tot += b - b1
        dtot -= 1
  else :
    i = bisect.bisect(sbs, b1)
    dtot = (i+i-len(sbs))
    tot = b1*dtot + cums[-1] - 2*cums[i]
  return tot,dtot


def der(k, n, h) :
  return ",".join((['x[%d] * d_h%d_x[%d]' % (k,h,i) for i in range(k)] +
                   ["h%d" % h] +
                   ["0"]*(n-k-1)))

def hs2ratio(h, hpar, hmin) :
  assert hpar >= h >= hmin
  
  dp = hpar - hmin
  dh = h - hmin
  if dp == 0 :
    return 0.5
  
  return min(max((h - hmin)/(hpar - hmin), 0), 1)
  
def _treeBranchAssignmentExprs(tree, clades, fctr, nodesMinHeight = None,
                               withDerivative = False, withInit = False,
                               withHeights = False, paranoid = False) :
  """ nodesMinHeight: minimum height (lower bound) for internal nodes
  paranoid: add internal consistency checks
  """
  
  allid = set(tree.all_ids())
  # taxa
  terms = set(tree.get_terminals())
  # internal nodes
  allint = allid - terms

  # Works for dated tips as well
  nh = nodeHeights(tree, allTipsZero = False)
  # Node before its descendants
  nInOrder = getPreOrder(tree, includeTaxa=False)
  # Reversed, child before parent
  nhInOrder = [(nh[x],x) for x in reversed(nInOrder)]

  # Mapping (per node) of minimum height of node, which is the max among all of
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
        htox.append("x[%d] = %.15g" % (i+1, hs2ratio(nh[k], nh[h], mh[k])))
    else :
      sr.append("h%d = x[%d] * h%d" % (k, i+1, h))
      if withInit:
        htox.append("x[%d] = %.15g" % (i+1,  hs2ratio(nh[k], nh[h], 0)))
      
    if withDerivative:
      sr.append("d_h%d_x = [%s]" % (k,der(i+1, len(nInOrder), h)))
      if mh[k] != 0 :
        sr.append("d_h%d_x[%d] -= %s" % (k,i+1,m))
        
    if paranoid:
      sr.append("assert h%d >= 0" % k)
    
  for r,k in enumerate(clades) :
    n = clades[k][0][1]                                ; assert n != tree.root
    p = tree.node(n).prev
    if n in terms:
      if mh[n] != 0 :
        sr.append("b%d=(h%d - %.15g) # %d" % (n, p, mh[n], n))
        if withHeights: sr.append("bs%d = %.15g" % (n, mh[n]))
      else :
        sr.append("b%d=(h%d) # %d" % (n, p, n)) 
        if withHeights: sr.append("bs%d = 0" % (n))
       
      if withDerivative:
        sr.append("d_b%d_x = d_h%d_x" % (n,p))
    else :
      sr.append("b%d=(h%d - h%d) # %d" % (n, p, n, n))
      if withHeights: sr.append("bs%d = h%d" % (n,n))
      
      if withDerivative:
        sr.append("d_b%d_x = [u-v for u,v in zip(d_h%d_x, d_h%d_x)]" % (n,p,n))
        
    if paranoid:
      sr.append("assert b%d >= 0, (%d,b%d,h%d,x)" % (n,r,r,p))

  if withInit:
    return sr, mh[tree.root], htox
  return sr, mh[tree.root]

def minDistanceTree(method, tree, trees, limit = scipy.inf, norm = True,
                    nodesMinHeight = None, withDerivative = True,
                    initMethod = "opt", factr=1e7, warnings = False, internals = False) :
  """ Find a branch length assignment for tree which minimizes the total
  distance to the set of trees.

  limit is an upper bound (presumably from a prior call here with another tree).
  If the distance is known to be larger, the optimization for this tree can be
  skipped.

  initMethod : "tree", "random" or "opt"  
  """
  
  # Not correct for tip/node lower bounds
  assert nodesMinHeight is None or (method == BRANCH_SCORE)

  # Do not change tree passed as argument. Copy tree and set branch lengths of
  # the copy.
  tree = copy.deepcopy(tree)

  usesHeights = (method not in [BRANCH_SCORE,BRANCH_SCORE_2])
                      
  treeParts = allPartitions(tree, [tree])

  hsOnly = method == HEIGHTS_ONLY
  if usesHeights:
    if hsOnly :
      func = lambda t,(n,h) : h
    else :
      func = lambda t,(n,h) : (h, t.node(n).data.branchlength)
  else :
    func = lambda t,n : t.node(n).data.branchlength
  
  posteriorParts = allPartitions(tree, trees, func = func,
                                 withHeights = usesHeights, withRoot = hsOnly)
  if hsOnly :
    posteriorParts, rhs = posteriorParts
    
  # For numerical stability in computing gradients, scale trees
  # to get a mean root height of 1 
  fctr = len(trees) / sum([treeHeight(x) for x in trees]) if norm else 1
  if verbose: print fctr
  
  # Expression (as text) for computing the total distance. Variables are 
  # branch lengths and/or heights. 
  ee = ""
  # Expression (as text) for computing the derivative 
  dee = ""
  # Expressions (as text) for parts with multiple use (efficiency measure) 
  pee = ""

  # used only by HEIGHTS_SCORE
  rootDone = False

  # Constant term of total distance (independent from branch lengths) 
  c0 = 0

  # Optimal length of branch according to BRANCH_SCORE. Used to get a quick and
  # rough tarting point.
  optBranches = dict()

  if hsOnly :
    targets = "("
  
  for r,k in enumerate(treeParts) :
    # Node id
    nn = treeParts[k][0][1]
    
    if k in posteriorParts :
      # A clade in target tree which is present in some trees of the set.

      # Branchs/Heights from tree set for this clade
      if usesHeights :
        if hsOnly :
          hs = [h*fctr for h in posteriorParts[k]]
        else :
          br = [(h*fctr,b*fctr) for h,b in posteriorParts[k]]
      else :
        br = [b*fctr for b in posteriorParts[k]]
      
      # Number of trees in set which do not have the clade 
      a1 = len(trees) - len(posteriorParts[k])

      if method == BRANCH_SCORE_2:
        # Expanded form of the non constant part of [sum_i (b_r - b_ri)**2], where
        # b_ri is the branch length in the i'th posterior tree if the clade
        # exists, 0 otherwise.

        ee += "+ %.15f * b%d**2 + %.15f * b%d" % (a1 + len(br), nn, -2*sum(br), nn)

        if withDerivative:
          dee += ("+ %.15f * 2 * b%d * d_b%d_x[k] + %.15f * d_b%d_x[k]" %
                  (a1 + len(br), nn, nn, -2*sum(br), nn))

        # The constant term contribution 
        c0 += sum([x**2 for x in br])

      elif method == BRANCH_SCORE:
        vlsas,vls = _prepare(br, False)
        if withDerivative :
          pee += "  ab%d,abd%d = _absDiffBranchDer(b%d,%s)\n" % (nn,nn,nn,vlsas)
          pee += "  abd%d += %d\n" % (nn,a1)

          dee += "+(abd%d * d_b%d_x[k])" % (nn,nn)        
          ee += "+ ab%d" % (nn,)
        else :
          # term for trees with the clade
          ee += "+ _absDiffBranch(b%d,%s)" % (nn,vlsas)

        # term for trees without clade
        ee += "+ %d * b%d" % (a1,nn)
        
      elif method == HEIGHTS_SCORE :
        if not tree.node(nn).data.taxon  :
          vlsas = _prepare([x[0] for x in br])
          if withDerivative :
            pee += "  v%d,dv%d = _absDiffBranchDer(bs%d,%s)\n" % (nn,nn,nn,vlsas)

            ee += "+ v%d" % (nn,)
            if a1 > 0 :
              dee += "+(dv%d * d_h%d_x[k] + %d * d_b%d_x[k])" % (nn,nn,a1,nn)
            else :
              dee += "+(dv%d * d_h%d_x[k])" % (nn,nn)        
          else :
            ee += "+ _absDiffBranch(bs%d,%s)" % (nn,vlsas)

          if a1 > 0 :
            ee += "+ %d * b%d" % (a1,nn)

      elif method == HEIGHTS_ONLY :
        
        if not tree.node(nn).data.taxon  :
          hTarget = median(hs)
          if withDerivative :
            pee += "  dv%d = 1 if bs%d > %.14f else -1\n" % (nn,nn,hTarget)
            dee += "+(dv%d * d_h%d_x[k])" % (nn,nn)        

          ee += "+ abs(bs%d - %.14f)" % (nn,hTarget)

          targets += "(%d,%.14f)," % (nn,hTarget)
        
      elif method == ROOTED_AGREEMENT :
        dt = '[' + ','.join(["(%.15g,%.15g)" % (h,h+b) for h,b in br]) + ']'
        if withDerivative :
          # value, d(value)/d(branch-start), d(value)/d(branch-end-height)
          # start/end as usual are reversed to the tree direction, start has
          # lower height.
          
          pee += "  val%d,dvdl%d,dvdh%d = _sumNonIntersectDer(bs%d,bs%d+b%d,%s)\n" % ((nn,)*6 + (dt,))

          # make low the sum of low+high, save operation in loops
          pee += "  dvdl%d += dvdh%d\n" % (nn,nn)

          # fold constant into high, same reason
          pee += "  dvdh%d += %d\n" % (nn,a1)

          ee += "+ val%d" % (nn,)

          # brevity: dx == dx_k
          #
          # dV/dx = dV(l,u)/dx = dV(l)/dx + dV(u)/dx = dV(l)/dl * dl/dx + dV(u)/du * du/dx

          if tree.node(nn).data.taxon :
            dee += "+(dvdh%d * d_b%d_x[k])" % (nn,nn)     
          else :
            # rearranged and simplified
            dee += "+(dvdl%d * d_h%d_x[k] + dvdh%d * d_b%d_x[k])" % (nn,nn,nn,nn)     

        else :
          ee += "+ _sumNonIntersect(bs%d,bs%d+b%d,%s)" % (nn,nn,nn,dt)

        if a1 > 0 :
          ee += "+ %d * b%d" % (a1,nn)
      else :
        raise RuntimeError("Invalid method %d" % method)

      if initMethod == "opt" and not hsOnly :
        if method != BRANCH_SCORE :
          brx = [x[1] for x in br] if usesHeights else br
          vls = _prepare(brx, False)[1]
        # else vls is ready (been computed already)
        f = lambda b : _absDiffBranch(b, vls) + a1*b
        optBranches[nn] = fminbound(f, 0, vls[0][-1])
          
    else :
      # A tree clade not appearing in posterior: contributes the full branch
      # length for each posterior tree.
      if method == BRANCH_SCORE_2:
        ee += "+(%d * b%d**2)" % (len(trees), nn)

        if withDerivative:
          dee += "+(%d * 2 * b%d * d_b%d_x[k])" % (len(trees), nn, nn)
      else :
        # all linear distances add the missing branches as is
        ee += "+(%d * b%d)" % (len(trees), nn)

        if withDerivative:
          dee += "+(%d * d_b%d_x[k])" % (len(trees), nn)

      if initMethod == "opt" :
        optBranches[nn] = 0
        
    # Heights score need special treatment to include the root term.
    if (method == HEIGHTS_SCORE or hsOnly) and \
           (tree.node(nn).prev == tree.root and not rootDone):
      if hsOnly :
        rTarget = median(rhs) * fctr
        if withDerivative :
          pee += "  dvroot = 1 if h0 > %.14f else -1\n" % rTarget
          dee += "+(dvroot * (k==0) )"

        ee += "+ abs(h0 - %.14f)" % (rTarget)
        targets += "(%d,%.14f)," % (tree.root,rTarget)
      else :  
        rhs = [(b+h)*fctr for h, b in posteriorParts[k]]
        vlsas,vls = _prepare(rhs, False)
        if withDerivative :
          pee += "  vroot,dvroot = _absDiffBranchDer(h0,%s)\n" % vlsas

          ee += "+ vroot"
          dee += "+(dvroot * (k==0) )"

        else :
          ee += "+ _absDiffBranch(h0,%s)" % vlsas

      rootDone = True

  # Total distance of branches terminating at a clade which is missing in tree.
  # This is (not necessarily good) lower bound on the total distance.
  z0 = 0
  if not hsOnly:
    for k in posteriorParts :
      if k not in treeParts:
        if method == BRANCH_SCORE_2:
          f = lambda b : (b * fctr - 0)**2
        elif usesHeights :
          f = lambda (h,b) : b * fctr
        else :
          f = lambda b : b * fctr

        a0 = sum([f(x) for x in posteriorParts[k]])
        c0 += a0
        z0 += a0

  # Not used anymore, save memory now
  del posteriorParts ; posteriorParts = None
  
  # Tuck in the constant
  ee = ("%.15g " % c0) + ee

  if z0 >= limit :
    return (None, 0.0)

  if initMethod == "opt":
    if hsOnly :
      print targets
      exec ("def ftargets():\n  return " + targets + ')') in globals()
      _setTreeHeightsForTargets(tree, ftargets)
    else : 
      _setTreeHeights(tree, optBranches)
  elif norm :
    for n in tree.all_ids() :
      tree.node(n).data.branchlength *= fctr
      
  # Get code which transforms the heights encoding to branch lengths
  # A descendant height is specified as a fraction in [0,1] of its ancestor
  # height (but leading number is the root height).

  ba,minRoot,htox = _treeBranchAssignmentExprs(tree, treeParts, fctr,
                                               nodesMinHeight = nodesMinHeight,
                                               withDerivative = withDerivative,
                                               withInit = True, withHeights = usesHeights)
    
  # Define the distance function on the fly.
  cod = ("def v1score(x):\n  " + "\n  ".join(ba) + "\n" + pee + "\n  return " +
         (('(' + ee + ", array([(" + dee + ") for k in range(nDer)]) )")
          if withDerivative else ee))
  exec cod in globals()

  if verbose: print cod

  # Number of variables (heights)
  nx = len(tree.get_terminals())-1

  # Function for obtaining the encoding of starting target tree
  xcod = "def htoxs():\n  x = [0]*%d\n  " % nx + "\n  ".join(htox) \
         + "\n  return x"
  exec xcod in globals()

  # Function to obtain branch lengths from encoding
  codb = "def code2branches(x):\n  " + "\n  ".join(ba) + "\n  " + \
       "return (" + ",".join(['(%d,b%d)' % ((treeParts[k][0][1],)*2) for k in treeParts]) + ")" 
  exec codb in globals()
  if verbose: print codb

  if verbose :
    print "@@",nx, minRoot, treeHeight(tree) * fctr
    print cod

  maxfun = 15000

  # small protection against disasters
  while True:
    if initMethod != "random" :
      x0 = htoxs()
    else :
      x0 = [1 if norm else treeHeight(tree)] + \
           [random.random() for k in range(nx-1)]

    initialVal = v1score(x0)     # assert x0[0] >= minRoot
    if withDerivative :
      initialVal = initialVal[0]
    if 1 :  
      sol = scipy.optimize.fmin_l_bfgs_b(v1score, x0,
                                        approx_grad=0 if withDerivative else 1,
                                        bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                        factr = factr,
                                        iprint=-1, maxfun=maxfun)
      if warnings and sol[2]['warnflag'] != 0 :
        print "WARNING:", sol[2]['task']
    else :
      sol = scipy.optimize.fmin_tnc(v1score, x0,
                                    approx_grad=0 if withDerivative else 1,
                                    bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                    maxfun=maxfun,
                                    messages=0)
      assert sol[2] == 1

    finaleVal = v1score(sol[0])
    if withDerivative :
      finaleVal = finaleVal[0]

    if finaleVal <= initialVal :
      break

    # try again from a random spot
    initMethod = "random"
    factr /= 10
    if factr < 1e6 :
      # failed, leave as is
      sol = htox()
      finaleVal = v1score(sol[0])
      if withDerivative :
        finaleVal = finaleVal[0]
      break

  brs = code2branches(sol[0])
  for nn,br in brs:
    # numerical instability : don't permit negative branches
    tree.node(nn).data.branchlength = max(br/fctr, 0)

  val = finaleVal
  ## if withDerivative :
  ##   val = val[0]

  # if norm under BRANCH_SCORE_2, there is no way to scale back
  return (tree, val/fctr if (norm and method != BRANCH_SCORE_2) else val) + \
         ((v1score, htoxs, code2branches, sol[0], fctr) if internals else tuple())






















if 0 :
# keep a little while for reference



  def minPosteriorDistanceTree(tree, trees, limit = scipy.inf, norm = True,
                               nodesMinHeight = None, withDerivative = False,
                               withInit = False, factr=10000000.0, warnings = True) :
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
      nn = treeParts[k][0][1]
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

        ee += "+ %.15f * b%d**2 + %.15f * b%d" % (a1 + len(br), nn, -2*sum(br), nn)

        if withDerivative:
          dee += ("+ %.15f * 2 * b%d * d_b%d_x[k] + %.15f * d_b%d_x[k]" %
                  (a1 + len(br), nn, nn, -2*sum(br), nn))

        # The constant term contribution 
        c0 += sum([x**2 for x in br])
      else :
        # A tree clade not appearing in posterior: contributes the full branch
        # length for each posterior tree.

        ee += "+(%d * b%d**2)" % (len(trees), nn)

        if withDerivative:
          dee += "+(%d * 2 * b%d * d_b%d_x[k])" % (len(trees), nn, nn)

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

    ba,minRoot,htox = _treeBranchAssignmentExprs(tree, treeParts, fctr,
                                                 nodesMinHeight = nodesMinHeight,
                                                 withDerivative = withDerivative,
                                                 withInit = True)

    # Define the posterior distance function on the fly.
    cod = ("def f(x):\n  " + "\n  ".join(ba) + "\n  return " +
           (('(' + ee + ", array([(" + dee + ") for k in range(nDer)]) )")
            if withDerivative else ee))
    exec cod

    # Number of variables (heights)
    nx = len(tree.get_terminals())-1

    xcod = "def htoxs():\n  x = [0]*%d\n  " % nx + "\n  ".join(htox) \
           + "\n  return x"
    exec xcod 

    if verbose :
      print "@@",nx, minRoot, treeHeight(tree) * fctr
      print cod

    maxfun = 15000
    if 1 :
      if withInit :
        x0 = htoxs()
        if norm:
          x0[0] = 1
      else :
        x0 = [1 if norm else treeHeight(tree)] + [random.random() for k in
                                                  range(nx-1)]

      assert x0[0] >= minRoot
      zz = scipy.optimize.fmin_l_bfgs_b(f, x0,
                                        approx_grad=0 if withDerivative else 1,
                                        bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                        factr = factr,
                                        iprint=-1, maxfun=maxfun)
      if warnings and zz[2]['warnflag'] != 0 :
        print "WARNING:", zz[2]['task']
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
    codb = "def code2branches(x):\n  " + "\n  ".join(ba) + "\n  " + \
         "return (" + ",".join(['(%d,b%d)' % ((treeParts[k][0][1],)*2) for k in treeParts]) + ")" 
    ## cod = "def b(x):\n  " + "\n  ".join(ba) + "\n  " + \
    ##      "return (" + ",".join(['b%d' % k for k in range(len(treeParts))]) + ")" 
    exec codb
    if verbose: print codb

    # Do not change tree passed as argument. Copy tree and set branch lengths of
    # the copy.
    ss = copy.deepcopy(tree)

    brs = code2branches(zz[0])
    for nn,br in brs:
      ss.node(nn).data.branchlength = br/fctr

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

  from cchelp import sumNonIntersect as _sumNonIntersect, sumNonIntersectDer as _sumNonIntersectDer

  def minPosteriorRADistanceTree(tree, trees, limit = scipy.inf, norm = True,
                                 nodesMinHeight = None, withDerivative = False,
                                 withInit = True, factr=10000000.0, warnings= True) :
    """ Find a branch length assignment for tree which minimizes the total
    distance to the set of trees undex Rooted Agreement Distance.

    limit is an upper bound (presumably from a prior call here with another tree).
    If the distance is known to be larger, the optimization for this tree can be
    skipped.
    """

    # tip/node lower bounds unimplemented
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
    pee = ""

    for r,k in enumerate(treeParts) :
      nn = treeParts[k][0][1]
      if k in posteriorParts :
        # A tree clade which is in some posterior trees

        # Branchs from posterior for this clade
        br = [(h*fctr,b*fctr) for h,b in posteriorParts[k]]

        # Number of posterior trees without the clade 
        a1 = (len(trees) - len(posteriorParts[k]))

        assert len(trees) == a1 + len(br)

        dt = '[' + ','.join(["(%.15g,%.15g)" % (h,h+b) for h,b in br]) + ']'
        if withDerivative :
          pee += "  val%d,dvdl%d,dvdh%d = _sumNonIntersectDer(bs%d,bs%d+b%d,%s)\n" % ((nn,)*6 + (dt,))
          # make low the sum of low+high, save operation in loops
          pee += "  dvdl%d += dvdh%d\n" % (nn,nn)
          # fold constant into high, same reason
          pee += "  dvdh%d += %d\n" % (nn,a1)

          ee += "+ val%d" % (nn,)

          # brevity: dx == dx_k
          #
          # dV/dx = dV(l,u)/dx = dV(l)/dx + dV(u)/dx = dV(l)/dl * dl/dx + dV(u)/du * du/dx

          if tree.node(nn).data.taxon :
            dee += "+(dvdh%d * d_b%d_x[k])" % (nn,nn)     
          else :
            dee += "+(dvdl%d * d_h%d_x[k] + dvdh%d * d_b%d_x[k])" % (nn,nn,nn,nn)     

          ## if tree.node(nn).data.taxon :
          ##   dee += "+((dvdh%d+%d) * d_b%d_x[k])" % (nn,a1,nn)     
          ## else :
          ##   dee += "+((dvdl%d+dvdh%d) * d_h%d_x[k] + (dvdh%d + %d) * d_b%d_x[k])" \
          ##          % (nn,nn,nn,nn,a1,nn)     
        else :
          ee += "+ _sumNonIntersect(bs%d,bs%d+b%d,%s)" % (nn,nn,nn,dt)

        if a1 > 0 :
          ee += "+ %d * b%d" % (a1,nn)

      else :
        # A tree clade not appearing in posterior: contributes the full branch
        # length for each posterior tree.

        ee += "+(%d * b%d)" % (len(trees), nn)

        if withDerivative:
          dee += "+(%d * d_b%d_x[k])" % (len(trees), nn)

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
    # height (but root height comes first).
    ba,minRoot,htox = _treeBranchAssignmentExprs(tree, treeParts, fctr,
                                                 nodesMinHeight = nodesMinHeight,
                                                 withDerivative = withDerivative,
                                                 withInit = True, withHeights = True)

    # Define the posterior distance function on the fly.
    cod = ("def rascore(x):\n  " + "\n  ".join(ba) + "\n" + pee + "\n  return " +
           (('(' + ee + ", array([(" + dee + ") for k in range(nDer)]) )")
            if withDerivative else ee))
    exec cod in globals()
    if verbose: print cod

    # Number of variables (heights)
    nx = len(tree.get_terminals())-1

    xcod = "def htox():\n  x = [0]*%d\n  " % nx + "\n  ".join(htox) \
           + "\n  return x"
    exec xcod

    # Function to get the branch lengths from optimized heights
    codb = "def code2branches(x):\n  " + "\n  ".join(ba) + "\n  " + \
         "return (" + ",".join(['(%d,b%d)' % ((treeParts[k][0][1],)*2) for k in treeParts]) + ")" 

    exec codb in globals()
    if verbose: print codb

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


    # Do not change tree passed in arguments. Set optimized branch lengths on
    # a copy.
    ss = copy.deepcopy(tree)

    brs = code2branches(zz[0])
    for nn,br in brs:
      ss.node(nn).data.branchlength = br/fctr

    val = finaleVal
    if withDerivative :
      val = val[0]
    return (ss, val/fctr if norm else val)




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
    pee = ""

    for r,k in enumerate(treeParts) :
      nn = treeParts[k][0][1]
      if k in posteriorParts :
        # A tree clade which is in some posterior trees

        # Branchs from posterior for this clade
        br = [b*fctr for b in posteriorParts[k]]

        # Number of posterior trees without the clade 
        a1 = (len(trees) - len(posteriorParts[k]))

        assert len(trees) == a1 + len(br)
        if withDerivative :
          pee += "  ab%d,abd%d = _absDiffBranchDer(b%d,%s)\n" % (nn,nn,nn,_prepare(br))
          pee += "  abd%d += %d\n" % (nn,a1)

          ee += "+ ab%d" % (nn,)
        else :
          ee += "+ _absDiffBranch(b%d,%s)" % (nn,_prepare(br))

        ee += "+ %d * b%d" % (a1,nn)

        if withDerivative:
          dee += "+(abd%d * d_b%d_x[k])" % (nn,nn)        
      else :
        # A tree clade not appearing in posterior: contributes the full branch
        # length for each posterior tree.

        ee += "+(%d * b%d)" % (len(trees), nn)

        if withDerivative:
          dee += "+(%d * d_b%d_x[k])" % (len(trees), nn)

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

    # Define the posterior distance function on the fly.
    cod = ("def v1score(x):\n  " + "\n  ".join(ba) + "\n" + pee + "\n  return " +
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

    # Function to get the branch lengths from optimized heights
    codb = "def code2branches(x):\n  " + "\n  ".join(ba) + "\n  " + \
         "return (" + ",".join(['(%d,b%d)' % ((treeParts[k][0][1],)*2) for k in treeParts]) + ")" 
    exec codb in globals()
    if verbose: print cod

    if verbose :
      print "@@",nx, minRoot, treeHeight(tree) * fctr
      print cod

    if 0 :
      cod8 = "def tv1score(tree):\n"
      for r,k in enumerate(treeParts) :
        cod8 += "  b%d = tree.node(%d).data.branchlength\n" % ((treeParts[k][0][1],)*2)
      cod8 += "  return " +  ee
      exec cod8

      nhs = nodeHeights(tree, allTipsZero = False)

      cod9 = "def tv1scoreh(hs):\n"
      hs = getPreOrder(tree, includeTaxa=False)
      for b in getPreOrder(tree, includeTaxa=True)[1:] :
        nd = tree.node(b)
        if not nd.succ:
          nh = "%.15g" % nhs[b]
        else :
          nh = "hs[%d]" % hs.index(b)
        cod9 += "  b%d = hs[%d] - %s\n" % (b, hs.index(tree.node(b).prev), nh)
      cod9 += "  return " +  ee
      exec cod9

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
                                   approx_grad=1,
                                   bounds = [[minRoot,None]] + [[0,1]]*(nx-1),
                                   maxfun=maxfun,
                                   messages=0)
      assert zz[2] == 1


    # Do not change tree passed as argument. Copy tree and set branch lengths of
    # the copy.
    ss = copy.deepcopy(tree)

    brs = code2branches(zz[0])
    for nn,br in brs:
      ss.node(nn).data.branchlength = br/fctr

    val = finaleVal
    if withDerivative :
      val = val[0]
    return (ss, val/fctr if norm else val)


  def minPosteriorHSDistanceTree(tree, trees, limit = scipy.inf, norm = True,
                                 nodesMinHeight = None, withDerivative = False,
                                 withInit = True, factr=10000000.0,
                                 warnings = True) :
    """ Find a branch length assignment for tree which minimizes the total
    distance to the set of trees.

    limit is an upper bound (presumably from a prior call here with another tree).
    If the distance is known to be larger, the optimization for this tree can be
    skipped.
    """

    #assert not withDerivative
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
    pee = ""

    for r,k in enumerate(treeParts) :
      nn = treeParts[k][0][1]
      if k in posteriorParts :
        # A tree clade which is in some posterior trees

        # Branchs from posterior for this clade
        br = [(h*fctr,b*fctr) for h, b in posteriorParts[k]]

        # Number of posterior trees without the clade 
        a1 = (len(trees) - len(posteriorParts[k]))

        assert len(trees) == a1 + len(br)

        if not tree.node(nn).data.taxon  :

          if withDerivative :
            pee += "  ab%d,abd%d = _absDiffBranchDer(b%d,%s)\n" % (nn,nn,nn,_prepare(br))
            pee += "  abd%d += %d\n" % (nn,a1)

            ee += "+ ab%d" % (nn,)
          else :
            ee += "+ _absDiffBranch(bs%d,%s)" % (nn,_prepare([x[0] for x in br]))

          ee += "+ %d * b%d" % (a1,nn)

          if withDerivative:
            dee += "+(abd%d * d_b%d_x[k])" % (nn,nn)        

      else :
        # A tree clade not appearing in posterior: contributes the full branch
        # length for each posterior tree.

        ee += "+(%d * b%d)" % (len(trees), nn)

        if withDerivative:
          dee += "+(%d * d_b%d_x[k])" % (len(trees), nn)

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

    # Define the posterior distance function on the fly.
    cod = ("def v1score(x):\n  " + "\n  ".join(ba) + "\n" + pee + "\n  return " +
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

    #global code2branches
    # Function to get the branch lengths from optimized heights
    codb = "def code2branches(x):\n  " + "\n  ".join(ba) + "\n  " + \
         "return (" + ",".join(['(%d,b%d)' % ((treeParts[k][0][1],)*2) for k in treeParts]) + ")" 
    exec codb in globals()
    if verbose: print cod
    #code2branches = code2branchesa
    if verbose :
      print "@@",nx, minRoot, treeHeight(tree) * fctr
      print cod

    if 0 :
      cod8 = "def tv1score(tree):\n"
      for r,k in enumerate(treeParts) :
        cod8 += "  b%d = tree.node(%d).data.branchlength\n" % ((treeParts[k][0][1],)*2)
      cod8 += "  return " +  ee
      exec cod8

      nhs = nodeHeights(tree, allTipsZero = False)

      cod9 = "def tv1scoreh(hs):\n"
      hs = _getNodeIDsDescendingHeight(tree, tree.root, 0)
      for b in _getNodeIDsDescendingHeight(tree, tree.root, 1)[1:] :
        nd = tree.node(b)
        if not nd.succ:
          nh = "%.15g" % nhs[b]
        else :
          nh = "hs[%d]" % hs.index(b)
        cod9 += "  b%d = hs[%d] - %s\n" % (b, hs.index(tree.node(b).prev), nh)
      cod9 += "  return " +  ee
      exec cod9

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
        #pdb.set_trace()
        #global mcalls, dcalls
        #mcalls, dcalls = 0,0

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
    for nn,br in brs:
      ss.node(nn).data.branchlength = br/fctr
    ## for r,k in enumerate(treeParts) :
    ##   ss.node(treeParts[k][0][1]).data.branchlength = brs[r]/fctr

    val = finaleVal
    if withDerivative :
      val = val[0]
    return (ss, val/fctr if norm else val)
