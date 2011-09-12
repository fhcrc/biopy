## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

from __future__ import division

import sys
       
from numpy import mean, median, corrcoef
import pylab

from collections import Counter, namedtuple, defaultdict

from biopy.genericutils import fileFromName

from biopy import INexus, beastLogHelper, demographic
from biopy.treeutils import toNewick, countNexusTrees, getTreeClades, nodeHeights


__all__ = ['drawTree', 'getSpacing', 'getTaxaOrder' , "descendantMean",
           "descendantWeightedMean", "descendantBetween"]

class descendantMean(object) :
  def __init__(self, descendants) :
    if isinstance(descendants, (float,int,long)) :
      self.xcenter = descendants
    else :
      self.xcenter = mean([x.xcenter for x in descendants])

  def center(self) :
    return self.xcenter
  
class descendantWeightedMean(object) :
  def __init__(self, descendants) :
    if isinstance(descendants, (float,int,long)) :
      self.xcenter = descendants
      self.count = 1
    else :
      self.count = sum([x.count for x in descendants])
      self.xcenter = sum([x.xcenter * x.count for x in descendants])/self.count

  def center(self) :
    return self.xcenter

class descendantBetween(object) :
  def __init__(self, descendants) :
    if isinstance(descendants, (float,int,long)) :
      self.xcenter = descendants
      self.right = self.left = self.xcenter
    else :
      self.left = min([x.left for x in descendants])
      self.right = max([x.right for x in descendants])
      self.xcenter = (descendants[0].right + descendants[1].left) / 2

  def center(self) :
    return self.xcenter

def drawTree(tree, nid, cladesDict = None, positioning = descendantMean,
             fill = True, color="lime", alpha = 0.05) :
  node = tree.node(nid)
  d = node.data.demographic
  if not node.succ:
    p = d.population(0)
    br = node.data.branchlength
    p1 = d.population(br)
    center = node.data.x
    
    return (positioning(center), 0, br,center + p/2, p, p1, 
            frozenset([node.data.taxon]) if cladesDict is not None else None)

  pr = []
  for si in node.succ :
    pr.append(drawTree(tree, si, cladesDict, positioning, fill, color, alpha))

  pr = sorted(pr, key = lambda x : x[0].center()) # sort by center
  centerc = positioning([x[0] for x in pr])
  pcenter = centerc.center()
                    
  clade = pr[0][6] | pr[1][6] if cladesDict is not None else None

  if clade is None or clade in cladesDict:
    for (center,h,hpar,xright,p,p1,c),orient in zip(pr,[0,1]):
      if c is not None and c not in cladesDict[clade]:
        continue

      if orient == 0:
        x = xright, xright - p, pcenter-p1, pcenter
      else:
        x = xright, xright - p, pcenter, pcenter + p1 

      y = h, h, hpar, hpar
      if fill :
        pylab.fill(x,y, color=color, alpha = alpha)
      else :
        x = x + (x[0],) ;  y = y + (y[0],)
        pylab.plot(x,y, color=color)

  p = sum([x[5] for x in pr])
  p1 = d.population(d.naturalLimit())
  h = pr[0][2]
  hpar = h + d.naturalLimit()
  xright = pcenter + pr[1][5]
  
  if nid == tree.root and \
     (cladesDict is None or all([x[6] in cladesDict[clade] for x in pr])):
    x = xright, xright - p, pcenter - p1/2, pcenter + p1/2
    y = h, h, hpar, hpar

    if fill :
      pylab.fill(x,y, color=color, alpha=alpha)
    else :
      x = x + (x[0],) ;  y = y + (y[0],)
      pylab.plot(x,y, color=color)

  if nid == tree.root :
    return hpar
  
  return (centerc, h, hpar, xright, p, p1, clade)


_Info = namedtuple('Info', 'width widthatend maxspace')

def treeMinSpace(tree, nid = -1) :
  if nid == -1 :
    nid = tree.root
  node = tree.node(nid)
  d = node.data.demographic
  if not node.succ:
    return _Info(d.population(0), d.population(node.data.branchlength), 0)

  l,r = [treeMinSpace(tree, x) for x in node.succ]
  spc = 2*max(l.widthatend - l.width, r.widthatend - r.width, 0)

  spc = max(spc, l.maxspace,r.maxspace)
  if nid == tree.root:
    return spc
  return _Info(spc +  l.width + r.width, d.population(d.naturalLimit()),spc)

def getSpacing(trees, additional = 0) :
  wd = max([mean(sorted([tree.node(x).data.demographic.population(0) for x in
                  tree.get_terminals()])[-2:]) for tree in trees])
  spc = max([treeMinSpace(tree) for tree in trees])

  assert additional >= 0
  return max((1+additional)*wd, spc)


def taxaDistance(tree, nid, distance) :
  node = tree.node(nid)
  if not node.succ:
    return (1, [node.data.taxon])
  else:
    l,r = [taxaDistance(tree, x, distance) for x in node.succ]
    n = l[0] + r[0]
    for xl in l[1] :
      for xr in r[1] :
        key = (xl,xr) if xl < xr else (xr,xl)
        distance[key] = n-1
    return (n, l[1] + r[1])


def getAllSingleFlipTaxaOrders(tree, nid) :
  node = tree.node(nid)
  if not node.succ:
    return [[node.data.taxon]]
  else:
    l,r = [getAllSingleFlipTaxaOrders(tree, x) for x in node.succ]
    l0,r0 = l[0],r[0]

    all = []
    all.append(l0+r0)
    all.append([r0+l0, node])
    for x,i in l[1:] :
      all.append([x + r0, i])
    for x,i in r[1:] :
      all.append([l0 + x, i])
    return all

def getCR(oo, dis) :
  x,y = [],[]
  n = len(oo)
  for i in xrange(n) :
    ti = oo[i]
    for j in xrange(i+1,n) :
      x.append(j-i)
      tj = oo[j]
      key = (ti,tj) if ti < tj else (tj,ti)
      y.append(dis[key])
  cr = corrcoef(x,y)[0,1]
  return cr


def getTaxaOrder(trees, refTree = None, reportTopologies = False,
                    nSample = 1, progress = False) :
  
  if progress: print >> sys.stderr, "getting topologies...",

  dtops = Counter(toNewick(tree, topologyOnly=1) for tree in trees)

  averagedPosteriorDistances = defaultdict(lambda : 0)

  if progress: print >> sys.stderr, ("taxa distance matrix (%d tops)..." % len(dtops)),
  
  for tx in dtops:
    tree = INexus.Tree(tx)
    distances = dict()
    taxaDistance(tree, tree.root, distances)
    for k in distances:
      averagedPosteriorDistances[k] += distances[k] * dtops[tx]/len(trees)
    
  if refTree is None:
    sdtops = sorted(dtops, key = lambda x : dtops[x])
    refTrees = [INexus.Tree(t) for t in sdtops[-nSample:]]
  else :
    refTrees = [refTree]

  overallMaxCR = -1
  
  for tree in refTrees :
    allOrders = getAllSingleFlipTaxaOrders(tree, tree.root)
    maximizingOrder = allOrders[0]
    treeMaxCR = getCR(maximizingOrder, averagedPosteriorDistances)

    if progress: print >> sys.stderr, "optimizing... (cr",treeMaxCR,")",
    nTries = 0

    while True:
      nTries += 1
      mnode = None

      for order,node in allOrders[1:] :
        cr = getCR(order, averagedPosteriorDistances)
        
        if cr > treeMaxCR:
          treeMaxCR = cr
          maximizingOrder = order
          mnode = node
          
      if progress: print >> sys.stderr, treeMaxCR,

      if mnode is not None:
        mnode.succ = [mnode.succ[1], mnode.succ[0]]
        allOrders = getAllSingleFlipTaxaOrders(tree, tree.root)
      else :
        break

    if treeMaxCR > overallMaxCR :
      overallMaxCR = treeMaxCR
      overallMaximizingOrder = maximizingOrder
      mtree = tree
      
    if progress: print >> sys.stderr, ("%d tries" % nTries)

  if progress: print >> sys.stderr, "done"

  if reportTopologies :
    dtops = dict([(top,[]) for top in dtops])
    for k,top in enumerate(toNewick(tree, topologyOnly=1) for tree in trees) :
      dtops[top].append(k)
      
    return (overallMaximizingOrder, mtree, dtops)
  
  return (overallMaximizingOrder, mtree)


#_Info1 = namedtuple('Info1', 'left right lleft lright rleft rright lchild')

# left right: x-position of leftest and rightest tips in clade
# lmid rmid:  x-position of left and right descentands to aim for
# lh rh:      y-position of left and right descentands to aim for
# lchild:     node id of the left child

_Info1 = namedtuple('Info1', 'left right lmid rmid lh rh lchild')

def _setLeftRight(tree, nid, nh) :
  node = tree.node(nid)
  if not node.succ:
    x = node.data.x
    h = nh[nid]
    node.data.cladeInfo = _Info1(x, x, x, x, h, h, None)
  else :
    chi = [_setLeftRight(tree, si, nh) for si in node.succ]
    r = max([x.right for x in chi])
    l = min([x.left for x in chi])
    m = [x.left + x.right for x in chi]

    ilchild = 0 if m[0] < m[1] else 1
    
    chil = chi[ilchild]
    chir = chi[1-ilchild]
    node.data.cladeInfo = _Info1(l, r,
                                 (chil.left + chil.right)/2,
                                 (chir.left + chir.right)/2,
                                 (chil.lh + chil.rh)/2,
                                 (chir.lh + chir.rh)/2,
                                 node.succ[ilchild])
    ## node.data.cladeInfo = _Info1(l, r,
    ##                              (chil.left + chil.right)/2,
    ##                              (chir.left + chir.right)/2,
    ##                              max(chil.lh , chil.rh),
    ##                              max(chir.lh , chir.rh),
    ##                              node.succ[ilchild])

  return node.data.cladeInfo

def _drawTreeOnly(tree, nid, xnode, nh, color, alpha, keepPositions) :
  node = tree.node(nid)
  if not node.succ:
    return
  cInfo = node.data.cladeInfo
  myh = nh[node.id]
  for si in node.succ :
    chh =  nh[si]
    if si == cInfo.lchild :
      a = (xnode - cInfo.lmid) / (myh - cInfo.lh)
      dx = a * (myh - chh)
      xchild = xnode - dx

    else :
      a = (cInfo.rmid - xnode) / (myh - cInfo.rh)
      dx = a * (myh - chh)
      xchild = xnode + dx
      
    # print node.id, xnode, si, si == cInfo.lchild, xchild, myh, chh, cInfo
    if keepPositions :
      node.data.x = xnode
    pylab.plot([xnode, xchild], [myh, chh], color = color, alpha = alpha) 
    _drawTreeOnly(tree, si, xchild, nh, color, alpha, keepPositions)


def drawTreeOnly(tree, color="green", alpha = 0.05, allTipsZero = True, xroot =
                 None, keepPositions = False) :
  # A mapping from node id to (positive) node height ( with !allTipsZero, tips
  # may have > 0 height, at least one tip has 0 height). 
  nh = nodeHeights(tree, allTipsZero = allTipsZero)

  # set auxiliary info per node
  rl = _setLeftRight(tree, tree.root, nh)

  # position root in the middle
  if xroot is None :
    xroot = (rl.right + rl.left)/2
  
  _drawTreeOnly(tree, tree.root, xroot, nh, color= color, alpha = alpha,
                keepPositions = keepPositions)
  
  for i in tree.all_ids() :
    del tree.node(i).data.cladeInfo


## def getTaxa(tree, nid, rand = False) :
##   node = tree.node(nid)
##   if not node.succ:
##     return [node.data.taxon]
##   else:
##     l,r = [getTaxa(tree, x, rand) for x in node.succ]
##     if rand :
##       if random.random() < .5:
##         l,r = r,l
##     return l + r
## def getCRold(oo, dis) :
##   x,y = [],[]
##   for k in dis:
##     tx1,tx2 = k
##     x.append(abs(oo.index(tx1) - oo.index(tx2)))
##     y.append(dis[k])

##   cr = corrcoef(x,y)[0,1]
##   return cr

## def getTaxaOrderOld(trees, refTree = None, reportTopologies = False, nSample = 1) :
##   print >> sys.stderr, "getting tops...",
##   tops = [toNewick(tree, topologyOnly=1) for tree in trees]
##   dtops = dict([(x,0) for x in tops])
##   for t in tops :
##     dtops[t] += 1
##   dis = dict()

##   print >> sys.stderr, ("taxa distance matrix (%d tops)..." % len(dtops)),
  
##   for tx in dtops:
##     tree = INexus.Tree(tx)
##     d = dict()
##     taxaDistance(tree, tree.root, d)
##     for k in d:
##       if k not in dis:
##         dis[k] = 0.0
##       dis[k] += d[k] * dtops[tx]/len(trees)
    
##   if refTree is None:
##     sdtops = sorted(dtops, key = lambda x : dtops[x])
##     refTrees = [INexus.Tree(t) for t in sdtops[-nSample:]]
##   else :
##     refTrees = [refTree]

##   mmcr = -1
  
##   for tree in refTrees :
##     moo = getTaxa(tree, tree.root)
##     mcr = getCRold(moo, dis)

##     print >> sys.stderr, "optimizing...",mcr,
##     nt = 0

##     while True:
##       nt += 1
##       mnode = None

##       for nid in tree.all_ids() :
##         node = tree.node(nid)
##         if node.succ:
##           node.succ = [node.succ[1], node.succ[0]]
##           oo = getTaxa(tree, tree.root)
##           cr = getCRold(oo, dis)

##           if cr > mcr:
##             mcr = cr
##             moo = oo
##             mnode = node
##           node.succ = [node.succ[1], node.succ[0]]
##       print >> sys.stderr, mcr,

##       if mnode is not None:
##         mnode.succ = [mnode.succ[1], mnode.succ[0]]
##       else :
##         break

##     if mcr > mmcr :
##       mmcr = mcr
##       mmoo = moo
##       mtree = tree
      
##     print >> sys.stderr, ("%d tries" % nt)

##   print >> sys.stderr, "done"

##   if reportTopologies :
##     for top in dtops:
##       dtops[top] = []
##     for k,top in enumerate(tops):
##       dtops[top].append(k)
      
##     return (mmoo, mtree, dtops)
  
##   return (mmoo, mtree)


