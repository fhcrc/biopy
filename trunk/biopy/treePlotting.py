## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

from __future__ import division

import sys, random
       
from numpy import mean, median, corrcoef, cumsum
import pylab

from collections import Counter, namedtuple, defaultdict

from genericutils import fileFromName

from biopy import INexus, beastLogHelper, demographic
from treeutils import toNewick, countNexusTrees, getTreeClades, nodeHeights, getPostOrder
from combinatorics import allPairs

__all__ = ['drawTree', 'getSpacing',
           "DescendantMean", "DescendantWeightedMean",
           "DescendantBetween", "DescendantBalanced", 
           'getTaxaOrder', 'getGTorder', "setSTspacing"]

class DescendantMean(object) :
  def __init__(self, node, descendants) :
    if descendants is not None :
      assert len(descendants) == 2
      self.xcenter = mean([x.xcenter for x in descendants])
    else :
      self.xcenter = node.data.x
    
  def center(self) :
    return self.xcenter

  def setOffset(self, offset) :
    self.xcenter += offset
    
  def widenBy(self, offset) :
    return self.xcenter + offset/2
  
class DescendantWeightedMean(object) :
  def __init__(self, node, descendants) :
    if descendants is not None :
      self.count = sum([x.count for x in descendants])
      self.xcenter = sum([x.xcenter * x.count for x in descendants])/self.count
      self.rcount = descendants[-1].count
    else :
      self.xcenter = node.data.x
      self.count = 1
      self.rcount = 0
      
  def center(self) :
    return self.xcenter

  def widenBy(self, extra) :
    return self.xcenter + extra*self.rcount/self.count

  def setOffset(self, offset) :
    self.xcenter += offset
    
class DescendantBetween(object) :
  def __init__(self, node, descendants) :
    if descendants is not None :
      self.left = min([x.left for x in descendants])
      self.right = max([x.right for x in descendants])
      self.xcenter = (descendants[0].right + descendants[1].left) / 2
      assert len(descendants) == 2
    else :
      self.xcenter = node.data.x
      self.right = self.left = self.xcenter

  def center(self) :
    return self.xcenter

  def widenBy(self, offset) :
    return self.xcenter + offset/2

  def setOffset(self, offset) :
    self.xcenter += offset

class DescendantBalanced(object) :
  def __init__(self, node, descendants) :
    if descendants is not None :
      assert len(descendants) == 2
      self.xcenter = (2*(descendants[0].xcenter + descendants[1].xcenter) +
                      descendants[0].p[1] - descendants[1].p[1])/4.0
      #print descendants[0].xcenter, descendants[0].p[1],  descendants[1].xcenter,  descendants[1].p[1]
    else :
      self.xcenter = node.data.x
      
    self.p = [node.data.demographic.population(z)
              for z in (0,node.data.branchlength)]

  def center(self) :
    return self.xcenter

  def widenBy(self, offset) :
    return self.xcenter + offset/2

  def setOffset(self, offset) :
    self.xcenter += offset

_Info2 = namedtuple('Info2', 'center hnode hpar xright wnode wpar clade')

def _drawBranch(x, y, fill, generalPlotAttributes, splitPoints) :
  if fill is not None:
    pylab.fill(x,y, **fill)
    
  if generalPlotAttributes is not None:
    if splitPoints is True :
      x = x + (x[0],) ;  y = y + (y[0],)
      pylab.plot(x,y, **generalPlotAttributes)
    else :
      pylab.plot(x[1:3],y[1:3], **generalPlotAttributes)
      pylab.plot([x[0],x[-1]],[y[0],y[-1]], **generalPlotAttributes)
      if isinstance(splitPoints, dict) :
        pylab.plot(x[0:2],y[0:2], **splitPoints)
        pylab.plot(x[2:3],y[2:3], **splitPoints)

def _drawTree(tree, nid, cladesDict, positioning, fill, generalPlotAttributes,
              splitPoints, keepAux) :
  node = tree.node(nid)
  d = node.data.demographic
  if not node.succ:
    p = d.population(0)
    br = node.data.branchlength
    p1 = d.population(br)
    center = node.data.x

    aux = _Info2(positioning(node, None), 0, br, center + p/2, p, p1, 
                 frozenset([node.data.taxon]) if cladesDict is not None else None)
    if keepAux :
      node.data.plotAux = aux
    return aux

  pr = [_drawTree(tree, si, cladesDict, positioning, fill, generalPlotAttributes,
                  splitPoints, keepAux)for si in node.succ]

  pr = sorted(pr, key = lambda x : x.center.center()) # sort by center
  centerc = positioning(node, [x.center for x in pr])
  pcenter = centerc.center()
                    
  clade = pr[0].clade | pr[1].clade if cladesDict is not None else None

  if clade is None or clade in cladesDict:
    for (center,h,hpar,xright,p,p1,c),orient in zip(pr,[0,1]):
      if c is not None and c not in cladesDict[clade]:
        continue

      if orient == 0:
        x = xright, xright - p, pcenter-p1, pcenter
      else:
        x = xright, xright - p, pcenter, pcenter + p1 

      y = h, h, hpar, hpar

      _drawBranch(x, y, fill, generalPlotAttributes, splitPoints)      

  pConst = d.naturalLimit() is None
  if pConst :
    p = d.population(0)
  else :
    p = sum([x.wpar for x in pr])

  br = d.naturalLimit() or node.data.branchlength

  p1 = d.population(br)
  h = pr[0].hpar
  hpar = h + br
  if pConst :
    #xright = pcenter + p/2  # + (pr[1].wnode - pr[0].wnode)/2 + p/2
    xright = pcenter + (p * pr[1].wnode)/(pr[1].wnode + pr[0].wnode)
  else :
    xright = pcenter + pr[1].wpar
  
  if nid == tree.root and \
     (cladesDict is None or all([x.clade in cladesDict[clade] for x in pr])):
    if pConst :
      if h == hpar :
        hpar = 1.1 * h
      pcenter = xright - p/2

    x = xright, xright - p, pcenter - p1/2, pcenter + p1/2
    y = h, h, hpar, hpar
    
    _drawBranch(x, y, fill, generalPlotAttributes, splitPoints)      

  aux = _Info2(centerc, h, hpar, xright, p, p1, clade)
  if keepAux :
    node.data.plotAux = aux
    
  if nid == tree.root :
    return hpar
  
  return aux

def drawTree(tree, nid = None, cladesDict = None, positioning = DescendantMean,
             fill = None, generalPlotAttributes = None, splitPoints = False,
             keepAux = False) :

  if fill is None and generalPlotAttributes is None :
    generalPlotAttributes = dict()

  if fill is not None and 'ec' not in fill :
    fill['ec'] = 'none'
    
  if generalPlotAttributes is not None and 'color' not in generalPlotAttributes :
    generalPlotAttributes['color'] = 'blue'
    
  if nid is None :
    nid = tree.root

  return _drawTree(tree, nid, cladesDict, positioning,
             fill, generalPlotAttributes, splitPoints, keepAux) 

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

  spc = max(spc, l.maxspace, r.maxspace)
  if nid == tree.root:
    return spc
  return _Info(spc +  l.width + r.width, d.population(d.naturalLimit()), spc)

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

def getGTorder(gtree, stree, nTries = 3, perTreeTries = 3) :
  """ get genes order to plot inside a species tree.

  preliminary version.
  """
  gtax = []
  gtx = defaultdict(lambda : [])
  for n in gtree.get_terminals() :
    gn = gtree.node(n)
    gtx[gn.data.snode.id].append(gn)
    gtax.append(gn)

  nh = nodeHeights(gtree)
  for x in gtree.all_ids() :
    gtree.node(x).data.ht = nh[x]
    
  ms, mp = getGTorderFixedStree(gtree, stree, gtax, gtx)
  sint = [stree.node(n) for n in stree.all_ids() if stree.node(n).succ]
  for nk in range(nTries) :
    cont = True
    while cont :
      random.shuffle(sint)
      cont = False
      for n in sint :
        for ll in range(perTreeTries) :
          sc = n.succ
          n.succ = [sc[1],sc[0]]
          tms, tmp = getGTorderFixedStree(gtree, stree, gtax, gtx)
          if tms < ms :
            ms, mp = tms, tmp
            cont = True
          else :
            n.succ = sc
  for g,o in zip(gtax,mp) :
    g.data.o = o
  return ms, mp, gtx

def getGTorderFixedStree(gtree, stree, gtax, gtx) :
  allsn = getPostOrder(stree, stree.root)
  # taxa in layout order
  stax = filter(lambda n : not n.succ, allsn)

  # all terminals contain data.snode 
  no = 0
  for i,n in enumerate(stax):
    for gn in gtx[n.id] :
      gn.data.grp = [i,i]
      gn.data.grpnds = [[gn],[gn]]
      gn.data.o = no
      no += 1
      gn.data.sz = 1

  gtpost = getPostOrder(gtree, gtree.root)
  for n in gtpost :
    if n.succ :
      sns = [gtree.node(x) for x in n.succ]
      l,r = zip(*[x.data.grp for x in sns])
      n.data.grp = min(l),max(r)
      lg,rg = zip(*[x.data.grpnds for x in sns])
      n.data.grpnds = reduce(lambda x,y : x+y, \
                             [y for x,y in zip(l,lg) if n.data.grp[0] == x]),\
                       reduce(lambda x,y : x+y, \
                              [y for x,y in zip(r,rg) if n.data.grp[1] == x])
      n.data.sz = sum([x.data.sz for x in sns])

  nint = filter(lambda n : len(n.succ), gtpost)
  def score(nint) :
    htot, tot = 0.0, 0
    for x in nint :
      l,r = [[z.data.o for z in u] for u in x.data.grpnds]
      dd = ((max(r) - min(l) + 1) - x.data.sz)
      if dd > 0 :
        tot += dd
        htot += dd * x.data.ht
    return tot, -htot

  # randomize order
  for kk in gtx:
    a = [n.data.o for n in gtx[kk]]
    random.shuffle(a)
    for n,i in zip(gtx[kk], a):
      n.data.o = i

  ms = score(nint)
  mp = [x.data.o for x in gtax]
  msLast = (sys.maxint, 0)
  while ms < msLast:
    msLast = ms
    for kk in gtx:
      gtxkk = gtx[kk]
      a = [n.data.o for n in gtxkk]
      for i0,i1 in allPairs(range(len(gtxkk))) :
        sw = [gtxkk[x].data.o for x in (i0,i1)]
        gtxkk[i1].data.o, gtxkk[i0].data.o = sw
        s = score(nint)
        if s < ms :
          ms = s
          mp = [x.data.o for x in gtax]
        else :
          gtxkk[i0].data.o, gtxkk[i1].data.o = sw
  return ms,mp

# left right: x-position of leftest and rightest tips in clade
# lmid rmid:  x-position of left/right descentands to aim for
# lh rh:      y-position of left/ right descentands to aim for
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

  return node.data.cladeInfo

# snode,x set in all gene taxa nodes
_Info3 = namedtuple('Info3', 'left right lmid rmid lh rh lchild snode')

# works for const(pop) and non const trees
def _tr1(v, snaux, nparaux, isLeft) :
  assert v <= snaux.xright
  c = nparaux.center.center()
  a = (snaux.xright - v) / snaux.wnode
  assert (0-1e-12) <= a <= (1+1e-12)
  if isLeft :
    vt = c - (c - (nparaux.xright - nparaux.wnode)) * a
  else :
    vt = c + (nparaux.xright - c) * (1-a)
  return vt
  
def _embSetLeftRight(gtree, nid, nh, stree, snh) :
  gnode = gtree.node(nid)
  if not gnode.succ:
    x = gnode.data.x
    h = nh[nid]
    gnode.data.cladeInfo = _Info3(x, x, x, x, h, h, None, gnode.data.snode)
  else :
    chi = [_embSetLeftRight(gtree, si, nh, stree, snh) for si in gnode.succ]
    myh = nh[nid]
    gnode.data.intermidiate = [[], []]
    for ns,x in enumerate(chi) :
      while x.snode.prev is not None and myh >= snh[x.snode.prev] :
        # need to move x info past species split
        sn = x.snode
        snaux = sn.data.plotAux
        npar = stree.node(sn.prev)
        nparaux = npar.data.plotAux
        isLeft = npar.succ[0] == sn.id
        u = [_tr1(v, snaux, nparaux, isLeft) for v in (x.left, x.right, x.lmid, x.rmid)]
        h = snh[x.snode.prev]
        u = u + [h,h,None,npar]
        u = _Info3(*u)
        im = gnode.data.intermidiate[ns]
        if len(im) == 0 :
          im.append(((x.left+x.right)/2, (x.lh + x.rh)/2))
          
        # they will come out in reverse order
        gnode.data.intermidiate[ns].insert(0,((u.left+u.right)/2, h))
        x = u
      chi[ns] = x
    
    r = max([x.right for x in chi])
    l = min([x.left for x in chi])
    m = [x.left + x.right for x in chi]

    ilchild = 0 if m[0] < m[1] else 1
    
    chil = chi[ilchild]
    chir = chi[1-ilchild]
    assert chil.snode == chir.snode
    gnode.data.cladeInfo = _Info3(l, r,
                                 (chil.left + chil.right)/2,
                                 (chir.left + chir.right)/2,
                                 (chil.lh + chil.rh)/2,
                                 (chir.lh + chir.rh)/2,
                                 gnode.succ[ilchild],
                                 chil.snode)

  return gnode.data.cladeInfo

def _drawTreeOnly(tree, nid, xnode, nh, plotLinesAttr, keepPositions, txt) :
  node = tree.node(nid)

  if txt:
    pylab.text(xnode, nh[node.id], str(nid), fontsize = 11, va='top', ha='center')
  
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
    pylab.plot([xnode, xchild], [myh, chh], **plotLinesAttr)

    #pylab.text(xnode, myh, str(nid), fontsize = 11, va='top', ha='center')
    
    _drawTreeOnly(tree, si, xchild, nh, plotLinesAttr, keepPositions, txt)

def _embDrawTreeOnly(tree, nid, xnodeOrig, nh, plotLinesAttr, keepPositions, coalAttr, txt) :
  node = tree.node(nid)

  if txt:
    pylab.text(xnodeOrig, nh[node.id], str(nid), fontsize = 11, va='top', ha='center')

  if coalAttr:
    pylab.plot([xnodeOrig], [nh[node.id]], '.', **coalAttr)
  
  if not node.succ:
    return
  data = node.data
  cInfo = data.cladeInfo
  myhOrig = nh[node.id]

  assert cInfo.lchild in node.succ
  
  for si,bp in zip(node.succ,data.intermidiate) :
    xnode = xnodeOrig
    myh = myhOrig
    chh = nh[si]
    if si == cInfo.lchild :
      x,h = cInfo.lmid, cInfo.lh
      if len(bp) :
        for ci in bp[:-1] :
          pylab.plot([xnode, ci[0]], [myh, ci[1]], **plotLinesAttr)
          xnode,myh = ci
        x,h = bp[-1]
        
      a = (xnode - x) / (myh - h)
      dx = a * (myh - chh)
      xchild = xnode - dx

    else :
      x,h = cInfo.rmid, cInfo.rh
      if len(bp) :
        for ci in bp[:-1] :
          pylab.plot([xnode, ci[0]], [myh, ci[1]], **plotLinesAttr)
          xnode,myh = ci
        x,h = bp[-1]
      
      a = (x - xnode) / (myh - h)
      dx = a * (myh - chh)
      xchild = xnode + dx
      
    # print node.id, xnode, si, si == cInfo.lchild, xchild, myh, chh, cInfo
    if keepPositions :
      data.x = xnode
    pylab.plot([xnode, xchild], [myh, chh], **plotLinesAttr)
    _embDrawTreeOnly(tree, si, xchild, nh, plotLinesAttr, keepPositions, coalAttr, txt)


def drawTreeOnly(tree, stree = None, plotLinesAttr = None, 
                 allTipsZero = True, xroot = None, keepPositions = False,
                 coalAttr = None, txt = False) :
  assert not stree or allTipsZero

  # A mapping from node id to (positive) node height ( with !allTipsZero, tips
  # may have > 0 height, at least one tip has 0 height). 
  nh = nodeHeights(tree, allTipsZero = allTipsZero)
  if stree is not None :
    snh = nodeHeights(stree)
    rl = _embSetLeftRight(tree, tree.root, nh, stree, snh)
    # position root in the middle
    if xroot is None :
      xroot = stree.node(stree.root).data.plotAux.center.center()

    _embDrawTreeOnly(tree, tree.root, xroot, nh, plotLinesAttr,
                    keepPositions, coalAttr, txt)
      
  else :
    # set auxiliary info per node
    rl = _setLeftRight(tree, tree.root, nh)
  
    # position root in the middle
    if xroot is None :
      xroot = (rl.right + rl.left)/2
  
      _drawTreeOnly(tree, tree.root, xroot, nh, plotLinesAttr, keepPositions, txt)
  
  for i in tree.all_ids() :
    del tree.node(i).data.cladeInfo


_Info5 = namedtuple('Info5', 'spacing extermaLeft extermaRight center height width')

def _detemineSpacing(tree, nid, epsi, po) :
  node = tree.node(nid)
  demo = node.data.demographic
  p0 = demo.population(0)
  br = node.data.branchlength

  if not node.succ :
    node.data.x = 0
    return _Info5((), [(0 - p0/2,0)], [(0 + p0/2,0)], po(node, None), br, demo.population(br))

  ch = [_detemineSpacing(tree, si, epsi, po) for si in node.succ]
  
  # get exterma points as function of spacing between the two clades

  # left - right wall of left son
  # left - left wall of right son
  exl,exr = ch[0].extermaRight, ch[1].extermaLeft
  hr = [h for x,h in exr]
  hl = [h for x,h in exl]
  myh = ch[0].height
  extermaHeights = sorted(set(hr + hl)) #  + [ch[0].height]))
  # center as a function of spacing
  S0 = sum(ch[0].spacing)
  ch[1].center.setOffset(S0)
  center = po(node, [x.center for x in ch])
  
  ihr = ihl = 0
  sval = 0
  for h in extermaHeights :
    while ihl+1 < len(hl) and not (hl[ihl] <= h < hl[ihl+1]) :
      ihl += 1
      
    if ihl+1 == len(hl) :
      # above last wall point. Line goes from last exterma point
      # to center + width
      x0,h0 = exl[-1]
      fl = lambda s,x0=x0,h0=h0 : x0 + (center.widenBy(s)-x0) * ((h-h0)/(myh-h0))
    else :
      a = (h - hl[ihl])/(hl[ihl+1]-hl[ihl])
      xl = exl[ihl][0] * (1-a) + exl[ihl+1][0] * a
      fl = lambda s,xl=xl : xl

    while ihr+1 < len(hr) and not (hr[ihr] <= h < hr[ihr+1]) :
      ihr += 1
      
    if ihr+1 == len(hr) :
      # above last wall point. Line goes from last exterma point
      # to center + width
      x0,h0 = exr[-1]
      fr = lambda s,x0=x0,h0=h0 : s + S0 + x0 + (center.widenBy(s)-(x0+S0+s)) * ((h-h0)/(myh-h0))
    else :
      a = (h - hr[ihr])/(hr[ihr+1]-hr[ihr])
      xr = exr[ihr][0] * (1-a) + exr[ihr+1][0] * a
      fr = lambda s,xr=xr : s + S0 + xr

    dx = lambda s : fr(s) - fl(s)
    dx0 = dx(0)
    u = dx(1) - dx0
    sval = max(sval,(epsi - dx0)/u)

  ch[1].center.setOffset(sval+S0)
  center = po(node, [x.center for x in ch])
  br = node.data.branchlength
  p1 = demo.population(br)
  return _Info5(ch[0].spacing + (sval,) + ch[1].spacing,
                ch[0].extermaLeft + [(center.center()-ch[0].width,myh)],
                [(x+S0+sval,y) for x,y in ch[1].extermaRight] +
                [(center.center()+ch[1].width,myh)],
                center, br+myh, p1)
    
def setSTspacing(tree, epsi, po) :
  """ Set species tree taxa positions (node.data.x) for plotting.

  The horizontal distance between extrama points (of branches) is at least epsi.
  
  """
  rs = _detemineSpacing(tree, tree.root, epsi, po)

  xs = cumsum((0,)+rs.spacing)
    
  for n,x in zip(filter(lambda n : not n.succ, getPostOrder(tree)), xs) :
    n.data.x = x
 
