## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

from __future__ import division

from numpy import mean, median, corrcoef
import pylab
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
    return hpar
  
  return (centerc, h, hpar, xright, p, p1, clade)


def treeMinSpace(tree, nid = -1) :
  # 0:node, 1:width, 2:widthatend, 3:maxspace
  if nid == -1 :
    nid = tree.root
  node = tree.node(nid)
  d = node.data.demographic
  if not node.succ:
    return (node, d.population(0),
            d.population(node.data.branchlength), 0)

  l,r = [treeMinSpace(tree, x) for x in node.succ]
  spc = 2*max(l[2] - l[1], r[2] - r[1], 0)

  spc = max(spc, l[3],r[3])
  if nid == tree.root:
    return spc
  return (node, spc +  l[1] + r[1], d.population(d.naturalLimit()),spc)

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
        distance[tuple(sorted([xl,xr]))] = n-1
    return (n, l[1] + r[1])

def getTaxa(tree, nid, rand = False) :
  node = tree.node(nid)
  if not node.succ:
    return [node.data.taxon]
  else:
    l,r = [getTaxa(tree, x, rand) for x in node.succ]
    if rand :
      if random.random() < .5:
        l,r = r,l
    return l + r

def getCR(oo, dis) :
  x,y = [],[]
  for k in dis:
    tx1,tx2 = k
    x.append(abs(oo.index(tx1) - oo.index(tx2)))
    y.append(dis[k])

  cr = corrcoef(x,y)[0,1]
  return cr

def getTaxaOrder(trees, refTree = None) :
  tops = [toNewick(tree, topologyOnly=1) for tree in trees]
  dtops = dict([(x,0) for x in tops])
  for t in tops :
    dtops[t] += 1

  dis = dict()
  for tx in dtops:
    tree = INexus.Tree(tx)
    d = dict()
    taxaDistance(tree, tree.root, d)
    for k in d:
      if k not in dis:
        dis[k] = 0.0
      dis[k] += d[k] * dtops[tx]/len(trees)
    
  if refTree is None:
    refTree = INexus.Tree(sorted(dtops, key = lambda x : dtops[x])[-1])

  tree = refTree
  moo = getTaxa(tree, tree.root)
  mcr = getCR(moo, dis)

  while True:
    mnode = None

    for nid in tree.all_ids() :
      node = tree.node(nid)
      if node.succ:
        node.succ = [node.succ[1], node.succ[0]]
        oo = getTaxa(tree, tree.root)
        cr = getCR(oo, dis)

        if cr > mcr:
          mcr = cr
          moo = oo
          mnode = node
        node.succ = [node.succ[1], node.succ[0]]
    if mnode is not None:
      mnode.succ = [mnode.succ[1], mnode.succ[0]]
    else :
      break

  return (moo, tree)

