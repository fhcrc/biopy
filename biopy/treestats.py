## This file is part of biopy.
## Copyright (C) 2013 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

"""
=============================================
Tree Statistics (for approximate phylogenies)
=============================================
"""

from __future__ import division

from math import log,exp
from collections import Counter
from itertools import count

from biopy.treeutils import *
from otus import clusterFromTree, nPairs

__all__ = ["LineagesStat", "Chao1Stat", "EntropyStat", "SimpsonHill2Stat",
           "treeStatsThroughTime",  "treeStatsThroughTimePerLoc",
           "element2Cluster", "pairCounts", "clusteringTreeError",
           "NodeInfo"]

class NodeInfo(object) :
  def __init__(self, n) :
    tx = n.data.taxon

    i = tx.find('size=')
    if i >= 0 :
      i += 5
      sz = int(tx[i:tx.index(';',i)])
    else :
      sz = 1

    i = tx.find('|')
    if i >= 0 :
      i += 1
      loc = tx[i:tx.index(';',i)]
    else :
      loc = None
      
    self.size = sz
    self.loc = loc


class LineagesStat(object) :
  def __init__(self, tree) :
    pass

  @staticmethod
  def startValue(nodesVals) :
    return sum(nodesVals)
  
  @staticmethod
  def init(node) :
    return node.data.info.size

  @staticmethod
  def combine(val, nodesVals) :
    return val - sum(nodesVals) + 1, 1

  @staticmethod
  def value(val) :
    return val

class Chao1Stat(object) :
  def __init__(self, tree) :
    pass

  @staticmethod
  def startValue(nodesVals) :
    a = sum([1 for x in nodesVals if x==1])
    b = sum([1 for x in nodesVals if x==2])
    return (len(nodesVals) + (a**2/(2*b) if b > 0 else 0), a, b)
  
  @staticmethod
  def init(node) :
    return node.data.info.size

  @staticmethod
  def combine(val, nodesVals) :
    v,a,b = val
    a -= sum([1 for x in nodesVals if x==1])
    b -= sum([1 for x in nodesVals if x==2])
    s = sum(nodesVals)
    if s == 2:
      b += 1
    return (v - 1, a, b),s

  @staticmethod
  def value(val) :
    v,a,b = val
    return v + (a**2/(2*b) if b > 0 else 0)

entropy = lambda p : p*log(p)

class EntropyStat(object) :
  def __init__(self) :
    pass

  @staticmethod
  def startValue(nodesVals) :
    N = sum(nodesVals)
    c = log(N)
    return ((sum([entropy(x) for x in nodesVals]) - c*N)/N,c,N)
  
  @staticmethod
  def init(node) :
    return node.data.info.size

  @staticmethod
  def combine(val, nodesVals) :
    s,c,N = val
    tot = sum(nodesVals)
    d = (sum([entropy(x) for x in nodesVals]) - tot*c)/N
    a = entropy(tot/N)
    return (s - d + a,c,N), tot

  @staticmethod
  def value(val) :
    return exp(-val[0])

class SimpsonHill2Stat(object) :
  def __init__(self) :
    pass

  @staticmethod
  def startValue(nodesVals) :
    N = sum(nodesVals)
    N2 = N**2
    return (sum([x**2 for x in nodesVals])/N2,N2)
  
  @staticmethod
  def init(node) :
    return node.data.info.size

  @staticmethod
  def combine(val, nodesVals) :
    s,N2 = val
    tot = sum(nodesVals)
    d = sum([x**2 for x in nodesVals])/N2
    a = tot**2/N2
    return (s - d + a,N2), tot

  @staticmethod
  def value(val) :
    return 1/val[0]

def treeStatsThroughTime(stats, tree, taxonInfo = NodeInfo, stp = None, ca = None) :
  if not ca :
    ca = CAhelper(tree)
  trms = set(tree.get_terminals())
  nint = [tree.node(i) for i in tree.all_ids() if i not in trms]
  snint = sorted(nint, key = lambda x : x.data.rh)

  tnext = stp if stp is not None else None

  for i in trms :
    n = tree.node(i)
    n.data.info = taxonInfo(n)
    n.data.val = [stat.init(n) for stat in stats]

  valsThroughTime = []
  v = [tree.node(i).data.val for i in trms]
  gval = [stat.startValue([x[k] for x in v]) for k,stat in enumerate(stats)]

  valsThroughTime.append((0,[stat.value(g) for stat,g in zip(stats,gval)]))

  for n in snint:
    chval = [tree.node(x).data.val for x in n.succ]
    gval,v = zip(*[stat.combine(g, [x[k] for x in chval]) for k,(stat,g)
                   in enumerate(zip(stats,gval))])
    n.data.val = v

    tm = n.data.rh
    if tnext is None or tm >= tnext :
      valsThroughTime.append((n.data.rh, [stat.value(g) for stat,g in zip(stats,gval)]))
      if tnext is not None :
        tnext += stp

  return valsThroughTime

def treeStatsThroughTimePerLoc(stats, tree, taxonInfo = NodeInfo, stp = None, ca = None) :
  if not ca :
    ca = CAhelper(tree)
  trms = set(tree.get_terminals())
  nint = [tree.node(i) for i in tree.all_ids() if i not in trms]
  snint = sorted(nint, key = lambda x : (x.data.rh, x.data.level))

  tnext = stp if stp is not None else None
  
  allloc = set()
  for i in trms :
    n = tree.node(i)
    d = n.data
    d.info = taxonInfo(n)
    d.val = [stat.init(n) for stat in stats]
    allloc.add(d.info.loc)
    
  allloc = sorted(allloc)
  loc2index = dict(zip(allloc,count()))

  ntrms = [tree.node(i) for i in trms]
  for n in ntrms :
    n.data.lvals = [None]*len(allloc)
    nl = loc2index[n.data.info.loc]
    n.data.lvals[nl] = n.data.val
    
  valsThroughTime = []
  v = [n.data.val for n in ntrms]
  gval = [stat.startValue([x[k] for x in v]) for k,stat in enumerate(stats)]

  glvals = []
  for nl,loc in enumerate(allloc) :
    vn = [n.data.lvals[nl] for n in ntrms]
    s = [stat.startValue([x[k] for x in vn if x is not None]) for k,stat in enumerate(stats)]
    glvals.append(s)

  xx = [[stat.value(g) for stat,g in zip(stats,glvals[nl])]
        for nl in range(len(allloc))]
  x = [stat.value(g) for stat,g in zip(stats,gval)]
  valsThroughTime.append((0,x,xx))

  for n in snint:
    chval = [tree.node(x).data.val for x in n.succ]
    gval,v = zip(*[stat.combine(g, [x[k] for x in chval])
                   for k,(stat,g) in enumerate(zip(stats,gval))])
    n.data.val = v

    n.data.lvals = [None]*len(allloc)
    #import pdb; pdb.set_trace()
    for nl,loc in enumerate(allloc) :
      chval = [tree.node(x).data.lvals[nl] for x in n.succ]
      gl = glvals[nl]

      #chv = [x[k] if x is not None else None  for x in chval]
      if all([x is not None for x in chval]) :
        gl,v = zip(*[stat.combine(g, [x[k] for x in chval])
                     for k,(stat,g) in enumerate(zip(stats,gl))])
        n.data.lvals[nl] = v
      elif not all([x is None for x in chval]) :
        n.data.lvals[nl] = chval[0] if chval[0] is not None else chval[1]
        
      glvals[nl] = gl
    
    x = [stat.value(g) for stat,g in zip(stats,gval)]
    xx = [[stat.value(g) for stat,g in zip(stats,glvals[nl])]
          for nl in range(len(allloc))]

    tm = n.data.rh
    if tnext is None or tm >= tnext :
      valsThroughTime.append((tm, x, xx))
      if tnext is not None :
        tnext += stp

  return valsThroughTime, allloc




def element2Cluster(clust) :
  d = dict()
  for k,c in enumerate(clust) :
    for x in c:
      d[x] = k
  return d

def pairCounts(cref, ctry) :
  rcref = sorted(cref, key = lambda x : len(x), reverse = 1)
  rctry = sorted(ctry, key = lambda x : len(x), reverse = 1)

  revtry = element2Cluster(rctry)
  revref = element2Cluster(rcref)

  m = [(sum([nPairs(k) for k in
             Counter([revtry[x] for x in c]).values()]), nPairs(c))
       for c in rcref]

  w = [(nPairs(c) - sum([nPairs(k) for k in
                       Counter([revref[x] for x in c]).itervalues()]),
        nPairs(c))
       for c in rctry]

  inin,inout = sum([x[0] for x in m]), sum([x[1]-x[0] for x in m])
  outin = sum([x[0] for x in w])

  r = ((inin,inout),(outin, nPairs(sum([len(x) for x in cref]))-(inin+inout+outin)))
  return r

def subTreeTL(cTree, nRoot, txs) :
  tx = [x.data.taxon for x in nRoot.data.terms]
  assert len(txs) > 0 and txs.issubset(tx)

  # single path from tip to node
  if len(txs) == 1 :
    return nRoot.data.rh

  # full house
  if nRoot.data.cladesize == len(txs) :
    return nRoot.data.tl
  
  assert nRoot.succ

  tl = 0
  for ch in [cTree.node(x) for x in nRoot.succ] :
    # taxa in child
    r = txs.intersection([x.data.taxon for x in ch.data.terms])
    if len(r) :
      tl += subTreeTL(cTree, ch, r) + ch.data.branchlength
  return tl

# clu is a partition of taxa in cTree.
# th is the height (distance/2)

def clusteringTreeError(cTree, clu, th, full = False, verbose = False) :
  ca = CAhelper(cTree)
  
  txToNode = dict()

  for t in cTree.get_terminals() :
    n = cTree.node(t)
    txToNode[n.data.taxon] = n
  assert len(txToNode) == len(cTree.get_terminals())

  if full :
    report = []
  tot = 0
  for nc,txs in enumerate(clu) :
    #txs = [tx[x] for x in c]
    n = ca.getCAi(txs)
    d = n.data

    if len(txs) == 1 :
      tot += th
      if full:
        report.append((nc, n.id, len(txs), 0, th, th))
    else :
      tl = subTreeTL(cTree, n, set(txs)) if len(txs)>1 else 0
      if full:
        #n.id, cTree.toNewick(n.id,topologyOnly=1), d.tl, th - d.rh, max(th - d.rh, 0))
        report.append((nc, n.id, len(txs), tl, th - d.rh, max(th - d.rh, 0)))
      if verbose and len(txs) > 1:
        # cTree.toNewick(n.id,topologyOnly=1)
        print nc, n.id, d.cladesize, tl, th - d.rh, max(th - d.rh, 0)
      tot += tl + max(th - d.rh, 0)

  otot = 0
  trueClust = clusterFromTree(cTree, th, ca)
  for n in trueClust :
    otot += n.data.tl + (th - n.data.rh)

  return (tot, otot, th*len(cTree.get_terminals())) + ((report,) if full else tuple())

