## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

import random, operator
from math import log
from numpy import array, cumsum, inf, mean

from scipy.optimize import brentq

from treeutils import nodeHeight, treeHeight, nodeHeights
import demographic


__all__ = ["GeneTreeSimulator"]


def _getXsList(ds, x0, x1) :
  """ Get X-axis points in [x0,x1] from all demographics in ds."""
  xs = []
  for d in ds:
    if isinstance(d, demographic.ConstantPopulation) :
      continue
    xs.extend([z for z in d.xvals if x0 < z < x1])

  xs = sorted([x0] + xs + [x1])
  # remove duplications
  xs = [x for k,x in enumerate(xs) if k == len(xs)-1 or x != xs[k+1]]
  return xs

def _isZero(demo, x0, dx) :
  
  if demo.population(x0) != 0 or demo.population(x0+dx) != 0 :
    return False

  for x in _getXsList([demo], x0, x0+dx) :
    if demo.population(x) != 0 :
      return False

  return True
  
def _integrateLinRatio(t, a,b,c,d) :
  """ Integrate (a+b*x)/(c+d*x) in [0,t] """
  return ((a*d - b*c)*log((d/c)*t + 1) + b*d*t)/d**2

def _integrateLin2Ratio(t, a1, b1, a2, b2, a3, b3) :
  """ Integrate (a1+b1*x)*(a2+b2*x)/(a3+b3*x) in [0,t] """

  # Quadratic term factored out
  x1,x2 = - a3**2*b1*b2/b3**2 + a1*a2, (-2*a3*b1*b2/b3 + a1*b2 + a2*b1)
  y1,y2 =   a3, b3

  p1 = ((b1*b2)/b3) * (t**2/2 + (a3/b3) * t)
  p2 = _integrateLinRatio(t, x1,x2,y1,y2)

  return p1+p2

def _dRatioIntegrate(d1, d2, x) :
  """ Integrate (d1(t) / d2(t)) in [0, x].

  d1,d2 are piecewise linear.
  """

  assert isinstance(d1, demographic.LinearPiecewisePopulation) \
         and isinstance(d2, demographic.LinearPiecewisePopulation)

  if x == 0 :
    return 0
  
  xs = [0] + sorted([z for z in d1.xvals + d2.xvals if z < x]) + [x]

  cum = 0
  for i in range(len(xs)-1) :
    x1, x2 = xs[i], xs[i+1]
    y11, y12 = [d1.population(z) for z in [x1,x2]]
    y21, y22 = [d2.population(z) for z in [x1,x2]]

    dx = x2-x1
    if dx > 0 :
      a,b = y11, (y12-y11)/dx
      c,d = y21, (y22-y21)/dx

      if d == 0 :
        #print a,b,c,d, dx, dx * (0.5 * dx * b + a) / c
        cum += dx * (0.5 * dx * b + a) / c
      else :
        cum += _integrateLinRatio(x2-x1, a,b,c,d)

  return cum

def _dRatio2Integrate(x0, x1, d1, d2, d3) :
  """ Integrate (d1(t)*d2(t) / d3(t)) in [x0, x1]. """
  assert x0 <= x1
  if x0 == x1 :
    return 0
  
  xs = _getXsList([d1,d2,d3], x0, x1)

  isStep = isinstance(d1, demographic.StepFunctionPopulation)
  
  cum = 0
  for i in range(len(xs)-1) :
    x1, x2 = xs[i], xs[i+1]

    if isStep :
      ii = d1.population((x1+x2)/2) * \
           (dRatioIntegrate(d2, d3, x2) - dRatioIntegrate(d2, d3, x1))
    else :
      y11, y12 = [d1.population(z) for z in [x1,x2]]
      y21, y22 = [d2.population(z) for z in [x1,x2]]
      y31, y32 = [d3.population(z) for z in [x1,x2]]

      dx = x2-x1
      a,b = y11, (y12-y11)/dx
      c,d = y21, (y22-y21)/dx
      e,f = y31, (y32-y31)/dx

      if f == 0 :
        ii = (((b*d * dx**3/3) + (a*d + b*c) * dx**2/2 + a*c*dx)/e)
      else :
        ii = _integrateLin2Ratio(dx, a,b,c,d,e,f)
    cum += ii
  return cum

def _cumMigrationB2AD(t1, t2, p_b2a, aPop, bPop) :
  """ Cumulative migration (integral(p_b2a(x)*b(x)/a(x) dx) in [t1,t2]"""
  assert t1 <= t2
  return _dRatio2Integrate(t1, t2, p_b2a, bPop, aPop)

def _timeToNextA2BmigrationD(t, nl_a, p_b2aD, pa_d, pb_d) :
  """ Random waiting time starting at t for an A lineage to 'jump' to B.
  In other words for the event of A lineage has an ancestor which is a migrant
  from B.
  
  Rate of B migration and a,b populations are parameters.
  """
  if nl_a == 0 :
    return inf
  
  r = -log(1 - random.random())
  f = lambda x : nl_a * _cumMigrationB2AD(t, x, p_b2aD, pa_d, pb_d) - r
  assert f(t) < 0
  h = t + 100
  if f(h) == -r :
    return inf
  
  while f(h) < 0 :
    h *= 2
  return brentq(f, t, h) - t


#def sampleTaxaName(species, k) :
#  return species + ("_%03d" % k)


def _simplify(dls, x1) :
  """ Create a demographic on [0,x1] eqvivalent to the sum of all dls's."""
  xs = _getXsList(dls, 0, x1)
  pops = [sum([d.population(x) for d in dls]) for x in xs]
  return demographic.LinearPiecewisePopulation(pops, xs[1:])

def _simplifyc(dls, x1) :
  """ Create a demographic on [0,x1] eqvivalent to the sum of all dls's."""
  if len(dls) == 1 :
    return dls[0]
  xs = _getXsList(dls, 0, x1)
  # for const/step, use middle 
  pops = [sum([d.population((x+y)/2) for d in dls]) for x,y in zip(xs[:-1], xs[1:])]
  if len(pops) == 1 or len(pops) == 2 and pops[0] == pops[1] :
    return demographic.ConstantPopulation(pops[0])
  return demographic.StepFunctionPopulation(pops, xs[1:-1])

def _aa(x, tls, dls) :
  assert x <= tls[-1]
  k = 0
  if x > tls[0] :
    k = 1
    while k < len(tls) :
      if tls[k-1] < x <= tls[k] :
        break
  return sum([d.population(x) for d in dls[k]])

def _dcombine(tls, dls, isStep) :
  """ 0 - tls[0]   sum(dls[0])
     tls[1] - tls[0]  sum(dls[1])
     ...
  """
  if isStep :
    #import pdb; pdb.set_trace()
    sdls = [_simplifyc(d, x1-x0) for x0,x1,d in zip([0] + tls[:-1], tls, dls)]
    xs = []
    ys = []
    for x0,x1,d in zip([0] + tls[:-1], tls, sdls) :
      if isinstance(d, demographic.ConstantPopulation) :
        ys.append(d.population(0))
      else :
        xs.extend([x0 + x for x in d.xvals])
        ys.extend(d.vals)
      xs.append(x1)
    xs = xs[:-1]
    if len(xs) :
      demo = demographic.StepFunctionPopulation(ys, xs)
    else :
      demo = demographic.ConstantPopulation(ys[0])
  else :
    sdls = [_simplify(d, x1-x0) for x0,x1,d in zip([0] + tls[:-1], tls, dls)]
    xs = []
    ys = [sdls[0].population(0)]
    for x0,x1,d in zip([0] + tls[:-1], tls, sdls) :
      xs.extend([x0 + x for x in d.xvals])
      ys.extend(d.vals[1:])

    demo = demographic.LinearPiecewisePopulation(ys, xs)

  #import pdb; pdb.set_trace()  
  return demo

class LineageInfo(object) :
  """ Gene subtrees of one sp (extant or ancestral) """
  def __init__(self, node, nh, trees) :
    self.node = node
    assert nh is not None
    self.nodeH = nh # nodeHeight(tree, node.id) # tree.nodeHeight(node)
    self.trees = trees

  def __str__(self) :
    ts = ["%g [%s]" % (h,t) for t,h in self.trees]
    return "%s %g # %s" % (self.node.id, self.nodeH, " # ".join(ts))
  
  def popAtAbsTime(self, t) :
    d = self.node.data.demographic
    # assert  d.naturalLimit() is not None, str(d)
    assert t >= self.nodeH and \
           (d.naturalLimit() is None or t <= self.nodeH + d.naturalLimit()), (t, str(d), self.nodeH,)
    
    return d.population(t - self.nodeH)

# list s-node: current trees are between node and nodes parent,
#      trees: list of lineages as trees text, each with root height
#  list of migration "pairs" (sa1,sa2,sa3,...), (sb1, sb2,...) with 2 rate
#      functions.

# Generate a coal waiting time for all species. Take the minimum.
# generate waiting times for all migrations.

# Should change it to use tree builder

class GeneTreeSimulator(object) :
  def __init__(self, tree) :
    self.tree = tree
    self.nhts = nodeHeights(tree)
    assert all([isinstance(tree.node(n).data.demographic,
                          (demographic.StepFunctionPopulation,
                           demographic.ConstantPopulation,
                           demographic.LinearPiecewisePopulation))
                for n in tree.all_ids()])
    
    nlp = sum([isinstance(tree.node(n).data.demographic,demographic.LinearPiecewisePopulation)
               for n in tree.all_ids()])
    assert nlp == 0 or nlp == len(tree.all_ids())
    self.isStepDemographic = nlp == 0

  def _startup(self, nodeId) :
    #def _startup(self, nodeId, tipNameFunc) :
    """ Return (linages, mp, species, demographics)

       linages - a sequence with linages information (LineageInfo) per species tree
       node (i.e. one sp).

       mp - sequence of ( (taxa-list-left,taxa-list-right) , (M-function-left,
         M-function-right), (total-popsize-left, total-popsize-right) )

       species - names (taxa) of all extant species for this ancestal species.

       demographics - total popsize, demographic from time 0 to node
         height. Value at time t is the sum of pop sizes
         for all subspecies in the clade.
    """

    tree = self.tree
    
    node = tree.node(nodeId)
    if node.data.taxon :
      tips = list(node.data.labels)
      #lins = [LineageInfo(tree, node, [list(x) for x in zip(tips,
      #[0.0,]*len(tips))])]
      lins = [LineageInfo(node, self.nhts[nodeId],
                          [list(x) for x in zip(tips, [0.0,]*len(tips))])]
      sps = [node.data.taxon]
      dms = node.data.demographic
      mp = []
    else :
      s1,s2 = [self._startup(child) for child in node.succ]
      lins = s1[0] + s2[0]
      sps = s1[2] + s2[2]
      nh = self.nhts[node.id]
      #from treeutils import toNewick
      #assert nh == nodeHeight(tree, node.id), (nh, node.id, nodeHeight(tree,
      # node.id), toNewick(tree),  nodeHeights(tree),
      # [(x,nodeHeight(tree,x),tree.node(x).data.branchlength) for x in tree.all_ids()])
      br = node.data.branchlength
      #print ([nh, nh+br], [[s1[1][-1][2], s2[1][-1][2]], [node.data.demographic]])
      #dms = dcombine([nh, nh+br], [[s1[1][-1][2], s2[1][-1][2]], [node.data.demographic]])
      dms = _dcombine([nh, nh+br], [ [s1[3], s2[3]], [node.data.demographic] ], self.isStepDemographic)
      #print dms
      mp = s1[1] + s2[1] + [((s1[2],s2[2]), node.data.ima, (s1[3], s2[3]))]

    if nodeId == tree.root:
      d1 = dict([(l.node.data.taxon, k) for k,l in enumerate(lins)])
      for k,m in enumerate(mp):
        mp[k] = (tuple([d1[x] for x in idl] for idl in m[0]),) + mp[k][1:]

    return (lins, mp, sps, dms)

  def simulateGeneTree(self) :
    tree = self.tree
    sps, mp = self._startup(tree.root)[:2]

    nIM = 0
    
    t = 0

    #lim = min([(tree.nodeHeight(tree.node(sp.node.prev)),nk) for nk,sp in
    #enumerate(sps)])
    #lim = min([(nodeHeight(tree, sp.node.prev),nk) for nk,sp in enumerate(sps)])
    lim = min([(self.nhts[sp.node.prev],nk) for nk,sp in enumerate(sps)])

    while not (len(sps) == 1 and len(sps[0].trees) == 1) :
      # Generate a coal waiting time for all species. 
      
      atimes = []
      for nk,sp in enumerate(sps) :
        nl = len(sp.trees)
        demog = sp.node.data.demographic
        wt = demog.timeToNextCoalescent(nl, t - sp.nodeH) if nl > 1 else inf
        atimes.append((wt, nk)) 

      # Take the minimum. We can optimize based on the fact that any migration
      # waiting time greater than that would be rejected
      cwt,cnk = min(atimes)
      
      itimes = []
      #import pdb; pdb.set_trace()
      for nk, ((linLeft,linRight), (pLeft, pRight), (dLeft, dRight)) in enumerate(mp) :

        if not _isZero(pRight, t, cwt) :
          # a left, b right
          nl_a = sum([len(sps[i].trees) for i in linLeft])

          wt1 = _timeToNextA2BmigrationD(t, nl_a, pRight, dLeft, dRight)
          # 0 means a (left) jumps into b (right) pool
          itimes.append((wt1, nk, 0))                     ; assert pRight.population(t + wt1) > 0
          
        if not _isZero(pLeft, t, cwt) : 
          nl_b = sum([len(sps[i].trees) for i in linRight])
          
          wt2 = _timeToNextA2BmigrationD(t, nl_b, pLeft, dRight, dLeft)
          # 1 means b (right) jumps into a (left) pool
          itimes.append((wt2, nk, 1))                     ; assert pLeft.population(t + wt2) > 0

      iwt,ink,idir = min(itimes) if len(itimes) else (inf,-1,-1)

      wt = min(iwt, cwt)
      
      if 0:
        print "t",t
        print "lim",lim
        print "atimes", atimes
        print "itimes",itimes
        print "wt", wt
        print
        
      if t + wt > lim[0]:
        # Advance to next speciation, not necessarily the one whose waiting time was exceeded 
        n1 = lim[1]
        # join the two lineages
        p = sps[n1].node.prev
        for nk,sp in enumerate(sps):
          if nk != n1 and sp.node.prev == p :
            break
        lins = sps[n1].trees + sps[nk].trees
        np = tree.node(p)
        lmerged, ljoins = min(n1,nk), max(n1,nk)
        #sps[lmerged] = LineageInfo(tree, np, lins)
        sps[lmerged] = LineageInfo(np, self.nhts[p], lins)        
        t = sps[n1].nodeH
        sps.pop(ljoins)

        lim = min([(self.nhts[sp.node.prev],nk)   # (tree.nodeHeight(tree.node(sp.node.prev)),nk)
                   if sp.node.prev is not None else (inf,nk)
                   for nk,sp in enumerate(sps)])

        #lim = min([(nodeHeight(tree, sp.node.prev),nk)   # (tree.nodeHeight(tree.node(sp.node.prev)),nk)
        #           if sp.node.prev is not None else (inf,nk)
        #           for nk,sp in enumerate(sps)])

        found = False
        for k,m in enumerate(mp) :
          if sorted(m[0][0] + m[0][1]) == [lmerged, ljoins] :
            mp.pop(k)
            found = True
            break

        assert found
        assert all([sorted(m[0][0] + m[0][1]) != [lmerged, ljoins] for m in mp])

        for k,m in enumerate(mp) :
          assert all([(ljoins in idl) == (lmerged in idl) for idl in m[0]])
          
          ni = lambda x : x if x < ljoins else x - 1
          mp[k] = (tuple([ni(x) for x in idl if x != ljoins] for idl in m[0]),) + mp[k][1:]
        
      elif cwt < iwt :
        wt, nk = cwt,cnk 
        sp = sps[nk]
        
        i,j = random.sample(range(len(sp.trees)), 2)
        l1,l2 = sp.trees[i],sp.trees[j]

        # height in tree axis
        th = t + wt

        # Build in newick format
        label = '(' + l1[0] + ':' + str(th - l1[1]) + ',' \
                + l2[0] + ':' + str(th - l2[1]) + ')'

        sp.trees.pop(max(i,j))
        sp.trees[min(i,j)] = [label, th]

        t = th
      else :
        #import pdb; pdb.set_trace()
        # iwt < cwt
        linLeft,linRight = mp[ink][0]
        # move the lineage
        # idir: 0 means a (left) jumps into b (right) pool
        
        linFrom, linTo = (linLeft,linRight) if idir == 0 else (linRight,linLeft)
        
        lcan = reduce(operator.add, [zip([i]*len(sps[i].trees), range(len(sps[i].trees)))
                                     for i in linFrom])
        i,k = random.choice(lcan)
        stree = sps[i].trees.pop(k)

        # advance time to event
        t += iwt
        pops = array([sps[i].popAtAbsTime(t) for i in linTo])
        aprobs = cumsum(pops / pops.sum())

        ito = (random.random() < aprobs).argmax()

        sps[linTo[ito]].trees.append(stree)

        nIM += 1

    return (sps[0].trees[0], nIM)


def setIMrates(stree, mSpec, sSpec, balanced = None, restrictToTree = True) :
  nh = nodeHeights(stree)
  
  # internals, leaves to root order
  internals = set(stree.all_ids()) - set(stree.get_terminals())
  internals = [z for y,z in sorted([(nh[x], x) for x in internals])]

  for x in internals :
    n = stree.node(x)
    pim = mSpec.sample()
    
    if sSpec is not None:
      h = nh[x]
      tOrig = t = sSpec.sample()
      if restrictToTree:
        # Clip migration "stop" times if they go beyond stop time of children 
        for ch in n.succ:
          nch = stree.node(ch)
          if not nch.data.taxon:
            assert nh[n.id] >= nh[nch.id]

            bound = nh[ch] - nch.data.imc
            if h - t < bound:
              t = h - bound

      if h > t :
        d = demographic.LinearPiecewisePopulation([0, 0, pim], [h-t, h])
        n.data.imc = t
      else :
        if h == t :
          assert tOrig > h
          t = tOrig

        d = demographic.LinearPiecewisePopulation([(1 - h/t)*pim , pim], [h])
        n.data.imc = h

      n.data.ima = (d, d)
    else: # sSpec is None
      d = demographic.ConstantPopulation(pim)
      if balanced :
        n.data.ima = (d, d)
      else :
        d1 = demographic.ConstantPopulation(mSpec.sample())
        n.data.ima = (d, d1)
      
if 0 :
  import randomDistributions

  def setIMrates(stree, pim = 0.05, spr = 0.15, bdGrowth = None,
                 restrictToTree = True) :
    # growth rate = lambda - mu
    # Distribution of times from split until complete separation
    if bdGrowth is not None :
      mm = spr*(0.5/bdGrowth)
    else :
      # Crude: take mean from tree 
      mm = mean([stree.node(n).data.branchlength
                 for n in stree.all_ids() if n != stree.root])

    e1 = randomDistributions.LogNormal(mm, .25)

    nh = nodeHeights(stree)

    # internals, leaves to root order
    internals = set(stree.all_ids()) - set(stree.get_terminals())
    internals = [z for y,z in sorted([(nh[x], x) for x in internals])]

    for x in internals :
      n = stree.node(x)
      h = nh[x]
      tOrig = t = e1.sample()
      if restrictToTree:
        # Clip migration "stop" times if they go beyond stop time of children 
        for ch in n.succ:
          nch = stree.node(ch)
          if not nch.data.taxon:
            assert nh[n.id] >= nh[nch.id]

            bound = nh[ch] - nch.data.imc
            if h - t < bound:
              t = h - bound

      if h > t :
        d = demographic.LinearPiecewisePopulation([0, 0, pim], [h-t, h])
        n.data.imc = t
      else :
        if h == t :
          assert tOrig > h
          t = tOrig

        d = demographic.LinearPiecewisePopulation([(1 - h/t)*pim , pim], [h])
        n.data.imc = h

      n.data.ima = (d, d)
