import unittest

from biopy.treeutils import getCommonAncesstor, treeDiameterInfo, \
     rootAtMidpoint, getPostOrder
from biopy.genericutils import flatten

# , _populateTreeWithNodeToTipDistances, _cleanTreeWithNodeToTipDistances

from biopy.combinatorics import allPairs
from biopy.parseNewick import parseNewick
from itertools import chain

def pathLenTo(tree, fr, to) :
  p = 0
  while fr != to :
    p += fr.data.branchlength
    fr = tree.node(fr.prev)
  return p

def allTipDistances(tree) :
  tr = [tree.node(x) for x in tree.get_terminals()]
  res = dict()
  for x,y in allPairs(tr) :
    a = getCommonAncesstor(tree, [x.id,y.id])
    a = tree.node(a)
    p = pathLenTo(tree, x, a) + pathLenTo(tree, y, a)
    res[frozenset((x.id,y.id))] = p
  return res

def maxTipPathTheHardWay(tree) :
  r = allTipDistances(tree)
  k = max(r, key = r.get)
  return (r[k],) + tuple([tree.node(y) for y in k])
  ## #ca = CAhelper(tree)
  ## tr = [tree.node(x) for x in tree.get_terminals()]
  ## mx = (-1, None)
  ## for x,y in allPairs(tr) :
  ##   a = getCommonAncesstor(tree, [x.id,y.id])
  ##   a = tree.node(a)
  ##   p = pathLenTo(tree, x,a) + pathLenTo(tree, y,a)
  ##   if p > mx[0] :
  ##     mx = (p, x, y)
  ## return mx

def _populateTipDistancesFromParent(tree, n, parDists) :
  if n.id != tree.root :
    assert not n.data.dtips[-1] and parDists
    n.data.dtips[-1] = [ [a[0],a[1] + n.data.branchlength] for a in parDists]
    parDists = n.data.dtips[-1]
  else :
    assert n.data.dtips[-1] and not parDists
    parDists = []

  for i in range(len(n.succ)) :
    d = flatten([n.data.dtips[j] for j in range(len(n.succ)) if j != i] + [parDists])
    _populateTipDistancesFromParent(tree, tree.node(n.succ[i]), d)

def _populateTreeWithNodeToTipDistances(tree) :
  for n in getPostOrder(tree) :
    if not n.succ:
      n.data.dtips = [[[n,0]],[],[]]
    else :
      ch = [tree.node(c) for c in n.succ]
      n.data.dtips = [[[a[0],a[1]+x.data.branchlength] for a in x.data.dtips[0]] +
                      [[a[0],a[1]+x.data.branchlength] for a in x.data.dtips[1]]
                      for x in ch]
      if n.id != tree.root :
         n.data.dtips.append([])

  _populateTipDistancesFromParent(tree, tree.node(tree.root), [])

    
def tipToInternalDistance(tree, nTip, nNode, dx) :
  """ As straighforward code as I can make it, still not exactly naive
  
  nTip - tip node
  nNode - internal node
  dx - distance from internal node towards parent (dx < branch
  """
  
  nodePath = [nNode.id]
  pathDistance = [-dx]
  n = nNode
  while n.id != tree.root :
    nodePath.append(n.prev)
    pathDistance.append( pathDistance[-1] + n.data.branchlength )
    n = tree.node(n.prev)
  
  pathDistance[0] = dx
  #print nodePath, pathDistance
  
  nid = nNode.id
  distance = 0
  while nTip.id not in nodePath :
    distance += nTip.data.branchlength
    nTip = tree.node(nTip.prev)

  #print distance, nTip.id
  return distance + pathDistance[nodePath.index(nTip.id)]


class TestTreeRooting(unittest.TestCase):
  def setUp(self) :
    self.trees = (
      '((a:1.0,b:2.0):1.0,c:5.0,d:3.0)',
      '(((a:1.0,b:2.0):4.0,c:5.0):1.0,d:3.0)',
      '((((((a1061:0.00337,a2565:0.00201):0.0643,a5868:0.06113):0.09594,((a2852:0.0,a6725:0.0):0.12791,a6754:0.097):0.03236):0.0116,(((a1271:0.00457,a1407:0.00576):0.02102,a7324:0.0):0.1118,(a3164:0.0,a5982:0.01176):0.14367):0.02905):0.00626,(((a3918:0.00371,a719:0.00792):0.08,(a4356:0.0,a5341:0.00756):0.09624):0.04958,a8397:0.0938):0.0406):0.0038,((((a2686:0.00404,a2939:0.00014):0.0954,a7018:0.05497):0.04348,(a4732:0.00772,a4856:0.0):0.14393):0.02,((a2073:0.00143,a5568:0.00343):0.12645,(a2418:0.00168,a7503:0.01774):0.12139):0.03171):0.01743)',
      )
    self.trees = [parseNewick(t) for t in self.trees]

  def test_diameter(self) :
    for t in self.trees:
      d,a,b,c = treeDiameterInfo(t)
      d1,a1,b1 = maxTipPathTheHardWay(t)
      self.assertEqual( set((a1.data.taxon,b1.data.taxon)),
                        set((a.data.taxon,b.data.taxon)) )
      self.assertEqual(d , d1)

  def test_mid(self) :
    for t in self.trees:
      tm = rootAtMidpoint(t)
      tm = parseNewick(tm)
      r1 = allTipDistances(t)
      u1 = dict([(frozenset([t.node(y).data.taxon for y in x]),r1[x]) for x in r1])
      r2 = allTipDistances(tm)
      u2 = dict([(frozenset([tm.node(y).data.taxon for y in x]),r2[x]) for x in r2])
      self.assertEqual(u1, u2)

  def test_tipDistances(self) :
    for t in self.trees:
      _populateTreeWithNodeToTipDistances(t)
      for nid in t.all_ids() :
        n = t.node(nid)
        e2 = dict([(tt,tipToInternalDistance(t, t.node(t.search_taxon(tt)), n, 0))
                   for tt in t.get_taxa()])
        ii = chain(*[d for d in n.data.dtips if d])
        v = sum([(x-y)**2 for x,y in zip(sorted(e2.values()), sorted([x[1] for x in ii]))])
        self.assertEqual( float("%f" % v),0.0)
      
if __name__ == '__main__':
  unittest.main()





if 0:
  import array
  from itertools import imap, chain
  
  def _populateTipDistancesFromParent(tree, n, parDists) :
    if n.id != tree.root :
      assert not n.data.dtips[-1] and parDists

      i = imap(lambda x : x + n.data.branchlength, parDists)
      parDists = n.data.dtips[-1] = array.array('f', i)
    else :
      assert n.data.dtips[-1] and not parDists
      parDists = []

    for i in range(len(n.succ)) :
      d = chain(*([n.data.dtips[j] for j in range(len(n.succ)) if j != i] + [parDists]))
      _populateTipDistancesFromParent(tree, tree.node(n.succ[i]), d)

  def _populateTreeWithNodeToTipDistances(tree) :
    for n in getPostOrder(tree) :
      if not n.succ:
        n.data.dtips = [array.array('f',[0]),[],[]]
      else :
        ch = [tree.node(c) for c in n.succ]
        n.data.dtips = [[a+x.data.branchlength for a in x.data.dtips[0]] +
                        [a+x.data.branchlength for a in x.data.dtips[1]]
                        for x in ch]
        if n.id != tree.root :
           n.data.dtips.append([])

    _populateTipDistancesFromParent(tree, tree.node(tree.root), [])

  def _cleanTreeWithNodeToTipDistances(tree) :
    for nid in tree.all_ids() :
      n = tree.node(nid)
      del n.data.dtips

if 0 :
  def _rootPointByTipVarianceOptimization(tree) :
    _populateTreeWithNodeToTipDistances(tree)
    minLoc = float('inf'),None,None

    for nid in tree.all_ids() :
      if nid == tree.root :
        continue
      n = tree.node(nid)
      ## pl,mn = [flatten([[a[1] for a in d] for d in n.data.dtips[:-1] if d])] + \
      ##         [[a[1] for a in n.data.dtips[-1]]]
      pl,mn = array.array('f',chain(*[d for d in n.data.dtips[:-1] if d])), n.data.dtips[-1]

      nl = len(mn)+len(pl)
      spl, smn = sum(pl), sum(mn)

      b,c = 2 * (spl - smn)/nl, sum([x**2 for x in pl + mn])/nl
      a1,b1 = (len(pl) - len(mn))/nl, (spl + smn)/nl

      ac,bc,cc = (1 - a1**2),  (b - (2 * a1 * b1)), (c - b1**2)

      dx = min(max(-bc / (2 * ac) , 0), n.data.branchlength)

      val = dx**2 * ac + dx * bc +  cc
      #print n.id,dx,val
      if val < minLoc[0] :
        minLoc = (val, n, dx)

    _cleanTreeWithNodeToTipDistances(tree)
    return minLoc
