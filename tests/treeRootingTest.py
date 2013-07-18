import unittest

from biopy.treeutils import getCommonAncesstor, treeDiameterInfo, rootAtMidpoint
from biopy.combinatorics import allPairs
from biopy.parseNewick import parseNewick

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
      
if __name__ == '__main__':
  unittest.main()
