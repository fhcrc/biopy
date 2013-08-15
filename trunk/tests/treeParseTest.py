import unittest
from collections import namedtuple

import biopy.parseNewick
from biopy.treeutils import getCommonAncesstor

Case = namedtuple('Case', ['newick', 'nTaxa', 'heights', 'attributes'])

class TestTreeParse(unittest.TestCase):

    def setUp(self):
        self.cases = [
          Case("(a,b)", 2, None, None),
          # blanks
          Case("(a, b)", 2, None, None), 
          Case("(a, b )", 2, None, None), 
          Case("((a,b),c)", 3, None, None), 
          Case("((a,b) ,c)", 3, None, None), 
          Case("((a,b) , c)", 3, None, None), 
          Case("( (a,b) , c)", 3, None, None), 
          Case("( (a , b) , c )", 3, None, None),
          # heights
          Case("(a:1,b:1)", 2, [ (('a','b'), 0), (('a',), 1)], None ),
          Case("(a :1,b:1)", 2, [ (('a','b'), 0), (('a',), 1)], None ),
          Case("(a : 1,b:1)", 2, [ (('a','b'), 0), (('a',), 1)], None ),
          Case("( a : 1, b:1)", 2, [ (('a','b'), 0), (('a',), 1)], None ),
          # attributes
          Case("(a[&v=1],b)", 2, None, [(('a',), (('v', '1'),)), ]), 
          Case("(a,b)[&v=1]", 2, None, [(('a','b'), (('v', '1'),)), ]), 
          Case("(a,b)[&v=1,u=2]", 2, None,
               [( ('a','b'), (('v', '1'), ('u', '2'))), ]),

          # Case("""(a[&v=1\]2],b)""", 2, None, [(('a',), (('v', """1\]2"""),)), ]), 
          
          # attributes + some blanks
          Case("(a [&v=1],b)", 2, None, [(('a',), (('v', '1'),)), ]),  
          Case("(a [&v=1] ,b)", 2, None, [(('a',), (('v', '1'),)), ]),  
          Case("(a [&v= 1] ,b)", 2, None, [(('a',), (('v', '1'),)), ]),
          Case("(a [&v =1] ,b)", 2, None, [(('a',), (('v', '1'),)), ]),
          Case("(a [& v=1] ,b)", 2, None, [(('a',), (('v', '1'),)), ]),   
          Case("(a [& v = 1] ,b)", 2, None, [(('a',), (('v', '1'),)), ]),
          
          Case("(a,b)[& v=1,u=2]", 2, None,
               [( ('a','b'), (('v', '1'), ('u', '2'))), ]),
          
          Case("(a,b)[& v=1, u= 2]", 2, None,
               [( ('a','b'), (('v', '1'), ('u', '2'))), ]), 

          Case("(a[&v=1]:1,b:1)", 2, [ (('a','b'), 0), (('a',), 1)],
               [(('a',), (('v', '1'),)), ]),
          
          Case("(a:[&v=1]1,b:1)", 2, [ (('a','b'), 0), (('a',), 1)],
               [(('a',), (('v', '1'),)), ]), 
          
          Case("(a[comment],b)", 2, None, None),
          ]

    def test_basic_parse(self):
      for case in self.cases:
        t = biopy.parseNewick.parseNewick(case.newick)
        self.assertEqual(len(t.get_terminals()), case.nTaxa)

    def test_c_parse(self):
      for case in self.cases:
        try :
          nodes1 = biopy.parseNewick.parsetree(case.newick)
        except :
          self.assertTrue(False, msg = case.newick)
          
        nodes2 = [] ; biopy.parseNewick._readSubTree(case.newick, nodes2)
        self.assertEqual(tuple(nodes1), tuple(nodes2))

    def test_detailed_parse(self):
      for case in self.cases:
        if case.heights is not None or case.attributes is not None:
          t = biopy.parseNewick.parseNewick(case.newick)
          if case.heights is not None:
            for c,h in case.heights :
              i = getCommonAncesstor(t, [t.search_taxon(x) for x in c])
              self.assertEqual(t.node(i).data.branchlength, h)
          if case.attributes is not None:
            for c,a in case.attributes :
              i = getCommonAncesstor(t, [t.search_taxon(x) for x in c])
              self.assertTrue( hasattr(t.node(i).data, "attributes"), msg=case.newick) 
              for n,v in a:
                self.assertEqual(t.node(i).data.attributes[n], v, msg=case.newick)

# add fails
# (a,b[&s=cal,])

if __name__ == '__main__':
  unittest.main()
