import treesset

def test00() :
  """
>>> ts = treesset.TreesSet()

>>> tree0txt = '(a,b)' ; i0 = ts.add(tree0txt)
>>> str(ts[i0]) == tree0txt
True

>>> tree1txt = '((a,b),c)'
>>> i1 = ts.add(tree1txt)
>>> str(ts[ i1 ]) == tree1txt
True

# Still valid
>>> str(ts[ i0 ]) == tree0txt == tree0txt
True

>>> ts = treesset.TreesSet()

# Non binary 
>>> tree2txt = '(a,b,c)'
>>> i2 = ts.add(tree2txt)
>>> str(ts[ i2 ]) == tree2txt
True

>>> ts = treesset.TreesSet()
>>> tree3txt = '((a,b),(c,d),e)'
>>> i3 = ts.add(tree3txt)
>>> str(ts[ i3 ]) == tree3txt
True

>>> ts = treesset.TreesSet()
>>> tree3atxt = '((a,b),(c,d))'
>>> i3a = ts.add(tree3atxt)
>>> str(ts[ i3a ]) == tree3atxt
True

>>> from biopy.parseNewick import parseNewick

>>> ts = treesset.TreesSet()
>>> tree4txt = '(a:1,b:1)'
>>> i4 = ts.add(tree4txt)
>>> str(ts[ i4 ]) == str(parseNewick(tree4txt))
True

>>> ts = treesset.TreesSet()
>>> tree5txt = '((a:1,b:1):2,c:3)'
>>> i5 = ts.add(tree5txt)
>>> str(ts[ i5 ]) == str(parseNewick(tree5txt))
True

>>> ts = treesset.TreesSet(compressed=False)
>>> i5 = ts.add(tree5txt)
>>> str(ts[ i5 ]) == str(parseNewick(tree5txt))
True

>>> ts = treesset.TreesSet(precision=8)
>>> i5 = ts.add(tree5txt)
>>> str(ts[ i5 ]) == str(parseNewick(tree5txt))
True

# fun with dated tips
>>> ts = treesset.TreesSet()
>>> tree6txt = '(a:1,b:2)'
>>> i6 = ts.add(tree6txt)
>>> str(ts[ i6 ]) == str(parseNewick(tree6txt))
True

>>> ts = treesset.TreesSet()
>>> tree7txt = '((a:1,b:2):3,c:3)'
>>> i7 = ts.add(tree7txt)
>>> str(ts[ i7 ]) == str(parseNewick(tree7txt))
True

>>> ts = treesset.TreesSet()
>>> tree8txt = '((a:1,b:2):3,(d:1,c:3):2)'
>>> i8 = ts.add(tree8txt)
>>> str(ts[ i8 ]) == str(parseNewick(tree8txt))
True

>>> ts = treesset.TreesSet()
>>> tree9txt = '(a:1,b:3,c:2)'
>>> i9 = ts.add(tree9txt)
>>> str(ts[ i9 ]) == str(parseNewick(tree9txt))
True

# Diffrent taxa sets
>>> ts = treesset.TreesSet()
>>> tree10txt,tree11txt = '(a,b)','(a,b,c)'
>>> i10 = ts.add(tree10txt)
>>> str(ts[ i10 ]) == str(parseNewick(tree10txt))
True
>>> i11 = ts.add(tree11txt)
>>> str(ts[ i11 ]) == str(parseNewick(tree11txt))
True
>>> str(ts[ i10 ]) == str(parseNewick(tree10txt))
True
"""
  pass

def test01() :
  """
>>> ts = treesset.TreesSet()
>>> i = ts.add('(a,b)', name = "tree0", rooted = True)
>>> t = ts[i]
>>> t.name
'tree0'
>>> t.rooted
True
"""
  pass
  
def basicFilterTest() :
  """
>>> ts = treesset.TreesSet()
>>> i = ts.add('(a,(b,c))')
>>> tsf = ts.filterTaxa('c')
>>> str(tsf[0])
'(a,b)'
>>> tsf = ts.filterTaxa('a')
>>> str(tsf[0])
'(b,c)'
>>> tsf = ts.filterTaxa('b')
>>> str(tsf[0])
'(a,c)'
>>> ts = treesset.TreesSet()
>>> i = ts.add('(a,b,c,d)')
>>> ts1 = ts.filterTaxa('c')
>>> str(ts1[0])
'(a,b,d)'
>>> ts = treesset.TreesSet()
>>> i = ts.add('((a:1,b:2):1,c:5)')
>>> ts1 = ts.filterTaxa('b')
>>> str(ts1[0])
'(a:2.0,c:5.0)'
"""
  pass

def basicLabelsTest() :
  """
>>> ts = treesset.TreesSet()
>>> i = ts.add('(a,b)E')
>>> t = ts[i]
>>> str(t)
'(a,b)E'
>>> t.node(0).data.taxon
'a'
>>> str(t)
'(a,b)E'
>>> i = ts.add('(a,(b,c))ABC')
>>> t = ts[i]
>>> str(t)
'((b,c),a)ABC'
>>> t.node(0).data.taxon
'a'
>>> str(t)
'((b,c),a)ABC'
>>> i = ts.add('((a,b),c)ABC')
>>> str(ts[i])
'((a,b),c)ABC'
  """
  pass

## ((((((10:0.036162075000000016,9:0.036162075000000016):0.06274895000000003,1:0.09891103000000001):0.026505180000000017,((13:0.014917999999999987,14:0.014917999999999987):0.03569254299999991,15:0.050610541999999814):0.07480567000000016):0.26405415,(4:0.032545126999999896,5:0.032545126999999896):0.3569252500000002):0.2710403200000002,7:0.6605106600000004):0.2432706699999998,(((16:0.024232836,6:0.024232836):0.009055312000000003,8:0.033288147):0.12789393999999998,3:0.16118209):2.3345778,((12:0.2212771,2:0.2212771):0.20966916000000002,11:0.43094626):0.47283506)

if __name__ == '__main__':
  import doctest
  doctest.testmod()
