## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

"""
Generic combinatorial utilities
===============================
"""


from __future__ import division
from math import log
import operator

from scipy.special import gammaln

__all__ = ["listcombinations", "permutations", "factorial", "nPairs",
           "allPairs",  "uniformvec", "lbinomial", "uniformLike", "choose",
           "lchoose"] 

def _prod(x) :
  """ Element-wise product of vector x.
  """
  
  return reduce(operator.mul, x, 1L)

def factorial(n) :
  """ n!

  :rtype: long
  """
  return _prod([long(x) for x in range(1, n+1)])

def choose(n,k) :
  """ Binomial: Number of ways to choose k of n.

  :rtype: long
  """
  return _prod(range(k+1,n+1)) // _prod(range(2,n-k+1))

def lchoose(n,k) :
  """ Log of choose(n,k). """

#  Uses scipy.special.gammaln, log of the Gamma function.

  return special.gammaln(n+1) - special.gammaln(k+1) - special.gammaln(n-k+1)

def nPairs(k) :
  """ choose(k,2) (Convenience).

  :rtype: int
  """
  return (k * (k-1)) // 2

# Hide recursion internal parameters

# http://jace.seacrow.com/archive/2007/02/15/generating-combinations-in-python

def listcombinations(listoflists) :
  """
    Generator that yields all possible combinations from a list of lists.

    >>> a = [[1, 2], [3, 4], [5, 6], [7, 8]]
    >>> for c in listcombinations(a): print c
    ...
    [1, 3, 5, 7]
    [1, 3, 5, 8]
    [1, 3, 6, 7]
    [1, 3, 6, 8]
    [1, 4, 5, 7]
    [1, 4, 5, 8]
    [1, 4, 6, 7]
    [1, 4, 6, 8]
    [2, 3, 5, 7]
    [2, 3, 5, 8]
    [2, 3, 6, 7]
    [2, 3, 6, 8]
    [2, 4, 5, 7]
    [2, 4, 5, 8]
    [2, 4, 6, 7]
    [2, 4, 6, 8]

    :author: Kiran Jonnalagadda `original
             source <http://jace.seacrow.com/archive/2007/02/15/generating-combinations-in-python>`_
    """
  for c in _listcombinations(listoflists):
    yield c

def _listcombinations(listoflists, curlist=[], parents=[]):
    if curlist == []:
        curlist = listoflists[0]
        remlist = listoflists[1:]
    else:
        remlist = listoflists
    for item in curlist:
        if len(remlist) > 0:
            for c in _listcombinations(remlist[1:], remlist[0], parents+[item]):
                yield c
        else:
            yield parents+[item]

def permutations(seq) :
  """ Iterator for all permutations of sequence seq.

  :param seq: Sequence of elements to permute
  :type seq: sequence
  """
  
  if len(seq) <= 1 :
    yield seq
  else:
    for perm in permutations(seq[1:]) :
      for i in range(len(seq)) :
        yield perm[:i] + seq[0:1] + perm[i:]


def allPairs(lst) :
  """ An iterator over all distinct pairs of elements in lst."""
  if len(lst) < 2:
    return
  
  x0 = lst[0]
  for x in lst[1:] :
    yield (x0, x)
  for p in allPairs(lst[1:]) :
    yield p


def lbinomial(k, n, p) :
  """ Log of probability of obtaining k successes out of n from a process with success
  probability p.
  """
  v = gammaln(n+1) - (gammaln(k+1) + gammaln(n-k+1))
  return v + k * math.log(p) + (n-k) * math.log(1-p)


def uniformvec(n,k) :
  """ Partition n to k equal sized integers (as near as possible)"""
  
  b = n / k
  r = n - b * k
  return (b+1,)*r + (b,)*(k-r)

#  return numpy.sum([log(scipy.stats.binom.pmf(x[i], sum(x[i:]),
#  1.0/(len(x)-i))) for i in range(len(x)-1)])

def uniformLike(x) :
  """ The non-normalized probability of getting the *bincont* in x.

  Assuming x was obtained by drawing an index uniformly from [0..n-1]
  (n=len(x)) and incrementing x[i] N times (N=sum(x)), what is the
  probability of obtaining x?

  Useful for comparing likelihood of diffrent x's or to the uniform
  vector from :py:func:`uniformvec`.
  """
  
  return sum([lbinomial(x[i], sum(x[i:]), 1.0/(len(x)-i)) for i in range(len(x)-1)])



