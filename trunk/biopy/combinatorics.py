## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

""" Generic combinatorial utilities.
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
  """ C{n} factorial.

  @rtype: long
  """
  return _prod([long(x) for x in range(1, n+1)])

def choose(n,k) :
  """ Binomial: Number of ways to choose C{k} of C{n}.

  @rtype: long
  """
  return _prod(range(k+1,n+1)) // _prod(range(2,n-k+1))

def lchoose(n,k) :
  """ Log of choose(C{n},C{k}).

  Via L{scipy.special.gammaln} log of the Gamma function.
  """
  return special.gammaln(n+1) - special.gammaln(k+1) - special.gammaln(n-k+1)

def nPairs(k) :
  """ choose(C{k},2) (Convenience).

  @rtype: int
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

    @author: Kiran Jonnalagadda U{original
             source<http://jace.seacrow.com/archive/2007/02/15/generating-combinations-in-python>} 
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

def permutations(t) :
  """ Iterator for all permutations of sequence C{t}.

  @param t: Sequence of elements to permute 
  @type t: sequence
  """
  
  if len(t) <= 1 :
    yield t
  else:
    for perm in permutations(t[1:]) :
      for i in range(len(t)) :
        yield perm[:i] + t[0:1] + perm[i:]


def allPairs(lst) :
  """ An iterator over all distinct pairs of elements in C{lst}.  """
  if len(lst) < 2:
    return
  
  x0 = lst[0]
  for x in lst[1:] :
    yield (x0, x)
  for p in allPairs(lst[1:]) :
    yield p


def lbinomial(k, n, p) :
  """ Log of probability of having C{k} successes out of C{n} with success
  probability C{p}.
  """
  v = gammaln(n+1) - (gammaln(k+1) + gammaln(n-k+1))
  return v + k * math.log(p) + (n-k) * math.log(1-p)


def uniformvec(n,k) :
  """ Partition n to k equal sized integers (as near as possible)

  @param n:
  @type n: int
  @param k:
  @type k: int
  """
  b = n / k
  r = n - b * k
  return (b+1,)*r + (b,)*(k-r)

#  return numpy.sum([log(scipy.stats.binom.pmf(x[i], sum(x[i:]),
#  1.0/(len(x)-i))) for i in range(len(x)-1)])

def uniformLike(x) :
  """ The non-normalized probability of getting the I{bincont} in C{x}.

  Assuming C{x} was obtained by drawing an index uniformly from [0..n-1]
  (C{n=len(x)}) and incrementing C{x[i]} N times (C{N=sum(x)}), what is the
  probability of obtaining C{x}?

  Useful for comparing likelihood of diffrent C{x}'s or to the "uniform"
  vector from L{uniformvec}.
  """
  
  return sum([lbinomial(x[i], sum(x[i:]), 1.0/(len(x)-i)) for i in range(len(x)-1)])



