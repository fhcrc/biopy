## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 
""" Birth-Death process utilities.

  Most math is taken from Gerhard (Stadler) U{thesis,
  http://deposit.ddb.de/cgi-bin/dokserv?idn=992078377&dok_var=d1&dok_ext=pdf&filename=992078377.pdf}.
"""

from __future__ import division
from math import log, exp

from random import expovariate
from numpy import cumsum
import scipy.special

from combinatorics import choose

__all__ = ["yuleTimes",
           "yuleHeights",
           "yuleHeightsConditionalOnRoot",
           "yuleHeightsLogLiklihood",
           "yuleHeightsCondRootLogLiklihood",
           "drawYuleTree",
           "drawBDTree"]
           
def yuleTimes(birthRate, n) :
  times = []

  for k in range(2, n+1) :
    times.append(expovariate(k * birthRate))

  return times
    
def yuleHeights(birthRate, n) :
  """ Draw a sequence of node heights of a binary tree with C{n} leaves under a
  birth process with birth rate C{birthRate}.

  The list contains C{n-1} sorted times, root height being the last.
  """
  
  times = yuleTimes(birthRate, n)

  return cumsum(list(reversed(times)))

# Taken from Tanja standler thesis, page 99-100
def _BDrootQ(r, lam, mu, n) :
  if lam == mu :
    return 1/(lam*(r**(-1/n)-1))
  return (1/(lam-mu)) * log((1 - (mu/lam) * r**(1/n))/(1 - r**(1/n)))

def _BDinternalQ(s, t, lam, mu, n) :
  if lam == mu :
    return s*t/(1 + lam*t*(1-s))
  h = lam-mu
  e1 = exp(-h*t)
  a1 = lam - mu*e1
  a2 = s * (1-e1)
  return (1/h) * log((a1 - mu*a2)/(a1-lam*a2))

def BDheights(lam, mu, n) :
  """ Draw C{n} node heights from a Birth-Death process with birth rate C{lam}
  and death rate C{mu}.

  Return sorted list of times."""
  
  assert 0 <= mu <= lam and lam > 0
  
  if mu == 0 :
    return yuleHeights(lam, n)
  
  r0 = random.random()
  t0 = _BDrootQ(r0, lam, mu, n)
  s = [_BDinternalQ(random.random(), t0, lam, mu, n) for k in range(n-1)]
  return sorted(s)

def BDtimes(lam, mu, n) :
  h = BDheights(lam, mu, n)
  return [h[k] - h[k-1] for k in range(n-2,0,-1)] + [h[0]]

from decimal import Decimal, getcontext
getcontext().prec = 108

def _f1(r,k,i) :
  return (1/(1-r)).ln() - \
         sum([(1-(1/(1-r))**j)*choose(k+i,j)*(-1)**j/j for j in range(1,k+i+1)])

def BDexpected(k,n,lam,mu) :
  """ Expected k'th speciation time for BD tree conditioned on n taxa.
  k == 1 is root. Have serious numerical issues for small mu, so using
  'Decimal', which unfortuanetly requires changing code :(
  """
  if mu == 0 :
    return sum([1/i for i in range(k+1, n+1)])/lam
  if mu == lam:
    return (n-k)/(lam*k)

  mu,lam = Decimal(str(mu)), Decimal(str(lam))
  r = mu/lam
  s1 = sum([_f1(r,k,i) * choose(n-k-1, i)*(1/r -1)**(k+i)/((k+i+1)*r) for i in range(n-k)])

  return float(Decimal(str(((k+1)/lam) * choose(n, k+1) * (-1)**k)) * s1)

from treeutils import TreeBuilder
import random

def treeFromTimes(times, order = None) :
  """ Build a biopython tree from intra-coalescent intervals in
  C{times}.

  Without an order, the joining is random. With C{order}, use the given ordering
  supplied as a sequence of pairs, each pair is two indices picked to join."""
  
  assert order is None or len(times) == len(order)

  n = len(times) + 1
  
  tb = TreeBuilder()
  taxa = [[tb.createLeaf("ID%d" % k), 0.0] for k in range(n)]
  height = 0.0
  for t in reversed(times) :
    i,j = order.pop(0) if order else \
          random.sample(range(len(taxa)), 2)
    l1,l2 = taxa[i],taxa[j]
    height += t
    nd = tb.mergeNodes(l1[0], height - l1[1], l2[0], height - l2[1])
    taxa.pop(max(i,j))
    taxa[min(i,j)] = [nd, height]

  return tb.finalize(taxa[0][0])

def drawYuleTree(birthRate, n, order = None) :
  times = yuleTimes(birthRate, n)
  return treeFromTimes(times, order)

def drawBDTree(birthRate, deathRate, n) :
  times = BDtimes(birthRate, deathRate, n)
  return treeFromTimes(times)

def _tanja08loglike(nodeHeights, birthRate, deathRate) :
  """ nodeHeights, internal only node heights, root is first"""
  x = nodeHeights

  b,d = birthRate,deathRate
  r = b-d
  a = d/b

  n = len(x) + 1
  
  e1 = -r*x[0]
  ee1 = exp(e1)
  
  p1 = scipy.special.gammaln(n+1) + e1 - log(1 - a*ee1) + \
      log(1-a)*(n) + (n-1) * log(r)
  p2 = sum( [-r*h - 2*log(1 - a * exp(-r*h)) for h in x])
  
  return p1 + p2

def yuleHeightsLogLiklihood(h, birthRate) :
  """ Log-Likelihood for sorted node heights list C{h} under a yule process with
  birth rate C{birthRate}.
  """
  
  return _tanja08loglike(list(reversed(h)), birthRate, 0)

def yuleTimesCondRoot(birthRate, tot, n, maxTries = 10000) :
  nTries = 0
  while nTries < maxTries :
    nTries += 1
    
    times = [None]

    for k in range(1, n-1) :
      times.append(expovariate(k * birthRate))

    t = sum(times[1:])
    if t < tot :
      times[0] = tot - t
      return times

  raise

def yuleHeightsConditionalOnRoot(birthRate, rootTime, n) :
  """ Draw a sequence of node heights of a binary tree with C{n} leaves under a
  birth process with birth rate C{birthRate}, Conditional on root height being
  C{rootTime}. 

  The list contains C{n-1} sorted times, root height being the last.
  """
  
  times = yuleTimesCondRoot(birthRate, rootTime, n)

  return cumsum(list(reversed(times)))


def yuleHeightsCondRootLogLiklihood(h, birthRate, normed=True) :
  """ Log-Likelihood for sorted node heights list 'h' under a yule process with
  birth rate 'birthRate', conditional on root height being h[-1].

  With not 'normed', the likelihood is unnormalized, that is up to a constant
  multiplier.
  """
  
  c0 = 0
  n = len(h)
  if normed:
    c0 = n*log(2.0) - log(n) - scipy.special.gammaln(n+2)
    
  rootTime = h[-1]
  c1 = log(birthRate / (1 - exp(-birthRate * rootTime)))
  
  ll = -birthRate * sum(h[:-1]) + (n-1) * c1 + c0
  return ll
