## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

""" 
"""

from __future__ import division

from numpy import prod

__all__ = ["ASD"]

def genBinomial(x, k) :
  """ Generalized binomial """
  return prod([(x-i)/(k-i) for i in range(k)])

def ASD(n, ne, mu, lam, withVariance = False) :
  """ Average Sequence Dissimilarity

  Return a pair of (within,between) ASD's for C{n} species, Yule species tree
  with birth rate of C{lam}, multispecies coalescence with constant population
  size C{ne} throughout, and mutation rate C{mu}.
  """
  
  mune = mu*ne
  dm = (8*mune+3)
  within = 6*mune/dm

  w = (8*mu)/(3*lam)
  between = within + (9*(n+1)/(4*(n-1)*(w-1)*dm)) * \
            ((w+1)/genBinomial(w + n, n) + ((w *  (n-1))/(n+1)) - 1)

  if withVariance :
    dm1 = (3 + 16 * mune)
    varWithin = within**2 * (3/dm1)
    within = (within, varWithin)

    v1 = (12 - 6*n*(n+1)/genBinomial(w + n, n-1)) / (dm * (w-1))
    v2 = (6 - 3*n*(n+1)/genBinomial(2*w + n, n-1)) / (dm1 * (2*w-1))
    ex2 = 9/16 * ( 1 -  (v1 - v2)/(n-1))
    print ex2

    between = (between, ex2 - between**2)
    
  return (within, between)
  
