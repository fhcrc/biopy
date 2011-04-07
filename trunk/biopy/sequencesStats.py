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

def ASD(n, ne, mu, lam) :
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

  return (within, between)
  
