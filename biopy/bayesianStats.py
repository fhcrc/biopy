## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

"""
===================
Bayesian statistics
===================

Calculate Bayesian statistics such as Highest Posterior Density (HPD) and
Effective Sample Size (BEAST interpretation).
"""

from __future__ import division

__all__ = ["hpd", "effectiveSampleSize"]

def hpd(data, level) :
  """ The Highest Posterior Density (credible) interval of data at level level.

  :param data: sequence of real values
  :param level: (0 < level < 1)
  """ 
  
  d = list(data)
  d.sort()

  nData = len(data)
  nIn = int(round(level * nData))
  if nIn < 2 :
    raise RuntimeError("not enough data")
  
  i = 0
  r = d[i+nIn-1] - d[i]
  for k in range(len(d) - (nIn - 1)) :
    rk = d[k+nIn-1] - d[k]
    if rk < r :
      r = rk
      i = k

  assert 0 <= i <= i+nIn-1 < len(d)
  
  return (d[i], d[i+nIn-1])

from cchelp import effectiveSampleStep

def effectiveSampleSize(data) :
  return len(data)/effectiveSampleStep(data)[0]

# import numpy

## def effectiveSampleSize(data, stepSize = 1) :
##   """ Effective sample size, as computed by BEAST Tracer.

##   :param data: sequence of real values
##   """
  
##   samples = len(data)

##   assert len(data) > 1,"no stats for short sequences"
  
##   maxLag = min(samples//3, 1000)

##   gammaStat = [0,]*maxLag
##   #varGammaStat = [0,]*maxLag

##   varStat = 0.0;

##   if type(data) != numpy.ndarray :
##     data = numpy.array(data)

##   normalizedData = data - data.mean()
  
##   for lag in range(maxLag) :
##     v1 = normalizedData[:samples-lag]
##     v2 = normalizedData[lag:]
##     v = v1 * v2
##     gammaStat[lag] = sum(v) / len(v)
##     #varGammaStat[lag] = sum(v*v) / len(v)
##     #varGammaStat[lag] -= gammaStat[0] ** 2

##     # print lag, gammaStat[lag], varGammaStat[lag]
    
##     if lag == 0 :
##       varStat = gammaStat[0]
##     elif lag % 2 == 0 :
##       s = gammaStat[lag-1] + gammaStat[lag]
##       if s > 0 :
##          varStat += 2.0*s
##       else :
##         break
      
##   # standard error of mean
##   # stdErrorOfMean = Math.sqrt(varStat/samples);

##   # auto correlation time
##   act = stepSize * varStat / gammaStat[0]

##   # effective sample size
##   ess = (stepSize * samples) / act

##   return ess
