## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

"""
====================
Random Distributions
====================
"""

from __future__ import division

from math import log, exp, sqrt, pi
import random
from scipy import special

__all__ = ["Uniform", "Exponential", "Gamma", "Poisson", "InvGamma",
           "LogNormal", "Normal",  "Delta", "parseDistribution"]

class Uniform(object) :
  def __init__(self, l, h) :
    self.low = float(l)
    self.high = float(h)
    self.p = 1 / (h - l)
    self.lp = log(self.p)
    
  def __repr__(self) :
    return "Uniform(%g,%g)" % (self.low,self.high)
  
  def pdf(self, x) :
    if self.low <= x <= self.high :
      return self.p

    return 0

  def logpdf(self, x) :
    if self.low <= x <= self.high :
      return self.lp

    raise OverflowError
    # 0 causes a OverflowError

  def domain(self) :
    return (self.low, self.high)
  
  def sample(self) :
    return random.uniform(self.low, self.high)
  
class Exponential(object) :
  def __init__(self, lam) :
    self.lam= lam
    self.lglam = log(self.lam)
    
  def __str__(self) :
    return "Exponential(1/%g)" % (1/self.lam)

  def __repr__(self) :
    return "Exponential(%g)" % self.lam

  def domain(self) :
    return (0, float('inf'))

  def set(self, lam) :
    self.__init__(1/lam)
    
  def pdf(self, x) :
    return self.lam * exp(-self.lam * x)

  def logpdf(self, x) :
    return -self.lam * x + self.lglam
  
  def sample(self) :
    return random.expovariate(self.lam)

class Gamma(object) :
  """ Gamma(shape, scale) mean = shape*scale var shape*scale^2"""
  def __init__(self, shape, scale) :
    self.shape = float(shape)
    self.scale = float(scale)
    self.shapetimeslscale = self.shape * log(self.scale)
    self.gammashape = special.gammaln(self.shape)
    self.denom = self.shapetimeslscale + self.gammashape
    
  def __repr__(self) :
    return "Gamma(%g,%g)" % (self.shape,self.scale)
  
  def domain(self) :
    return (0, float('inf'))

  def pdf(self, x) :
    if x == 0 : return 0
    return exp(self.logpdf(x))

  def logpdf(self, x) :
    return (self.shape-1) * log(x) - (x/self.scale) - self.denom
  
  def sample(self) :
    return random.gammavariate(self.shape, self.scale)
  
class Poisson(object) :
  def __init__(self, lam) :
    self.lam= lam
    self.loglam = log(self.lam) 

  def __repr__(self) :
    return "Poisson(%g)" % self.lam
  
  def pdf(self, k) :
    return exp(-self.lam) * self.lam**k / special.gamma(k+1)

  def logpdf(self, k) :
    return -self.lam + k*self.loglam - special.gammaln(k+1)

class Jeffeys(object) :
  def __init__(self) :
    pass

  def __str__(self) :
    return "Jeffreys"
  
  def __repr__(self) :
    return "dists.Jeffreys()"

  def domain(self) :
    return (0, float('inf'))
  
  def pdf(self, x) :
    return 1.0/x

  def logpdf(self, x) :
    return -log(x)

from scipy.stats import invgamma,lognorm
from scipy.special import gamma

class InvGamma(object) :
  def __init__(self, alpha, beta) :
    """ """
    self.alpha = alpha
    self.beta = beta
    self.dist = invgamma(alpha, scale = beta)
    self.co = log((beta**alpha)/ gamma(alpha))
    
  def __str__(self) :
    return self.__repr__()
    
  def __repr__(self) :
    return "invgamma(%g,%g)" % (self.alpha,self.beta)

  def domain(self) :
    return (0, float('inf'))

  def pdf1(self, x) :
    return self.dist.pdf(x)

  def pdf(self, x) :
    beta,alpha = self.beta, self.alpha
    return ((beta**alpha)/ gamma(alpha)) * (x**(-alpha-1) * exp(-beta/x))

  def _logpdf(self, x) :
    return log(self.pdf(x))

  def logpdf(self, x) :
    return self.co + (-self.alpha-1) * log(x) -self.beta/x

  def sample(self) :
    return self.dist.rvs()
  
class LogNormal(object) :
  def __init__(self, mirs, sil) :
    """ Mean in real Space , Std in log space"""
    self.meanInRealSpace = mirs
    self.stdInLog = sil
    s = mirs * exp(-sil**2/2.)
    self.dist = lognorm(sil, scale = s)
    
  def __str__(self) :
    return self.__repr__()
    
  def __repr__(self) :
    return "LogNormal(%g,%g)" % (self.meanInRealSpace,self.stdInLog)

  def domain(self) :
    return (0, float('inf'))

  def pdf(self, x) :
    return self.dist.pdf(x)

  def pdf1(self, x) :
    s =  self.stdInLog
    scale = self.dist.kwds['scale']
    return 1/(s*x*sqrt(2*pi)) * exp(-1/2*((log(x)-log(scale))/s)**2)

  def __logpdf(self, x) :
    return log(self.pdf(x))

  def logpdf(self, x) :
    s =  self.stdInLog
    scale = self.dist.kwds['scale']
    return -log(s*sqrt(2*pi)) -log(x) - ((log(x/scale)/s)**2)/2

  def sample(self) :
    return float(self.dist.rvs())
  

# lognormal
# scipy.stats.lognorm.rvs(x, scale=exp(y))
#  y == E scipy.mean(scipy.log(vv))
#  x == E scipy.std(scipy.log(vv))
#  exp(y + x**2/2.) ==  E scipy.mean(vv)
#  scipy.stats.lognorm.rvs(x, scale=z*exp(-x*x/2.)) has mean z in real space

from scipy.stats import norm

class Normal(object) :
  def __init__(self, m, s) :
    self.dist = norm(m, s)
    self.m = m
    self.s = s

  def domain(self) :
    return (-float('inf'), float('inf'))
    
  def pdf(self, x) :
    return self.dist.pdf(x)

  def logpdf(self, x) :
    return -log(self.s) - 0.5 * log(2*pi) - ((x - self.m)/self.s)**2/2.0;

  def sample(self) :
    return float(self.dist.rvs())

class Delta(object) :
  def __init__(self, v) :
    self.value = v

  def __repr__(self) :
    return "Delta(%g)" % self.value
  
  def domain(self) :
    return (self.value, self.value)
    
  def pdf(self, x) :
    return 1 if self.value == x else 0

  def logpdf(self, x) :
    return 0 if self.value == x else -float('inf')

  def sample(self) :
    return self.value

def parseDistribution(txt) :
  parts = txt.split(",")

  if parts[0].lower() == 'u' :
    if len(parts) == 3 :
      return Uniform(float(parts[1]), float(parts[2]))
  elif parts[0].lower() == 'e' :
    if len(parts) == 2 :
      return Exponential(float(parts[1]))
  elif parts[0].lower() == 'l' :
    if len(parts) == 3 :
      return LogNormal(float(parts[1]), float(parts[2]))
  elif parts[0].lower() == 'g' :
    if len(parts) == 3 :
      return Gamma(float(parts[1]), float(parts[2]))
  elif parts[0].lower() == 'i' :
    if len(parts) == 3 :
      return InvGamma(float(parts[1]), float(parts[2]))
  elif parts[0].lower() == 'n' :
    if len(parts) == 3 :
      return Normal(float(parts[1]), float(parts[2]))
  elif len(parts) == 1 :
    return Delta(float(parts[0]))
  
  raise RuntimeError("malformed distribution specification:" + txt)
