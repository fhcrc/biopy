## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

""" Demographic functions (population size over time).
"""

from __future__ import division

from math import exp, log, sqrt, pi, cos

import random, operator
from scipy.integrate import quad
from numpy import array, arctan, tan

__all__ = ["ConstantPopulation", "StepFunctionPopulation",
           "SinusoidalPopulation", "LinearPiecewisePopulation",
           "ExponentialGrowthPopulation", "ScaledDemographic"]

class Demographic(object) :
  @staticmethod
  def parseDemographic(spec) :
    tp = spec[0].upper()
    if tp == 'C' or tp.isdigit() :
      if tp == 'C' :
        spec = spec[2:]
      return ConstantPopulation(float(spec))
    
    spec = spec[2 if spec[1] == ',' else 1:].split(',')
    vals = [float(x) for x in spec[::2]]
    xvals = [float(x) for x in spec[1::2]]
    if tp == 'S' :
      return StepFunctionPopulation(vals, xvals)
    if tp == 'L' :
      return LinearPiecewisePopulation(vals, xvals)

  def getDT(self, N) :
    """ A suitable value of deltaT to be used as time increment when simulating"""
    return 1.0/N
  
  def NCEbySimulation(self, nLineages, tStart, N = 50000) :
    """ Draw, by simulation, the time until any 2 lineages from the initial
    C{nLineages} coalesce.

    Start at time C{tStart}, and advance by adding small time increments, based
    on number of steps C{N}. For each interval, the probability of no
    coalescence is prob(no event in t,t+dt) = exp(-lambda(t) * dt)
     """
    ko2 = (nLineages*(nLineages-1))//2
    dt = self.getDT(N)
    t = tStart
    
    while True:
      pop = self.population(t + dt/2)
      lam = ko2/pop

      probin = 1 - exp(-lam * dt)
      if random.random() < probin :
        return t + dt
      t += dt

  def NCE_SerialBySimulation(self, lineagesInfo, tStart, N = 50000) :
    """ Draw, by simulation, the time until any 2 lineages from the initial set
    of serial samples in C{lineagesInfo} coalesce.
    """
    lineagesInfo = list(lineagesInfo)
    nLineages,t0 = lineagesInfo.pop(0)

    while len(lineagesInfo) and lineagesInfo[0][1] < tStart :
      n,t0 = lineagesInfo.pop(0)
      nLineages += n
    ko2 = (nLineages*(nLineages-1))//2
    dt = self.getDT(N)
    t = tStart

    while True :
      pop = self.population(t + dt/2)
      lam = ko2/pop

      probin = 1 - exp(-lam * dt)
      if random.random() < probin :
        return t + dt
      t += dt

      while len(lineagesInfo) and lineagesInfo[0][1] < t :
        n,t0 = lineagesInfo.pop(0)
        nLineages += n
        ko2 = (nLineages*(nLineages-1))//2

  def timeToNextCoalescentSerial(self, lineagesInfo, tStart, done = 0) :
    lineagesInfo = list(lineagesInfo)
    nLineages,t0 = lineagesInfo.pop(0)

    while len(lineagesInfo) and lineagesInfo[0][1] < tStart :
      n,t0 = lineagesInfo.pop(0)
      nLineages += n

    assert nLineages > done
    nLineages -= done
    
    t = tStart
    if nLineages == 1:
      n,t = lineagesInfo.pop(0)
      nLineages += n

    assert nLineages > 1
    
    while True :
      t += self.timeToNextCoalescent(nLineages, t)

      if len(lineagesInfo) == 0 or t < lineagesInfo[0][1] :
        return t - tStart
      n,t = lineagesInfo.pop(0)
      nLineages += n
      
  # tests and validations
  
  def numerical(self, low, high) :
    """ Integrate 1/N(t) numerically in the range [low,high]."""
    return quad(lambda x : 1/self.population(x) ,low ,high)[0]

  def test1(self, nTimes = 1000) :
    """ Test difference between numerical integration and exact """

    xRange = self.naturalLimit()
    # If no natural limit pick arbitrary number. Extend limit to test
    # behaviour past boundary point
    xRange = xRange * 1.1 if xRange else 100.0
    
    xs = [random.uniform(0, xRange) for k in range(nTimes)]
    f = lambda x : 1/self.population(x)
    
    return array([self.integrate(high) - quad(f, 0 ,high)[0] for high in xs]).max()

  def graphPoints(self, xmax) :
    raise NotImplementedError
  
  def scale(self, factor) :
    raise NotImplementedError
  
class ConstantPopulation(Demographic) :
  """ Constant population size."""
  def __init__(self, N) :
    self.pop = float(N)
    
  def __str__(self) :
    return "%g over all interval," % self.pop

  def __repr__(self) :
    return "demographic.ConstantPopulation(%g)" % self.pop

  def graphPoints(self, xmax) :
    return ((0, xmax), (self.pop, self.pop))
  
  def population(self, t) :
    return self.pop

  def integrate(self, e) :
    """ Integrate 1/N(t) over 0 <= t <= e """

    return e/self.pop

  def integrateExpression(self, e) :
    """ Integration expression for 1/N(t) over 0 <= t <= e """

    return e + "/" + str(self.pop)

  def timeToNextCoalescent(self, nLineages, t) :
    """ Randomly draw the time interval where any 2 lineages from the initial
    C{nLineages} coalesce, starting from time C{t}."""
    
    return random.expovariate((nLineages * (nLineages-1)) / (2.0 * self.population(t)))

  def naturalLimit(self) :
    return None

  def scale(self, factor) :
    return ConstantPopulation(self.population(0) * factor)    
    
class StepFunctionPopulation(Demographic) :
  """ Stepwise Constant population size."""

  def __init__(self, vals, xvals) :
    """ population = vals[k] over xvals[k-1],xvals[k], where xvals[-1] is implicitly 0"""
    assert len(vals) == len(xvals) + 1
    
    self.vals = [float(x) for x in vals]
    self.xvals = [float(x) for x in xvals]

  def __str__(self) :
    if len(self.xvals) == 0 :
      r = "%g on all interval" % self.vals[0]
    else :
      r = "%g between 0 and %g" % (self.vals[0], self.xvals[0])
      for k in range(1,len(self.xvals)) :
        r = r + (", %g between %g and %g" % (self.vals[k], self.xvals[k-1], self.xvals[k]))
      r = r + (", %g onwards." % self.vals[-1])
      
    return r

  def __repr__(self) :
    return "StepFunctionPopulation(" + repr(self.vals) + "," + repr(self.xvals) + ")";
  
  def graphPoints(self, xmax) :
    # might be problem if xmax is less than max(x) of target
    zbxs = [0,] + self.xvals + [xmax,]
    pxs = reduce(operator.add, zip(zbxs[:-1], zbxs[1:]))
    pys = reduce(operator.add, zip(self.vals, self.vals))
    return (pxs, pys)

  def naturalLimit(self) :
    return max(self.xvals) if len(self.xvals) else None

  def population(self, t) :
    k = 0
    # locate interval
    while k < len(self.xvals) and self.xvals[k] < t :
      k += 1
      
    return self.vals[k]
  
  def integrate(self, xHigh) :
    """ Integrate 1/N(t) over  0 <= t <= xHigh """
    x = 0
    k = 0
    v = 0.0
    while x < xHigh :
      if k == len(self.xvals) or xHigh <= self.xvals[k] :
        v += (xHigh - x) / self.vals[k]
        x = xHigh
      else :
        v += (self.xvals[k] - x) / self.vals[k]
        x = self.xvals[k]
        k += 1
    return v

  def timeToNextCoalescent(self, nLineages, t) :
    """ Time to next coalescent event for nLineages lineages starting at t. """
    
    # find interval of t
    k = 0
    while k < len(self.xvals) and self.xvals[k] < t :
      k += 1

    # amount to deduct from each interval to set starting point at t
    if k > 0 :
      t0 = t - self.xvals[k-1]
    else :
      t0 = t
    xvals = [x - t0 for x in self.xvals[k:]]
    vals = self.vals[k:]

    #print t, xvals, vals
    ko2 = (nLineages*(nLineages-1))/2.0

    return StepFunctionPopulation(vals, xvals).solve(random.random(), ko2)
    
  def solve(self, v, ko2) :
    ival = -log(1-v)/ko2

    k = 0
    x = 0
    if len(self.xvals) :
      inc = (self.xvals[k] - x) * (1.0/self.vals[k])
      while ival > inc :
        ival -= inc
        x = self.xvals[k]
        k += 1
        if k == len(self.xvals) :
          break
        inc = (self.xvals[k] - x) * (1.0/self.vals[k])

    x += ival/(1.0/self.vals[k])
    return x

  def scale(self, factor) :
    yv = [factor * x for x in self.vals]
    return StepFunctionPopulation(yv, self.xvals)    

#from scipy import weave
from cchelp import demoLPintegrate, demoLPpopulation

class LinearPiecewisePopulation(Demographic) :
  # population = vals[k] at xvals[k-1] and linear between.
  # xvals[-1] is implicitly 0
  def __init__(self, vals, xvals) :
    assert len(vals) == len(xvals) + 1
    self.vals = [float(x) for x in vals]
    self.xvals = [float(x) for x in xvals]

  def __str__(self) :
    r = ""
    # "from %g to %g on [0 - %g]" % (self.vals[0], self.xvals[0])
    for k in range(0, len(self.xvals)) :
      if k > 0 :
        x0 = self.xvals[k-1]
      else :
        x0 = 0
      r = r + ("from %g to %g on [%g,%g], " % \
               (self.vals[k], self.vals[k+1], x0, self.xvals[k]))
    r = r + ("%g onwards." % self.vals[-1])
    return r

  def __repr__(self) :
    r = "LinearPiecewisePopulation(" + \
        repr(self.vals) + "," + repr(self.xvals) + ")"
    return r

  def graphPoints(self, xmax) :
    pxs = [0,] + self.xvals
    if xmax > max(pxs) :
      pxs = pxs + [xmax,]
    pys = [self.population(x) for x in pxs]
    return (pxs, pys)

  def naturalLimit(self) :
    return max(self.xvals) if len(self.xvals) else None

  def populationPython(self, t) :
    k = 0
    while k < len(self.xvals) and self.xvals[k] < t :
      k += 1

    if k == len(self.xvals) :
      return self.vals[k]
    
    # make t dt relative to start
    if k > 0 :
      t -= self.xvals[k-1]
      width = self.xvals[k] - self.xvals[k-1]
    else :
      width = self.xvals[k]
      
    return self.vals[k] + (t/width) * (self.vals[k+1] - self.vals[k])

  def population(self, t) :
    return demoLPpopulation(self.vals, self.xvals, t)
  
  # Integrate 1/N from 0 to xHigh
  def integratePython(self, xHigh) :
    x = 0.0
    k = 0
    v = 0.0
    while x < xHigh :
      if k == len(self.xvals) :
        v += (xHigh - x) * (1.0/self.vals[k])
        x = xHigh
      else :
        pop0 = self.vals[k]
        pop1 = self.vals[k+1]
        x1 = self.xvals[k] 
        dx = x1 - x
        
        if xHigh < x1 :
          ndx = xHigh-x
          pop1 = pop0 + (ndx/dx)*(pop1-pop0)
          dx = ndx
          
        if pop0 == pop1 :
          v += dx/pop0
        else :
          m = dx / (pop1 - pop0)
          v += m * log(pop1/pop0)

        x = x1
        k += 1

    return v


  def integrate(self, xHigh) :
    return demoLPintegrate(self.vals, self.xvals, xHigh)
    
  def timeToNextCoalescent(self, nl, t) :
    ko2 = (nl*(nl-1))//2

    k = 0
    while k < len(self.xvals) and self.xvals[k] < t :
      k += 1

    popt = self.population(t)
    if k > 0 :
      t0 = t - self.xvals[k-1]
    else :
      t0 = t
    xvals = [x - t0 for x in self.xvals[k:]]
    vals = self.vals[k:]
    vals[0] = popt

    if len(xvals) == 1 and xvals[0] == 0.0 :
      xvals = []
      vals = vals[:1]

    return LinearPiecewisePopulation(vals, xvals).solve(random.random(), ko2)
    
  def solve(self, v, ko2) :
    ival = -log(1-v)/ko2

    k = 0
    x = 0
    while 1 :
      if k == len(self.xvals) :
        break
      x1 = self.xvals[k]
      pop0 = self.vals[k]
      pop1 = self.vals[k+1]

      if pop0 != pop1 :
        inc = ((x1 - x) / (pop1 - pop0)) * log(pop1/pop0)
      else :
        inc = (x1 - x) / pop0
      
      if ival <= inc :
        break
      else :
        ival -= inc
        x = x1
        k += 1

    if k == len(self.xvals) :
      x += ival/(1.0/self.vals[k])
    else :
      x1 = self.xvals[k]
      pop0 = self.vals[k]
      pop1 = self.vals[k+1]
      if pop0 != pop1 :
        dx = x1 - x
        a = (exp(ival * (pop1-pop0)/dx) - 1) * (pop0/(pop1-pop0))

        x += a * dx
      else:
        x += pop0 * ival

    return x

  def scale(self, factor) :
    yv = [factor * x for x in self.vals]
    return LinearPiecewisePopulation(yv, self.xvals)    


class SinusoidalPopulation(Demographic) :
  """ Sinusoidal between amplitude and amplitude*base (both positive)
  N(t) = amp * ((1-base) * ((1+cos(t0 + t * x))/2) + base)
  max is on -t0/tscale, min on (pi-t0)/tscale """
  
  def __init__(self, amplitude, base, tscale = 1, t0 = 0) :
    
    self.amplitude = float(amplitude)
    self.base = float(base)
    self.tscale = float(tscale)
    self.t0 = float(t0)

    self.cycle = 2 * pi / tscale
    self.offset = self.t0 / self.tscale

    self.highlim = -self.offset + self.cycle/2

    self.ion0 = self.integral(0)
    i1 = self.integral(self.highlim)
    self.ifullcycle =  i1 - self.integral(-self.offset + -self.cycle/2)
    self.izerotomin = i1 - self.ion0

    self.ionxstart = self.integral(-self.offset - self.cycle/2)

  def __repr__(self) :
    return "SinusoidalPopulation(%g %g %g %g)" \
           % (self.amplitude, self.base, self.tscale, self.t0)

  def graphPoints(self, xmax, n = 100) :
    dt = self.cycle / n
    np = int(xmax/dt)
    xmax = float(xmax)
    x = [(xmax * k) / np for k in range(np+1)]
    return (x,[self.population(z) for z in x])

  def population(self, t) :
    x = self.t0 + self.tscale * t
    return self.amplitude * (((1 - self.base) * (1+cos(x))/2.0) + self.base)

  def integral(self, x) :
    """ Indefinite, valid in [(-pi-t0)/scale,(pi-t0)/t0]"""
    
    a = self.base
    b = self.tscale
    t0 = self.t0
    
    assert -pi <= x*b + t0 <= pi, (x,b,t0)
    
    c = self.amplitude

    a1 = sqrt(a)

    # Thanks Wolfram! (http://integrals.wolfram.com/index.jsp)
    return 2 * arctan(a1 * tan(0.5 * (t0 + b * x)))/(a1 * b * c)

  def integrate(self, x) :
    """ value of integral between 0 and x """
    assert x >= 0

    b = self.tscale

    # break into full cycles and remainder
    iturns = int(x / self.cycle)
    x = x - iturns * self.cycle

    if x <= self.highlim :
      v = self.integral(x) - self.ion0
    else :
      v = self.izerotomin + self.integral(x - self.cycle) - self.ionxstart

    return v + iturns * self.ifullcycle

    
  def timeToNextCoalescent(self, nl, t) :
    return self.solve(random.random(), nl*(nl-1)/2, t) - t

  def solve(self, v, ko2, t = 0) :
    """ x so that v = prob(wait time for coalescent starting at t >= x | (k over 2), N).
    Technically v = 1 - exp( (k over 2) int_0^x 1/N(t) dt
    """
    # invIntegral solves from 0. Add excess now
    v0 = self.integrate(t)
    ival = v0 + -log(1-v)/ko2

    return self.invIntegral(ival)
  
  def invIntegral(self, ival) :
    """ x so that integrate(x) = ival """
    ifullcycle = self.ifullcycle
    nturns = int(ival / ifullcycle)

    ival -= nturns * ifullcycle

    a = self.base
    b = self.tscale
    c = self.amplitude
    t0 = self.t0
    a1 = sqrt(a)
    
    if ival <= self.izerotomin :
      ival += self.ion0
      x = 0
    else :
      ival += -self.izerotomin + self.ionxstart
      x = self.cycle

    x += (arctan(tan(ival * ((a1 * b * c)/ 2)) / a1) * 2 - t0) / b
    return x + nturns * self.cycle
      

  def getDT(self, N) :
    return self.cycle/N

  def numerical(self, low, high) :
    """ integrate 1/N(t) numerically between low and high"""
    return quad(lambda x : 1/self.population(x) ,low ,high)[0]


  def test1(self,  nTimes = 1000) :
    """ test difference between numerical integration and exact inside main
    cycle interval. (inside -pi/b pi/b) """
    #assert (-math.pi/self.tscale) <= low <= high <= (math.pi/self.tscale)
    lim = self.cycle / 2
    t0 = self.t0/self.tscale
    
    xlow = [random.uniform(-lim-t0, lim-t0) for k in range(nTimes)]
    xhigh = [random.uniform(xl, lim-t0) for xl in xlow]
    
    return array([self.numerical(l, h)  -
                  (self.integral(h) - self.integral(l)) for l,h in zip(xlow,xhigh)]).max()

  
  def test2(self, nTimes = 1000, xRange = None) :
    """ test numerical vs exact integration [0,x] for x in xRange=10 cycles """
    if xRange is None :
      xRange = 10 * self.cycle

    xs = [random.uniform(0, xRange) for k in range(nTimes)]
    f = lambda x : 1/self.population(x)
    
    return array([self.integrate(high) - quad(f, 0 ,high)[0] for high in xs]).max()

  def test3(self, nTimes = 1000, xRange = None) :
    """ test solving """
    if xRange is None :
      xRange = 10 * self.cycle
      
    xs = [random.uniform(0, xRange) for k in range(nTimes)]
    return array([abs(self.invIntegral(self.integrate(r)) - r) for r in xs]).max()


class ExponentialGrowthPopulation(Demographic) :
  """ N0 * exp(-rt)"""
  def __init__(self, N, rate) :
    self.pop0 = float(N)
    self.rate = rate
    
  def __repr__(self) :
    return "%g exp(-%g t)" % (self.pop0,self.rate)

  def graphPoints(self, xmax) :
    N = 100
    xs = [k * xmax/100.0 for k in range(N)]
    return (xs, [self.population(x) for x in xs])
  
  def population(self, t) :
    return self.pop0 * exp(-self.rate * t)

  def integrate(self, e) :
    """ Integrate 1/N(t) for  0 <= t <= e """

    return ( exp(self.rate * e) - 1 ) / (self.pop0 * self.rate)

  def timeToNextCoalescent(self, nLineages, t) :
    """ time to next coalescent event for nLineages lineages starting at t. """
    ko2 = (nLineages*(nLineages-1))/2.0
    ival = -log(1-random.random())/ko2

    # shift of x axis by t is the same as setting the value at the new "0"
    p0 = self.population(t)
    e = log(1 + ival * p0 * self.rate)/self.rate

    return e

  def naturalLimit(self) :
    return None

  def scale(self, factor) :
    return ExponentialGrowthPopulation(self.population(0) * factor, self.rate)    

class ScaledDemographic(object) :
  def __init__(self, dfunc, factor) :
    """ a Demographic function equal to dfunc * factor"""
    self.demo = dfunc
    self.factor = factor
    
  def population(self, t) :
    return self.demo.population(t) * self.factor

  def integrate(self, t) :
    return self.demo.integrate(t) / self.factor

