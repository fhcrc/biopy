## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

""" Sample from the coalescent and calculate the likelihood of a tree under the
coalescent."""

from __future__ import division

import random
from math import log

from treeutils import TreeBuilder

__all__ = ["getArrivalTimes", "sampleCoalescentTree", "coalLogLike"]

def getArrivalTimes(demog, lineageInfo) :
  """ Get coalescent arrival times. When lineageInfo is an integer (n), get
  n-1 consecutive times which represent coalescent times from n lineages to 1.
  Those are equivalent to arrival times of the non-homogeneous Poisson
  process with rate function 1/demog(t).

  Otherwise lineageInfo is a sequence of (nSamp,t) pairs, where nSamp samples
  are added at time t.
  """

  times = []
  if isinstance(lineageInfo, int) :
    t = 0
    for k in range(lineageInfo, 1, -1) :
      t += demog.timeToNextCoalescent(k, t)
      times.append(t)
    return times

  # Serial samples
  tot = sum([n for n,t in lineageInfo])
  t = 0
  while len(times) + 1 < tot :
    t += demog.timeToNextCoalescentSerial(lineageInfo, t, done = len(times))
    times.append(t)
  return times

def sampleCoalescentTree(demog, labels) :
  """ Build a tree under the coalescent model using given demographic."""
  
  if isinstance(labels, int) :
    labels = ["tip%03d" % k for k in range(int(labels))]
  else :
    assert all([isinstance(s, str) for s in labels])

  n = len(labels)
  times= [0,]*n
  lineageInfo = len(labels)
      
  # times = reduce(operator.add,  [[t,]*n for n,t in lineageInfo])

  #assert len(labels) == len(times)

  at = getArrivalTimes(demog, lineageInfo)

  tb = TreeBuilder()
  nodes = [tb.createLeaf(l) for l in labels]
  strees = zip(nodes, times)

  for t in at :
    # find 2 labels with times < t and coalesce
    available = [k for k in range(len(strees)) if strees[k][1] < t]

    assert len(available) > 1,(available,t,at,strees)
    
    i,j = random.sample(available, 2)

    l1,l2 = strees[i],strees[j]
    n = tb.mergeNodes([[x[0], t - x[1]] for x in (l1,l2)])
  
    strees.pop(max(i,j))
    strees.pop(min(i,j))
    strees.insert(0, (n, t))

  assert len(strees) == 1
  return tb.finalize(strees[0][0])


def coalLogLike(demog, times, condOnTree = False) :
  """ Log-Likelihood of coalescent times with demographic function demog.
  times is a sorted list of (time,coal), where time goes backwards from zero
  and coal is false for tips, true for a coalescent."""

  like = 0.0
  nl = 0

  lastDemoIntegral = 0
  for time,coal in times :
    if coal :
      if not nl > 1:
        raise RuntimeError("Error in time sequence")
      
      f = demog.integrate(time)
      interval = f - lastDemoIntegral
      lastDemoIntegral = f

      nl2 = (nl * (nl-1))/2
      like -= nl2 * interval

      pop = demog.population(time)
      like += log((1 if condOnTree else nl2) / pop)

      nl -= 1
    else :
      nl += 1

  if nl != 1 :
    raise RuntimeError("Error in time sequence")

  return like


def ologlike(times, demog) :
  """ Log-Likelihood of coalescent times given demographic function demog.
  When times is a list of length 2 with second element a list, it is assumed to
  be serial data tip information."""

  #verbose = 1
  
  like = 0.0
  if len(times) == 2 and hasattr(times[1], '__iter__') :
    times, lineagInfo = times
    
    nl = list(lineagInfo)
    nLineages,t= nl.pop(0)
    assert t == 0.0

    # integral of 1/demog(t) between at [0,prevTime]
    lastDemoIntegral = 0

    # allow a single sample to start with
    if nLineages == 1 :
      n1,t = nl.pop(0)
      nLineages += n1
      lastDemoIntegral = demog.integrate(t)
      
    for time in times :
      while len(nl) and time >= nl[0][1] :
        n1,t = nl.pop(0)
        # no coalescent between previous time and t
        f = demog.integrate(t)
        interval = f - lastDemoIntegral
        lastDemoIntegral = f

        nLineageOver2 = (nLineages * (nLineages-1))/2
        like += -nLineageOver2 * interval

        # prevtime = t
        nLineages += n1
        
      assert nLineages > 1,time
      
      f = demog.integrate(time)
      interval = f - lastDemoIntegral
      lastDemoIntegral = f

      #if verbose: print "nl", nLineages, "interval", interval,

      nLineageOver2 = (nLineages * (nLineages-1))/2
      like += -nLineageOver2 * interval

      pop = demog.population(time)
      like += log(1.0 / pop) ; # log(nLineageOver2 / pop)

      nLineages -= 1
      #if verbose: print "pop", pop, "like", like  
  else :
    nLineages = len(times)+1

    lastDemoIntegral = 0

    for time in times :
      f = demog.integrate(time)
      interval = f - lastDemoIntegral
      lastDemoIntegral = f

      #if verbose: print "nl", nLineages, "time", time, "i", time - ltime, "interval", interval,
      #ltime = time

      nLineageOver2 = (nLineages * (nLineages-1))/2
      like += -nLineageOver2 * interval

      pop = demog.population(time)
      like += log(nLineageOver2 / pop)

      nLineages -= 1

      #if verbose: print "pop", pop, "like", like

  assert nLineages == 1
  return like
