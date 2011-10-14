## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

from __future__ import division

import random

from treeutils import TreeBuilder

__all__ = ["getArrivalTimes", "sampleCoalescentTree"]

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
  
  if isinstance(labels[0], str) :
    n = len(labels)
    times= [0,]*n
    lineageInfo = len(labels)
  else :
    raise ""
    if tipNameFunc is None :
      tot = sum([n for n,t in lineageInfo])
      labels = ["tip%03d" % k for k in range(tot)]
    else :
      labels = reduce(operator.add, [tipNameFunc(l) for l in lineageInfo])
      
    times = reduce(operator.add,  [[t,]*n for n,t in lineageInfo])

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
