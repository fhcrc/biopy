## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

from __future__ import division

__all__ = ["getArrivalTimes"]

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
