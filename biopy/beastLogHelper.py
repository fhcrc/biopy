## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

""" Handling of BEAST log and tree files."""

from __future__ import division
import sys, re

__all__ =  ["readTraces", "setDemographics"]

def readTraces(beastFile, traces, report = False) :
  """ Read traces from a BEAST log file.

  Return a sequence of traces, each a pair (name, sequence-of-values).
  
  @param beastFile: BEAST log file name
  @type beastFile: str
  @param traces: A trace name, or a sequence of trace names.
  """
  values = []

  if isinstance(beastFile, str) :
    beastFile = file(beastFile)
    
  if isinstance(traces, str) :
    traces = [traces,]

  assert isinstance(traces, (list, tuple))

  iTraces = []

  for line in beastFile:
    if line[0] == '#' :
      continue
    if line[:5] == 'state' or line[:6] == "Sample" :
      cols = line.strip().split()

      for t in traces:
        if t in cols :
          iTraces.append(cols.index(t))
        else :
          i = [k for k,c in enumerate(cols) if re.search(t, c) ]
          if len(i) :
            iTraces.extend(i)
            if report:
              print >> sys.stderr, "adding", " ".join([cols[k] for k in i])
          else :
            raise t + " not found."
      break

  if not len(iTraces) :
    raise "columns not found",",".join(traces)
  
  for line in beastFile:  
    l = line.strip().split()
    values.append( [float(l[x]) for x in iTraces] )

  return [ (cols[iTraces[k]],[v[k] for v in values]) for k in range(len(iTraces)) ]

from treeutils import convertDemographics

def setDemographics(trees) :
  """ Convert demographics for all trees """

  hasAll = []
  for tree in trees :
    missing = convertDemographics(tree)
    hasAll.append(missing == 0)

  return hasAll

