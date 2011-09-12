## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$
""" Handling of BEAST log and tree files."""

from __future__ import division
import sys, re

__all__ =  ["readTraces", "setDemographics"]

def readTraces(beastFile, traces) :
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
    if line[:5] == 'state' :
      cols = line.strip().split()

      for t in traces:
        if t in cols :
          iTraces.append(cols.index(t))
        else :
          i = [k for k,c in enumerate(cols) if re.search(t, c) ]
          if len(i) :
            iTraces.extend(i)
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

from demographic import LinearPiecewisePopulation as LPP

def _tod(xt, yt) :
  xt = [float(x) for x in xt.split(',') if len(x)]
  yt = [float(x) for x in yt.split(',')]
  return LPP(yt, xt)

def _toDemog(dtxt) :
  return _tod(*dtxt.split('|'))

def _toDemog1(dmt, dmv) :
  return _tod(dmt if dmt is not None else "", dmv)

def setDemographics(trees) :
  """ Convert demographic function stored in BEAST trees attributes to biopy
  demographic."""
  
  dmf = "dmf"
  dmv,dmt = "dmv", "dmt"

  hasAll = []
  # Read demographics in
  for tree in trees :
    has = True
    for i in tree.all_ids() :
      data = tree.node(i).data
      if hasattr(data, "attributes") :
        if dmf in data.attributes:
          dtxt = data.attributes[dmf]
          data.demographic = _toDemog(dtxt)
        elif dmv in data.attributes:
          data.demographic = _toDemog1(data.attributes.get(dmt), data.attributes.get(dmv))
        else :
          has = False
    hasAll.append(has)

  return hasAll
