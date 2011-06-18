#!/usr/bin/python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
from __future__ import division

import optparse, sys, os.path
from scipy.optimize import fmin_powell

from biopy.sequencesStats import ASD

parser = optparse.OptionParser(sys.argv[0] +
                               """ [OPTIONS] starbeast-file.xml mutation-rate[,rate2,,...]

  Heuristically estimate the birth rate and effective population size from a
  multispecies data in a *BEAST XML file.

  The alignments are grouped by ploidy (as specified in the *BEAST file), and
  one mutation rate is required per group. (the groups are ordered by increasing
  ploidy size), Mutation rates which are left empty will be estimated. For
  example 'starbeast.xml 0.1,' would use a fixed rate of 0.1 for the first
  group, and would estimate the rate for the second group.

  The estimate is based on Average Sequences Distance (ASD) - see
  http://arxiv.org/abs/1104.0727 for more details.
  """)

parser.add_option("-a", "--average-only", dest="meanOnly",
                  help="""Use only ASD. The default""" +
                  """is to use both mean and variance.""",
                  action="store_true", default = False)

options, args = parser.parse_args()

if len(args) != 2 :
  parser.print_help(sys.stderr)
  sys.exit(1)

meanOnly = bool(options.meanOnly)
xmlFileName = args[0]
mus = [(float(x) if x else None) for x in args[1].split(',')]

from biopy.beastXMLhelper import readBeastFile

try :
  d = {'DNA':None, 'species':None}

  readBeastFile(xmlFileName, d)

except Exception,e:
  # report error
  print >> sys.stderr, "Error:", e.message
  sys.exit(1)


def allPairs(lst) :
  if len(lst) < 2:
    return
  
  x0 = lst[0]
  for x in lst[1:] :
    yield (x0, x)
  for p in allPairs(lst[1:]) :
    yield p

def seqSim(seq1, seq2) :
  ne = sum([x!=y and x not in '-?' and y not in '-?'
            and x.upper() in "AGCT"
            and y.upper() in "AGCT"
            for x,y in zip(seq1, seq2)])
  return ne/len(seq1)

alignments = d['DNA']
species = d['species']['species']
genes = d['species']['genes'];
pl = sorted(dict.fromkeys([genes[x]['ploidy'] for x in genes]).keys())

if len(pl) != len(mus) :
  print >> sys.stderr, "**Error: please specify %d mutation rates." % len(pl)
  sys.exit(1)
  
between = dict([ (pr,[]) for pr in allPairs(species.keys())])
within = dict([ (sp,[]) for  sp in species.keys()])

for alk in alignments:
  al = alignments[alk]
  for spk in within:
    w = []
    sp = species[spk]
    for s1,s2 in allPairs(sp) :
      if s1 in al and s2 in al:
        w.append(seqSim(al[s1], al[s2]))
    if len(w) :
      within[spk].append([alk, w])

  for spi,spj in between:
    b = []
    iIndiv,jIndiv = species[spi],species[spj]
    for x1 in iIndiv:
      if x1 in al:
        for y1 in jIndiv:
          if y1 in al:
             b.append(seqSim(al[x1], al[y1]))
    
    if len(b) :
      between[(spi,spj)].append([alk, b])

ws = dict()
bs = dict()
for p in pl:
  aw = [[(sum([z**2 for z in x[1]]), sum(x[1]),len(x[1]))
         for x in within[s] if genes[x[0]]['ploidy'] == p]
        for s in within]
  ws[p] = (sum([sum([z[0] for z in x]) for x in aw]),
           sum([sum([z[1] for z in x]) for x in aw]) ,
           sum([sum([z[2] for z in x]) for x in aw]))

  ab = [[(sum([z**2 for z in x[1]]), sum(x[1]),len(x[1]))
         for x in between[s] if genes[x[0]]['ploidy'] == p]
        for s in between]
  bs[p] = (sum([sum([z[0] for z in x]) for x in ab]),
           sum([sum([z[1] for z in x]) for x in ab]), \
           sum([sum([z[2] for z in x]) for x in ab]) )

a1 = lambda a,b,c : b/c
v1 = lambda a,b,c : a/c - (b/c)**2

def combineMus(fmus, emus) :
  mus = []
  i = 0
  for k,z in enumerate(fmus) :
    if z is None :
      mus.append(emus[i])
      i += 1
    else :
      mus.append(z)
  return mus

def solve(fmus, avgs, ns) :
  def dis(x) :
    x = [abs(z) for z in x]
    ne, lam = x[:2]
    if len(x) > 2 :
      mus = combineMus(fmus, x[2:])
    else :
      mus = fmus
    s = 0.0
    for (p,(ws,bs)),mu in zip(avgs,mus) :
      w,b = ASD(ns, p*ne, mu, lam, True)
      s +=  (ws[-1]+bs[-1]) * ( (w[0] - a1(*ws))**2 + (b[0] - a1(*bs))**2 )
      if not meanOnly:
        s +=  (ws[-1]+bs[-1]) * ( (w[1]**0.5 - v1(*ws)**0.5)**2 + \
                                  (b[1]**0.5 - v1(*bs)**0.5)**2 )
    return s
  
  extra = sum([x is None for x in fmus])
  mu0 = [x for x in fmus if x is not None][-1]
  return fmin_powell(dis, [1,1] + [mu0]*extra, disp=0, full_output=1)

fsl = solve(mus, [[x , [ws[x],  bs[x]]] for x in ws], len(species))
sl = [abs(x) for x in fsl[0]]
print "Ne= %g Lambda= %g" % tuple(sl[:2])

if fsl[5] != 0 :
  print "***WARNING: failed to converge"

def toPer(x,y) :
  return "%.3g%%" % (100 * (1 - x/y))

amus = combineMus(mus, sl[2:])

for (p,(w,b)),mu in zip([[x , [(a1(*z),v1(*z)) for z in [ws[x],  bs[x]]]]
                         for x in ws], amus):
  xw,xb = ASD(len(species), p*sl[0], mu, sl[1], True)
  print "ploidy=%.3g mu= %.5g\tASD(w)= %g ASD(b)= %g error= %s,%s" \
        % (p, mu, w[0], b[0], toPer(xw[0],w[0]), toPer(xb[0],b[0]))
  if not meanOnly:
    sw, sb, sxw, sxb = [x**0.5 for x in (w[1], b[1], xw[1], xb[1])]
    print "\t\t\tSdSD(w)= %g SdSD(b)= %g error= %s,%s" % (sw, sb, toPer(sxw,sw), toPer(sxb,sb))
