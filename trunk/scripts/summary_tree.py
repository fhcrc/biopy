#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

from __future__ import division

import optparse, sys, os.path
from time import time

from biopy.genericutils import fileFromName
from biopy.treeutils import toNewick, countNexusTrees
from biopy.treesPosterior import minPosteriorDistanceTree
from biopy import INexus


parser = optparse.OptionParser(sys.argv[0] +
                               """ [OPTIONS] posterior-trees-nexus 

  Generate a single summary tree for a set of posterior trees. """)

parser.add_option("-b", "--burnin", dest="burnin",
                  help="Burn-in amount (percent, default 10)", default = "10")

parser.add_option("-e", "--every", dest="every",
                  help="""thin out - take one tree for every 'e'. Especially \
useful if you run out of memory (default all, i.e. 1)""", default = "1")

parser.add_option("-n", "--ntops", dest="ntops",
                  help="""Use the top 'n' topologies from the posterior (default\
 10)""", default = "10")

parser.add_option("-l", "--limit", dest="limit",
                  help="""run at most 'l' seconds, trying out topologies in""" + \
                  """ decreasing order of support and in random order for""" + \
                  """ topologies with equal support (default -1, i.e no""" + \
                  """ time limits). Note that I have not tested timing code under""" + \
                  """ M$windows or OSX""", 
                  default = "-1") 

parser.add_option("-p", "--progress", dest="progress",
                  help="Print out progress messages to terminal (standard error)",
                  action="store_true", default = False)

options, args = parser.parse_args()

nexusTreesFileName = args[0]
try :
  nexFile = fileFromName(nexusTreesFileName)
except Exception,e:
  # report error
  print >> sys.stderr, "Error:", e.message
  sys.exit(1)

progress = options.progress
burnIn =  float(options.burnin)/100.0
every = int(options.every)
ntops = int(options.ntops)
limit = int(options.limit)

if progress:
  print >> sys.stderr, "counting trees ...,",
nTrees = countNexusTrees(nexusTreesFileName)

# establish trees

nexusReader = INexus.INexus()

if progress:
  print >> sys.stderr, "reading %d trees ...," % int((nTrees * (1-burnIn) / every)),

trees = list(nexusReader.read(nexFile, slice(int(burnIn*nTrees), -1, every)))

if 1:
  if progress:
    print >> sys.stderr, "collect topologies and sort ...,",
  topology = dict()
  for tree in trees :
    k = toNewick(tree, None, topologyOnly=True)
    if k not in topology :
      topology[k] = [tree,]
    else :
      topology[k].append(tree)

  allt = topology.items()
  # Sort by amount of support
  allt.sort(reverse=1, key = lambda l : len(l[1]))

if limit <= 0 :
  if progress:
    pPost = sum([len(x[1]) for x in allt[:ntops]]) / len(trees)
    print >> sys.stderr, """using top %d topologies out of %d, covering %.1f%% of\
  posterior topologies...,""" % \
          (min(ntops, len(allt)), len(allt), 100*pPost),

  k = ntops-1
  lLast = len(allt[k][1])
  while k < len(allt) and lLast == len(allt[k][1]) :
    k += 1
  if k > ntops :
     print >> sys.stderr
     print >> sys.stderr, """*** WARNING ***:  %d additional topologies have the \
same support as the %dth one (%d trees, %.3f%%)""" % (k - ntops, ntops, lLast, lLast/len(trees))
       
  candidates = [x[1][0] for x in allt[:ntops]]
else :
  import random

  candidates = []
  while len(candidates) < len(allt) :
    l = []
    k = len(candidates)
    lFirst = len(allt[k][1])
    while k < len(allt) and lFirst == len(allt[k][1]) :
      l.append(allt[k][1][0])
      k += 1
    random.shuffle(l)
    # print len(candidates), len(l), lFirst
    candidates.extend(l)

  #for x in candidates:
  #  print toNewick(x, None, topologyOnly=1)
  print >> sys.stderr, """trying in order %d topologies (time permitting) ...""" \
        % len(candidates),
  
bestTree, bestScore = None, float('infinity')

if progress:
  print >> sys.stderr, "searching ...,",

nCandidatesTried = 0

if limit > 0 :
  startTime = time()
  
for tree in candidates:
  tr, score = minPosteriorDistanceTree(tree, trees)
  if score < bestScore:
    bestScore = score
    bestTree = tr
  nCandidatesTried += 1
  if limit > 0:
    if time() - startTime >= limit:
      print >> sys.stderr, "time limit reached ...,",
      break

if limit > 0:
  print >> sys.stderr, "examined %d topologies in %.1f seconds," \
        % (nCandidatesTried, time() - startTime),
  
if progress or limit > 0:
  print >>  sys.stderr, "done." 

print toNewick(bestTree)
