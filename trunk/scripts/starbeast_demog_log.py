#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

from __future__ import division

from biopy import INexus, beastLogHelper
from biopy.treeutils import getTreeClades, toNewick, countNexusTrees

import optparse, sys, os.path
parser = optparse.OptionParser(os.path.basename(sys.argv[0]) +
                               """ [OPTIONS] beast-species-tree output-log-file

  Generate a tracer-readable file with population size information per clade, to
  truly assess the ESS's on population sizes. """)

parser.add_option("-t", "--threshold", dest="threshold",
                  help="""Exclude all clades with posterior probability lower"""
                  + """than threshold (percent, default 5)""", default = "5") 

parser.add_option("-b", "--burnin", dest="burnin",
                  help="Burn-in amount (percent, default 10)", default = "10")

parser.add_option("-r", "--report", dest="report",
                  help="""Output a sorted list of clade percentages to """
                  + """standard output.""",
                  action="store_true", default = False)

options, args = parser.parse_args()

nexusTreesFileName = args[0]
burnIn =  float(options.burnin)/100.0
tops = True

threshold = float(options.threshold)/100.0

if len(args) > 1 :
  if os.path.exists(args[1]) :
    print >> sys.stderr, "not overwriting", args[1]
    sys.exit(1)
    
  logFile = file(args[1],'w')
else :
  logFile = sys.stdout
  
cladesDict = dict()
treeTopologies = list()

tree = INexus.INexus().read(nexusTreesFileName).next()
for n in tree.get_terminals():
  data = tree.node(n).data
  cladesDict[frozenset([data.taxon])] = []

nTrees = countNexusTrees(nexusTreesFileName)

statesNos = []
nexusFile = INexus.INexus()
for stateNo,tree in enumerate(nexusFile.read(nexusTreesFileName,
                                             slice(int(burnIn*nTrees), -1, 1))):
  beastLogHelper.setDemographics([tree])
  if tops:
    treeTopologies.append(toNewick(tree, topologyOnly=True))
  
  clades = getTreeClades(tree)
  for c,node in clades:
    cladeSet = frozenset(c)
    l = cladesDict.get(cladeSet)
    if l is None :
      l = cladesDict[cladeSet] = []
    demo = node.data.demographic
    b = node.data.branchlength if node.id != tree.root else demo.naturalLimit()
    if b is None :
      # constant root
      b = 0
    l.append( (stateNo, (demo.population(0), demo.population(b))) )
  for n in tree.get_terminals():
    data = tree.node(n).data
    demo = data.demographic
    l = cladesDict[frozenset([data.taxon])]
    l.append( (stateNo, (demo.population(0),
                         demo.population(node.data.branchlength))) )
  statesNos.append(tree.name[6:])

stateNo += 1
if bool(options.report) :
  cs1 = sorted([(len(cladesDict[c])/stateNo,c) for c in cladesDict], reverse=True)
  for p,c in cs1:
    print "%s%.2f \t%s" % ("*" if p < threshold else " ",p*100, '-'.join(sorted(c)))
  
torem = []
if threshold < 1.0 :
  for c in cladesDict:
    if len(cladesDict[c])/stateNo  < threshold:
      torem.append(c)
for c in torem:
  del cladesDict[c]
      
print >> logFile, "state",
if tops:
  print >> logFile, "\ttopology",
  
vals = []
for c in sorted(cladesDict, key = lambda c : (len(c), '-'.join(sorted(c)))):
  vals.append(cladesDict[c])
  k = '-'.join(sorted(c))
  print >> logFile, '\t' + k + "_b" + '\t' + k + "_e",
print >> logFile

def dist1(c1,c2) :
  d = 0
  for x in c1:
    if x not in c2:
      d += 1
  for x in c2:
    if x not in c1:
      d += 1
  return d

# Assigning unique numbers to topologies is a hard problem. I use a very rough
# heuristic here: first sort topologies by posterior probability (higher to
# lower), then order them by starting with the first and picking as next the
# first one with the smallest distance (Robinson-Foulds) to the last one.

if tops:
  uniqTops = dict()
  for nt,top in enumerate(treeTopologies):
    if top in uniqTops:
      uniqTops[top].append(nt)
    else:
       uniqTops[top] = [nt]

  sTops = sorted([(len(l),top) for top,l in uniqTops.items()], reverse=True)
  trees = [INexus.Tree(x[1]) for x in sTops]
  cTrees = [dict.fromkeys([tuple(sorted(c)) for c,n in getTreeClades(t)])
            for t in trees]
  order = []
  k = 0

  while True:
    s = sTops.pop(k)
    tree,ctree = trees.pop(k),cTrees.pop(k)
    order.append((s,tree,ctree))

    if len(cTrees) == 0 :
      break
    d1 = [dist1(ctree, t) for t in cTrees]
    k = d1.index(min(d1))

  topsIndex = dict([(x[0][1],k) for k,x in enumerate(order)])
  
cur = [list(x[0][1]) for x in vals]
sno = 0

while sno < stateNo:
  print >> logFile, statesNos[sno],
  if tops:
    print >> logFile, '\t', topsIndex[treeTopologies[sno]],
    
  for x,y in cur:
    print >> logFile, '\t',x,'\t',y,
  print >> logFile

  for cv,vl in zip(cur,vals):
    if len(vl) and vl[0][0] == sno :
      cv[0] = vl[0][1][0]
      cv[1] = vl[0][1][1]
      vl.pop(0)
  sno += 1

logFile.close()


