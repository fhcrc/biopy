## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

from __future__ import division

from biopy import INexus, beastLogHelper
from biopy.treeutils import getTreeClades

import optparse, sys, os.path
parser = optparse.OptionParser(sys.argv[0] +
                               """ [OPTIONS] beast-species-tree output-log-file

  Generate a tracer-readable file with population size information per clade, to
  really judge the ESS's on population sizes. """)

parser.add_option("-t", "--threshold", dest="threshold",
                  help="""Exclude all clades with posterior probability lower"""
                  + """than threshold""", default = "5") 

parser.add_option("-r", "--report", dest="report",
                  help="Report all clades percentages",
                  action="store_true", default = False)

options, args = parser.parse_args()

threshold = float(options.threshold)/100.0

if len(args) > 1 :
  if os.path.exists(args[1]) :
    print >> sys.stderr, "not overwriting", args[1]
    sys.exit(1)
    
  logFile = file(args[1],'w')
else :
  logFile = sys.stdout
  
cladesDict = dict()

tree = INexus.INexus().read(args[0]).next()
for n in tree.get_terminals():
  data = tree.node(n).data
  cladesDict[frozenset([data.taxon])] = []

statesNos = []
nexusFile = INexus.INexus()
for stateNo,tree in enumerate(nexusFile.read(args[0])):
  beastLogHelper.setDemographics([tree])
  
  clades = getTreeClades(tree)
  for c,node in clades:
    cladeSet = frozenset(c)
    l = cladesDict.get(cladeSet)
    if l is None :
      l = cladesDict[cladeSet] = []
    demo = node.data.demographic
    b = node.data.branchlength if node.id != tree.root else demo.naturalLimit()
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
vals = []
for c in sorted(cladesDict, key = lambda c : (len(c), '-'.join(sorted(c)))):
  vals.append(cladesDict[c])
  k = '-'.join(sorted(c))
  print >> logFile, '\t' + k + "_b" + '\t' + k + "_e",
print >> logFile

cur = [list(x[0][1]) for x in vals]
sno = 0

while sno < stateNo:
  print  >> logFile, statesNos[sno],
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


