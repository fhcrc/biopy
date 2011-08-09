#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

import optparse, sys, os.path

import pylab

from biopy.genericutils import fileFromName

from biopy import INexus, beastLogHelper, treePlotting
from biopy.treeutils import getTreeClades, countNexusTrees

parser = optparse.OptionParser(os.path.basename(sys.argv[0]) +
                               """ [OPTIONS] posterior-trees.nexus figfilename
        """)

parser.add_option("-b", "--burnin", dest="burnin",
                  help="Burn-in amount (percent, default 10)", default = "10")

parser.add_option("-e", "--every", dest="every",
                  help="""thin out - take one tree for every 'e'. Especially \
useful if you run out of memory (default all, i.e. 1)""", default = "1")

parser.add_option("-m", "--mean-tree", dest="mainTree", action="append",
                  help="""Show tree on top of river (the tree from \
starbeast_posterior_popsizes for instance)""", default = None)

parser.add_option("-s", "--single", dest="singleTopology",
                  help="Restrict to a single topology (either the most frequent \
topology, or the one from the mean-tree when given""",
                  action="store_true", default = False)

parser.add_option("-c", "--color", dest="color",
                  help=""" """,
                  action="store_true", default = False)

parser.add_option("", "--yclip", dest="yclip",
                  help="""Clip display of root demographic using upper HPD \
(0 < level < 1)""", default = 1) 

parser.add_option("", "--positioning", dest="positioning",
                  type='choice', choices=['mean', 'taxonmean', 'between'],
                  help="Method of positioning internal nodes based on clade \
taxa (mean [default] | taxonmean | between).", default = 'mean') 

parser.add_option("-p", "--progress", dest="progress",
                  help="Print out progress messages to terminal (standard error)",
                  action="store_true", default = False)

options, args = parser.parse_args()

progress = options.progress
every = int(options.every)
burnIn =  float(options.burnin)/100.0
colorTops = options.color

if options.positioning == 'mean' :
  positioning = treePlotting.descendantMean
elif options.positioning == 'taxonmean':
   positioning = treePlotting.descendantWeightedMean
elif options.positioning == 'between':
   positioning = treePlotting.descendantBetween
else:
  raise

if len(args) != 2 :
  parser.print_help(sys.stderr)
  sys.exit(1)

mainTree = []
if options.mainTree :
  mainTree = [INexus.Tree(t) for t in options.mainTree]
  beastLogHelper.setDemographics(mainTree)
  
nexusTreesFileName = args[0]

try :
  nexFile = fileFromName(nexusTreesFileName)
except Exception,e:
  # report error
  print >> sys.stderr, "Error:", e.message
  sys.exit(1)

if progress:
  print >> sys.stderr, "counting trees ...,",
  
nTrees = countNexusTrees(nexusTreesFileName)

nexusReader = INexus.INexus()

if progress:
  print >> sys.stderr, "reading %d trees ...," % int((nTrees * (1-burnIn) / every)),

trees = []
for tree in nexusReader.read(nexFile, slice(int(burnIn*nTrees), -1, every)):
  beastLogHelper.setDemographics([tree])
  trees.append(tree)
  
if progress:
  print >> sys.stderr, "preparing ...,",

wd = treePlotting.getSpacing(trees, .3)
ref = mainTree[0] if mainTree else None
taxOrder = treePlotting.getTaxaOrder(trees, refTree = ref,
                                     reportTopologies = colorTops) 
oo,refTree = taxOrder[:2]
xs = dict([(t,k*wd) for k,t in enumerate(oo)])

cd = None
if options.singleTopology:
  d1 = dict([(n.id,(c,n)) for c,n in getTreeClades(refTree, True)])
  cd = dict([(frozenset(c), (frozenset(d1[n.succ[0]][0]), frozenset(d1[n.succ[1]][0]))
              if n.succ else None)  for c,n in d1.values()])

if progress:
  print >> sys.stderr, "plotting ...,",

pylab.ioff()
fig = pylab.figure()

heights = []
alpha = min(10./len(trees), .3)

if colorTops :
  import colorsys
  tops = taxOrder[2].items()
  tops.sort(key = lambda x : len(x[1]))
  mainHue = 100./256
  c0 = (mainHue + len(tops[-1][1])/len(trees)/2) % 1
  for top in tops:
    topPercent = len(top[1])/len(trees)
    c2 = (c0 + topPercent/2) % 1
    col = colorsys.hsv_to_rgb(c2,1.0,184/256.)
    #print top[0],len(top[1]), c2, [x*256 for x in col]
    for k in top[1]:
      tree = trees[k]
      for x in tree.get_terminals() :
        tree.node(x).data.x = xs[tree.node(x).data.taxon]
      h = treePlotting.drawTree(tree, tree.root, cd, positioning = positioning,
                                color=col, alpha = alpha)
      heights.append(h)
    c0 = (c0 + topPercent) % 1
    
    
else :
  dtrees = trees

  for tree in dtrees:
    for x in tree.get_terminals() :
       tree.node(x).data.x = xs[tree.node(x).data.taxon]
    h = treePlotting.drawTree(tree, tree.root, cd, positioning = positioning,
                              color="lime", alpha = alpha)
    heights.append(h)

from biopy.bayesianStats import hpd

yc = float(options.yclip)
yclip = 0 < yc < 1
if yclip :
  hmin = hpd(heights, yc)[1]
else :
  hmin = 0
  
for tree,col in zip(mainTree, ["red", "blue", "green"]):
  for x in tree.get_terminals() :
     tree.node(x).data.x = xs[tree.node(x).data.taxon]
  h = treePlotting.drawTree(tree, tree.root, None,
                            positioning = positioning, fill=False, color=col)
  hmin = max(h, hmin)

for ni in tree.get_terminals():
  node = tree.node(ni)
  t = pylab.text(xs[node.data.taxon], -0.1, node.data.taxon, fontsize = 12,
                 va='top', ha='center')
  
ymax = hmin if yclip else  pylab.ylim()[1]

pylab.ylim((-0.5*ymax/10, ymax))

pylab.savefig(args[1], dpi=300)

if progress:
  print >>  sys.stderr, "done." 
