#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

import optparse, sys, os.path
parser = optparse.OptionParser("""%prog [OPTIONS] species-trees-file

Compute probabilities of all possible gene trees evolving according to the
multispecies coalescent inside species tree.

Due to combinatorical explosion, this is good only for *very* small cases.

Population size and tip information is provided as attributes of species tree
nodes.

For example:

((A[&dmv=1,labels=a1 a2]:1,B[&dmv=1,labels=b]:1)[&dmv=1]:1,C[&dmv=1,labels=c]:2)[&dmv=1]

'dmv' gives the population size, 'labels' the individuals per species. 
""")

parser.add_option("-o", "--nexus", dest="nexfile",
                  help="Print trees in nexus format to file", default = None)

options, args = parser.parse_args()

nexusTreesFileName = args[0]

from biopy import INexus, speciesTreesGeneTrees, beastLogHelper, \
     treeCombinatorics, __version__
from biopy.treeutils import TreeLogger

try :
  tree = INexus.Tree(args[0])
except Exception,e:
  # report error
  print >> sys.stderr, "Error:", e.message
  sys.exit(1)

has = beastLogHelper.setDemographics([tree])
if not has[0] :
  print >> sys.stderr, "Error: Missing demographic(s)"
  sys.exit(1)

for i in tree.get_terminals() :
  has = False 
  data = tree.node(i).data
  if hasattr(data, "attributes") :
    if "labels" in data.attributes:
      l = data.attributes["labels"].split(' ')
      data.labels = l
      has = len(l) > 0

  if not has :
    print >> sys.stderr, \
          "Error: Missing or invalid labels for taxon '%s'." % data.taxon
    sys.exit(1)

tlog = TreeLogger(options.nexfile, argv = sys.argv, version = __version__)
  
compat = None
trees = speciesTreesGeneTrees.compatibleGeneTreesInSpeciesTree(tree, compat)

for count,t in enumerate(sorted(trees, key = lambda x : x[0], reverse=True)) :
  tlog.outTree(treeCombinatorics.toNewick(t[1]), {'R' : None,
                                                  'W' : "%0.14f" % pt[0]})

tlog.close()
