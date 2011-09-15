#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

import optparse, sys, os.path
parser = optparse.OptionParser(usage = """ %prog [OPTIONS] species-tree

Generate gene trees under a multispecies with migration model. 'species-tree' is
either a nexus file or a tree in Newick format. Population sizes are specified
with the tree. Migration rates are either specified in the tree or generated
using a stochastic model.
""")

parser.add_option("-n", "--ntrees", dest="ngenetrees", metavar="N",
                  help="""Number of gene trees per species tree """
                  + """(default %default)""", default = "1") 

parser.add_option("-m", "--migration", dest="migration",
                  help="""migration model """
                  + """(default %default)""", default = "1") 

parser.add_option("-e", "--redraw", dest="redraw", metavar="E",
                  help="""Change model based parameters every E"""
                  + """ (default %default) trees.""", default = "1") 

parser.add_option("-t", "--per-species", dest="ntips", metavar="N",
                  help="""Number of individuals per species. Ignored if tree
                  contains labels (default %default).""", default = "2") 

parser.add_option("-o", "--nexus", dest="nexfile", metavar="FILE",
                  help="Print trees in nexus format to FILE", default = None)

parser.add_option("-l", "--log", dest="sfile", metavar="FILE",
                  help="Log species trees in nexus format to FILE", default = None)

#parser.add_option("", "--total", dest="total", metavar="N",
#                  help="""Stop after processing N species trees.""",
#                  default = None) 

# options for M and S. more scenarios?

options, args = parser.parse_args()

nGeneTrees = int(options.ngenetrees) ; assert nGeneTrees > 0
nTips = int(options.ntips)           ; assert nTips > 0

# nTotal = int(options.total) if options.total is not None else -1           
reDraw = int(options.redraw)


nexusTreesFileName = args[0]

from biopy.treeutils import toNewick, TreeLogger, setLabels, \
     convertDemographics, revertDemographics
from biopy import INexus, beastLogHelper, __version__
from biopy.geneTreeWithMigration import GeneTreeSimulator, setIMrates

tlog = TreeLogger(options.nexfile, argv = sys.argv, version = __version__)
if options.sfile is not None :
  slog =  TreeLogger(options.sfile, argv = sys.argv, version = __version__)
else :
  slog = None
  
if os.path.isfile(nexusTreesFileName) :
  i = INexus.INexus().read(fileFromName(nexusTreesFileName))
  trees = [tree for tree in i]
else :
  trees = [INexus.Tree(nexusTreesFileName)]
  
hasPops = beastLogHelper.setDemographics(trees)
if not all(hasPops) :
  print >> sys.stderr, "Error: Missing demographic(s)"
  sys.exit(1)

tipNameTemplate = "%s_tip%d"

for tree in trees:
  terms = tree.get_terminals()
  if not setLabels([tree]) :
    for i in terms:
      data = tree.node(i).data
      if not hasattr(data, "labels") :
        labels = [tipNameTemplate % (data.taxon, k) for k in range(nTips)]

  misPimr = convertDemographics(tree, formatUnited = None,
                                formatSeparated = ("imrv", "imrt"),
                                dattr = "pimr")
  hasPimr = misPimr == 1 and not hasattr(tree.node(tree.root).data, "pimr")
  if hasPimr:
    for i in tree.all_ids() :
      node = tree.node(i)
      if node.succ :
        node.data.ima = [tree.node(ch).data.pimr for ch in node.succ]
        
  counter = 0
  for k in range(nGeneTrees) :
    if counter == 0 :
      if not hasPimr:
        setIMrates(tree) 
      s = GeneTreeSimulator(tree)
    gt,nim = s.simulateGeneTree()
    tlog.outTree(gt, {'N' : str(nim)})
    if slog is not None :
      # prepare demographic attributes
      revertDemographics(tree)
      # prepare ima attributes - move to children nodes
      for i in tree.all_ids() :
        node = tree.node(i)
        if node.succ :
          for m,ch in zip(node.data.ima, node.succ):
            tree.node(ch).data.pimr = m

      revertDemographics(tree, dattr="pimr", formatSeparated = ("imrv", "imrt"))
      slog.outTree(toNewick(tree, attributes="attributes"))
                   
    counter = (counter + 1) % reDraw

  #nTotal -= 1
  #if nTotal == 0 :
  #  break

tlog.close()

if slog is not None :
  slog.close()
  
