#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

import optparse, sys, os.path
parser = optparse.OptionParser(usage = """ %prog [OPTIONS] taxa

Draw a species tree with population sizes using birth/death and simple
strategies for population sizes.

Taxa: either a number, a comma separated list of labels, or a {template,n}
pair. The last specifies 'n' taxa given by applying 0,1,...,n-1 to template.""")

parser.add_option("-n", "--ntrees", dest="ntrees",
                  help="""Number of trees to generate """
                  + """(default %default)""", default = "1", metavar="N")

parser.add_option("-b", "--birth", dest="birth",
                  help="""Birth rate (default %default)""", default = 1.0)

parser.add_option("-d", "--death", dest="death",
                  help="""Death rate (default %default)""", default = 0.0)

parser.add_option("-p", "--population-distribution", dest="popdist", metavar="DIST",
                  help=\
"""Draw starting populations from this distribution (default constant %default).

A number by itself specifies a constant. Otherwise, a comma separated list where the
first character specifies the distribution.

Supportd distributions:
 (1) u,l,h : uniform between l and h (Example 'u,2,4').
 (2) e,r   : exponential with rate r (mean 1/r) (Example 'e,2').
 (3) l,m,s : log-normal with mean m (in real space) and std s (in log space).
 (4) g,s,l : gamma with shape s and scale l.
 (5) i,a,b : inverse-gamma with shape(alpha) a and scale(beta) b.
 """,
                  default = "1")

parser.add_option("-c", "--continuous", dest="continuous", metavar="DIST",
                  help=\
"""Population size functions are linear over branch and continuous at divergence
points (population splits into two parts). The argument specifies a distribution
for the rate of change in population.  The population at the start of the branch
(i.e at divergence) is P_0 2^{-r b}, where P_0 is the population at the end of
the branch, r is the rate and b is the branch length, as a fraction of total
tree height.""",
                  default = None)

parser.add_option("-o", "--nexus", dest="nexfile", metavar="FILE",
                  help="Print trees in nexus format to FILE", default = None)

options, args = parser.parse_args()

def errExit(msg = None) :
  if msg is not None :
    print >> sys.stderr, msg
  parser.print_help(sys.stderr)
  sys.exit(1)

nTrees = int(options.ntrees)                    ; assert nTrees > 0

if len(args) < 1 :
  errExit("No taxa?")

taxTxt = args[0]
if taxTxt.isdigit() :
  taxa = ["s%d" % k for k in range(int(taxTxt))]
else:
  if "," not in taxTxt:
    errExit("Invalid taxa spesification.")
    
  t = taxTxt.split(',')
  if len(t) == 2 and t[1].isdigit() and "%" in t[0] :
    taxa = [t[0] % k for k in range(int(t[1]))]
  else :
    taxa = t

from collections import namedtuple
from biopy import __version__
from biopy.treeutils import toNewick, treeHeight, TreeLogger
from biopy.birthDeath import drawBDTree
from biopy.randomDistributions import parseDistribution

birthRate = float(options.birth)            ; assert birthRate > 0.0
deathRate = float(options.death)            ; assert 0 <= deathRate <= birthRate

popStartDist = parseDistribution(options.popdist)

tlog = TreeLogger(options.nexfile, argv = sys.argv, version = __version__)

SetPops = namedtuple("SetPops", "pop data")

def setPops(tree, nid, height, popStartDist, rateDist) :
  node = tree.node(nid)

  if not node.succ :
    node.data.pop0 = popStartDist.sample()
  else :
    p = [setPops(tree, x, height, popStartDist, rateDist) for x in node.succ]

    popEnd = [x.pop * 2**(-rateDist.sample() * (x.data.branchlength/height)) for x in p]
    node.data.pop0 = sum(popEnd)
    for x,v in zip(p, popEnd) :
      x.data.pope = v
      
  return SetPops(node.data.pop0, node.data)


for nt in range(nTrees) :
  t = drawBDTree(birthRate, deathRate, len(taxa))

  for n,tx in zip(t.get_terminals(),taxa) :
    t.node(n).data.taxon = tx

  for n in t.all_ids() :
    t.node(n).data.attributes = dict()

  if options.continuous is None :
    for n in t.all_ids() :
      t.node(n).data.attributes['dmv'] = popStartDist.sample()
  else :
    rate = parseDistribution(options.continuous)
    setPops(t, t.root, treeHeight(t), popStartDist, rate)
    for n in t.all_ids() :
      d = t.node(n).data
      if n != t.root :
        d.attributes['dmv'] = "{%f,%f}" % (d.pop0,d.pope)
        d.attributes['dmt'] = "%f" % d.branchlength
      else :
        d.attributes['dmv'] = "%f" % d.pop0
    
  tlog.outTree(toNewick(t, attributes="attributes") )

tlog.close()
