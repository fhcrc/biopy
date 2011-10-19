#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

import optparse, sys, os.path
parser = optparse.OptionParser(usage = """ %prog [OPTIONS] taxa

Draw a gene tree using the coalescent with a given population size function.

Taxa: either a number, a comma separated list of labels, or a {template,n}
pair. The last specifies 'n' taxa given by applying 0,1,...,n-1 to template.""")

parser.add_option("-n", "--ntrees", dest="ntrees",
                  help="""Number of trees to generate """
                  + """(default %default)""", default = "1", metavar="N")

parser.add_option("-p", "--population", dest="demospec", metavar="SPEC",
                  help= """The population function. A comma separated list: the
                  first element is the population size at time zero (now), the
                  rest come in pairs""",
                  default = "1")

parser.add_option("-b", "--basename", dest="base", metavar="NAME",
                  help= """Tree name prefix""", default = None)

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
  taxa = ["x%d" % k for k in range(int(taxTxt))]
else:
  if "," not in taxTxt:
    errExit("Invalid taxa spesification.")
    
  t = taxTxt.split(',')
  if len(t) == 2 and t[1].isdigit() and "%" in t[0] :
    taxa = [t[0] % k for k in range(int(t[1]))]
  else :
    taxa = t

from biopy import __version__
from biopy.treeutils import toNewick, TreeLogger
from biopy.coalescent import sampleCoalescentTree
from biopy.demographic import Demographic

demog = Demographic.parseDemographic(options.demospec)

tlog = TreeLogger(options.nexfile, argv = sys.argv, version = __version__)

for nt in range(nTrees) :
  t = sampleCoalescentTree(demog, taxa)
  tlog.outTree(toNewick(t, attributes="attributes"),
               name = options.base + ("_%d" % nt) if options.base else None)

tlog.close()
