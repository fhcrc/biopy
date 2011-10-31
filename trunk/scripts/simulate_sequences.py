#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

import optparse, sys, os.path
parser = optparse.OptionParser("""%prog [OPTIONS] trees-file-or-tree
""")

parser.add_option("-n", "--seqlen", dest="seqlen", metavar="N",
                  help="""Alignment length."""
                  + """(default %default)""", default = 400) 

parser.add_option("-m", "--model", dest="model",
                  help="""Substitution model."""
                  + """(default %default)""", default = "JC,1") 

parser.add_option("-a", "--annotate", dest="annotate",
                  help="""Output trees annotated with sequences. Default is to
                  generate an NEXUS alignment file per tree.""",
                  action="store_true", default = False)

parser.add_option("-o", "--nexus", dest="nexfile", metavar="FILE",
                  help="Export trees in nexus format to file", default = None)


from biopy import INexus, submodels, __version__
from biopy.treeutils import toNewick, TreeLogger

options, args = parser.parse_args()

if os.path.isfile(args[0]) :
  trees = INexus.INexus().read(file(args[0]))
else :
  try :
    trees = [INexus.Tree(args[0])]
  except Exception,e:
    print >> sys.stderr, "*** Error: failed to parse tree: ",\
          e.message,"[",args[0],"]"
    sys.exit(1)

seqLen = int(options.seqlen)

def parseFreqs(freqs) :
  pi = [0]*4
  qMat = [1]*7
  
  while len(freqs) >= 2 :
    code,f = freqs[:2]

    try:
      f = float(f)
    except :
      raise ValueError(f)
    
    freqs = freqs[2:]
    
    if len(code) == 1 :
      c = submodels.SubstitutionModel.NUCS.index(code)
      if not 0 <= c < 4 :
        raise RuntimeError("Can't parse nucleotide frequency" + ",".join(freqs))
      if not 0 < f < 1 :
        raise RuntimeError("Illegal nucleotide frequency: " + f)
      pi[c] = f
    elif len(code) == 2 :
      c1,c2 = sorted([submodels.SubstitutionModel.NUCS.index(c) for c in code])
      if c1 == c2 :
        raise RuntimeError("Illegal nucleotide specification: ", code)
      if f <= 0 :
        raise RuntimeError("Illegal rate: " + f)
      qMat[3*c1 + c2-1] = f
      
  notAssigned = sum([x == 0 for x in pi])
  if notAssigned > 0 :
    if sum(pi) >= 1 :
      raise("Illegal nucleotide stationary frequencies: sum to more than 1")
    p = (1 - sum(pi))/notAssigned
    for k in range(4) :
      if pi[k] == 0 :
        pi[k] = p
  else :
    if sum(pi) != 1 :
      print >> sys.stderr, "Nucleotide frequencies do not sum to 1 - normalizing..."
      pi = [x/sum(pi) for x in pi]

  return pi,[x/qMat[-1] for x in qMat][:-1]

try :
  model = options.model.split(',')
  modelName = model[0].upper()
  if modelName == 'JC' :
    mu = float(model[1]) if model[1] else 1
    smodel = submodels.JCSubstitutionModel(mu)
  elif modelName == 'HKY':
    mu = float(model[1]) if model[1] else 1
    kappa = float(model[2]) if model[2] else 1
    freqs = model[3:]
    pi = parseFreqs(freqs)[0]
    smodel = submodels.HKYSubstitutionModel(mu = mu, kappa = kappa, pi = pi)
  elif modelName == 'GTR':
    mu = float(model[1]) if model[1] else 1
    freqs = model[2:]
    pi,qMat = parseFreqs(freqs)
    smodel = submodels.StationaryGTR(mu = mu, m = qMat, pi = pi)
except RuntimeError as r:
  print >> sys.stderr, "error in parsing model." + r.message
  sys.exit(1)

tlog = None
if options.annotate :
  tlog = TreeLogger(options.nexfile, argv = sys.argv, version = __version__)
else :
  from biopy.INexus import exportMatrix
  
for count,tree in enumerate(trees):
  smodel.populateTreeSeqBySimulation(tree, seqLen)

  for n in tree.get_terminals() :
    data = tree.node(n).data
    if not hasattr(data, "attributes") :
      data.attributes = dict()
    data.attributes['seq'] = smodel.toNucCode(data.seq)

  if tlog :
    treeTxt = toNewick(tree, attributes = 'attributes')
    tlog.outTree(treeTxt, name = tree.name)
  else :
    d = dict()
    for n in tree.get_terminals() :
      data = tree.node(n).data
      d[data.taxon] = data.attributes['seq']

    fname = tree.name if tree.name else ("alignment%d" % count)
    f = file(fname + ".nex", "w")
    exportMatrix(f, d)
    
if tlog:
  tlog.close()
