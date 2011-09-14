#! /usr/bin/env python
## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division
import random

from biopy.combinatorics import _prod as prod
from biopy import INexus, __version__
from biopy.birthDeath import yuleTimes
from biopy.treeutils import toNewick, TreeLogger

def NHforest(z) :
  """ Number of labeled history forrests for this forrest."""
  u = [(n,k) for n,k in z if n > 0]
  n,k = sum([x[0] for x in u]), prod([x[1] for x in u])
  return prod(range(2,n+1))//k

def wchoice(pr) :
  """ Pick one option based on probablities in pr (sum(pr) == 1.0 """
  z = random.random()
  w = 0.0
  for k,x in enumerate(pr):
    w += x
    if z < w:
      return k
  raise "error"
    
def prepareTree(t, nid) :
  node = t.node(nid)
  if node.data.taxon :
     lr = (0,1)
  else :
    l = [prepareTree(t, x) for x in node.succ]
    s = sum([x[0] for x in l])+1
    lr = (s, s * prod([x[1] for x in l]))

  node.data.lrsize = lr
  return lr


def getOrder(t) :
  # internal descendants of root
  fnodes = [t.node(x) for x in t.node(t.root).succ if t.node(x).succ]
  lrnodes = [x.data.lrsize for x in fnodes]
  order = []

  while len(fnodes) :
    if len(fnodes) == 1:
      inext = 0
    else:
      w = []
      for k,fn in enumerate(fnodes):
        z = [t.node(x).data.lrsize for x in fn.succ]
        w.append(NHforest(z + lrnodes[:k] + lrnodes[k+1:]))
      s = sum(w)
      w = [x/s for x in w]
      inext = wchoice(w)

    order.append(fnodes[inext])
    n = fnodes.pop(inext)
    lrnodes.pop(inext)
    for a in [t.node(x) for x in n.succ]:
      if a.succ:
        fnodes.append(a)
        lrnodes.append(a.data.lrsize)

  return order

def setBranches(tr, order, intervals) :
  r = len(order)+1
  tr.node(tr.root).data.h = r
  for k,x in enumerate(order):
    p = tr.node(x.prev)
    x.data.h = r-1-k
    x.data.branchlength = sum(intervals[x.data.h:p.data.h])

  for t in tr.get_terminals():
    x = tr.node(t)
    p = tr.node(x.prev)
    x.data.branchlength = sum(intervals[:p.data.h])
    

import optparse, sys, os.path

prog = os.path.basename(sys.argv[0])
parser = optparse.OptionParser(prog + """ [options] tree

Sample a tree uniformly from the space of all ranked trees having the same
unranked topology as 'tree'. Both coalescent (constant size) and pure-birth
(Yule) are supported.""")

parser.add_option("-n", "--ntrees", dest="ntrees",
                  help="Number of samples (trees) to generate (=1)", default = "1")

parser.add_option("", "--nexus", dest="nexfile",
                  help="Print trees in nexus format to file", default = None)

parser.add_option("-p", "--population", dest="popSize",
                  help="Population size for coialscent model (constant, =1)",
                  default = 1)

parser.add_option("-b", "--birth-rate", dest="birthRate",
                  help="Birth rate for model yule model", default = None)

options, args = parser.parse_args()

try :
  tree = INexus.Tree(args[0])
except Exception,e:
  print >> sys.stderr, "*** Error: failed in reading tree: ",\
        e.message,"[",args[0],"]"
  sys.exit(1)

if not all([len(tree.node(x).succ) in [0,2] for x in tree.all_ids()]) :
  print >> sys.stderr, "*** Error: tree is not binary or not fully resolved"
  sys.exit(1)

nTrees = int(options.ntrees)
if nTrees <= 0 :
  sys.exit(1)

if options.birthRate is not None:
  
  if float(options.popSize) != 1.0:
    print >> sys.stderr, "*** Error: choose exactly one between coalescent and yule"
    sys.exit(1)
    
  birthRate = float(options.birthRate)
  assert birthRate > 0
  def getIntervals(n) :
    return yuleTimes(birthRate, n)

else :
  pop = float(options.popSize)
  def getIntervals(n) :
    return [random.expovariate(((k*(k-1))//2)/pop) for k in range(n, 1, -1)]

tlog = TreeLogger(options.nexfile, argv = sys.argv, version = __version__)

prepareTree(tree, tree.root)

for k in range(nTrees) :
  order = getOrder(tree)
  intervals = getIntervals(len(order)+2)
  setBranches(tree, order, intervals)
  tlog.outTree(toNewick(tree))

tlog.close()
