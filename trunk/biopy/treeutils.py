## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

""" Small helper in building BioPython trees. """

import operator
from genericutils import fileFromName

# Bio.Nexus.Tree stuff

from Bio.Nexus import Trees, Nodes

__all__ = ["TreeBuilder", "getClade", "getTreeClades", "getCommonAncesstor",
           "countNexusTrees", "toNewick", "nodeHeights", "nodeHeight", "treeHeight"]


class TreeBuilder(object) :
  """ A basic helper for building BioPython trees.

  Use:
   - tb = TreeBuilder()
   - Create some terminals via leaf = tb.createLeaf(name)
   - use mergeNodes to successively merge nodes.
   - Call tb.finalize(root-node) to get the tree
  """
  
  def __init__(self) :
    self.t = Trees.Tree()

  def createLeaf(self, name) :
    nd = Trees.NodeData()
    nd.taxon = name
    leaf = Nodes.Node(nd)
    self.t.add(leaf, None)
    return leaf

  def newNode(self) :
    nd = Trees.NodeData()
    node = Nodes.Node(nd)
    self.t.add(node, None)
    return node
  
  def mergeNodes(self, n1, h1, n2, h2) :
    nd = Trees.NodeData()
    node = Nodes.Node(nd)
    self.t.add(node, None)

    n1.set_prev(node.id)
    n1.data.branchlength = h1
    n2.set_prev(node.id)
    n2.data.branchlength = h2

    node.add_succ([n1.id, n2.id])
    return node

  def finalize(self, n) :
    t = self.t
    t.node(t.root).set_succ(n.succ)
    for p in n.succ :
      t.node(p).set_prev(t.root)
    t.kill(n.id)
    return t

def getClade(tree, nodeId) :
  n = tree.node(nodeId)
  if n.data.taxon :
    return [nodeId,]
  return reduce(operator.add, [getClade(tree, x) for x in n.succ])

def getCommonAncesstor(tree, taxaIds) :
  return reduce(tree.common_ancestor, taxaIds)

def _getTreeClades_i(tree, nodeID) :
  node = tree.node(nodeID)
  if node.data.taxon :
    return ([node.data.taxon], [])

  cl = [_getTreeClades_i(tree, ch) for ch in node.succ]
  allt = apply(operator.add, [x[0] for x in cl])
  clades = [(allt,node)]
  for cl1 in [x[1] for x in cl if x[1]] :
    clades.extend(cl1)
  return (allt, clades)

def getTreeClades(tree, trivialClades = False):
  """ Clades of subtree as a list of (taxa-list, tree-node).

  taxa-list is a list of strings.
  """

  c = _getTreeClades_i(tree, tree.root)[1]

  if trivialClades :
    for n in tree.get_terminals():
      nd = tree.node(n)
      c.append(([nd.data.taxon], nd))
    
  return c


def countNexusTrees(nexFileName) :
  """ Number of trees in a nexus file."""
  nexFile = fileFromName(nexFileName)
  c = 0
  for l in nexFile:
    if l.startswith("tree ") :
      c += 1
  return c

def toNewick(tree, nodeId = None, topologyOnly = False, attributes = None) :
  """ BioPython tree or sub-tree to unique newick format.

  Children nodes are sorted (via text), so representation is unique and does not
  depend on arbitrary children ordering.
  """
  
  if not nodeId :
    nodeId = tree.root
  node = tree.node(nodeId)
  data = node.data
  if data.taxon :
    rep = data.taxon
  else :
    reps = [toNewick(tree, n, topologyOnly, attributes) for n in node.succ]
    reps.sort()
    rep = "(" + reps[0] + "," + reps[1] + ")"

  if attributes is not None :
    attrs = getattr(data, attributes, None)
    if attrs is not None and len(attrs) :
      s = '[&'
      for a in attrs:
        s += (a+'='+str(attrs[a])+",")
      s = s[:-1] + ']'
      rep += s
    
  if not topologyOnly and nodeId != tree.root :
    rep = rep + (":%r" % data.branchlength)
  return rep


def _getNodeHeight(tree, n, heights, w = 0):
  if n.id in heights:
    return heights[n.id]
  
  if len(n.succ) == 0 :
    heights[n.id] = 0.0
    return 0.0

  i = n.succ[w]
  if i not in heights:
    if n.succ[1-w] in heights:
      w = 1-w
      i = n.succ[w]
      
  c = tree.node(i)
  h = _getNodeHeight(tree, c, heights, 1-w) + c.data.branchlength
  
  heights[n.id] = h
  return h

def nodeHeights(tree, nodes) :
  """ Return a mapping from node ids to node heights"""

  w = 0
  
  heights = dict()
  for n in nodes :
    _getNodeHeight(tree, n, heights, w)
    w = 1-w
    
  return heights

def nodeHeight(tree, node) :
  if node.data.taxon :
    h = 0
  else :
    ch = [tree.node(i) for i in node.succ]

    h = max([n.data.branchlength + nodeHeight(tree, n) for n in ch])
    
  return h

def treeHeight(tree) :
  return nodeHeights(tree, [tree.node(tree.root)])[tree.root]
