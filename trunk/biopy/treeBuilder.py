## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#
# $Id:$ 

from Bio.Nexus import Trees, Nodes

""" Small helper in building BioPython trees. """

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
