## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division

"""
   Tree representation based on the ideas of Mau et all.
   
   Bayesian Phylogenetic Inference via Markov Chain Monte Carlo Methods
   `{PDF<http://citeseer.ist.psu.edu/rd/27056960%2C5592%2C1%2C0.25%2CDownload/
   http://citeseer.ist.psu.edu/cache/papers/cs/5768/
   ftp:zSzzSzftp.stat.wisc.eduzSzpubzSznewtonzSzlastfinal.pdf/mau98bayesian.pdf>`_)

   Provide transformation to and from BioPython Nexus tree and strings in newick
   format. 
"""

__all__ = ["mauCanonical", "mau2Tree"]

import random
from numpy import array
from treeutils import TreeBuilder, nodeHeights


def mauCanonical(tree, internalNodes = True) :
  """ Convert tree to Mau representation.

  Result is a tuple containing two lists [taxa-names, internal-nodes-heights].
  if 'internalNodes' the tuple contains a third list containing the node id of
  each internal node and a boolean indication if the childrens have been swaped
  of not. 

  >>> mauCanonical(Trees.Tree('((chimp:1,human:1):1,gorilla:2)')) in \
            [[['gorilla', 'chimp', 'human'], [2.0, 1.0]], \
            [['chimp', 'human', 'gorilla'], [1.0, 2.0]], \
            [['human', 'chimp', 'gorilla'], [1.0, 2.0]]]
  True

  @param tree:
  @type tree:  ITrees.Tree
  @param internalNodes: when True, return low level information about the
  location of internal nodes and if their siblings has been swapped.
"""

  nh = nodeHeights(tree)
  return mauCanonicalSub(tree, tree.node(tree.root), nh, internalNodes)

def _mauCanonicalSub(t, node, nodeHeights, internalNodes = False) :
  """ Return representation of subtree of 't' rooted at node."""
  
  if node.data.taxon :
    r = [[node.data.taxon,], []]
    if internalNodes :
      r.append([])
    return r
  
  h = nodeHeights[node.id]

  s = list(node.succ)
  swap = bool(random.randint(0,1))
  if swap:
    s = [s[1],s[0]]
  lf,rt = [_mauCanonicalSub(t, t.node(x), nodeHeights, internalNodes) for x in s]
  r = [lf[0] + rt[0], lf[1] + [h,] + rt[1]]
  if internalNodes :
    r.append(lf[2] + [(node.id,swap),] + rt[2])
  return r


def _mau2treeInternal(leaves, heights, internals) :
  """ Reconstruct a tree from a mau representation.

  Arguments are the elements returned from mauCanonical, where 'heights' is
  assured to be an array.
  """
  
  if len(heights) == 0 :
    return ([leaves[0], None],0, None)

  mp = heights.argmax()
  mpp1 = mp+1
  
  lf, lh, li = _mau2treeInternal(leaves[0:mpp1], heights[0:mp], internals[0:mp])
  rt, rh, ri = _mau2treeInternal(leaves[mpp1:], heights[mpp1:], internals[mpp1:])
  h,(i,swap) = heights[mp], internals[mp]
  t = [[lf, h -lh, li], [rt, h - rh, ri]]
  if swap :
    t = t[1],t[0]
  return (t, h, i)

def mau2Tree(mauTree) :
  """ Convert a tree in Mau representation to a BioPython Nexus Tree object."""
  t = TreeBuilder()
  rep, h, nid = _mau2treeInternal(mauTree[0], array(mauTree[1]), mauTree[2])
  r = _mau2NexusInternal((rep, 0), t)

  return t.finalize(r[0])

def _mau2NexusInternal(rep, treeb) :
  rep, bl = rep[:2]

  if not rep[1]:
    # a leaf
    n = treeb.createLeaf(rep[0])
  else :
    l1 = _mau2NexusInternal(rep[0], treeb)
    l2 = _mau2NexusInternal(rep[1], treeb)
    n = treeb.mergeNodes([l1,l2])
    
  return n,bl


if __name__ == '__main__':
  import doctest
  doctest.testmod()
