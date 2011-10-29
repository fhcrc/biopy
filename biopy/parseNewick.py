## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

from __future__ import division
from treeutils import TreeBuilder

__all__ = ["parseNewick"]

from cchelp import parsetree

# Reference implementation in python - supposed to be the same as the C version
# in cchelp.

def _getStuff(s, sep) :
  e = 0
  while s[e] != sep or s[e-1] == '\\' :
    e += 1
  return e

def _findIndex(s, ch, stopAt) :
  e = 0

  while s[e] != ch:
    if s[e] in stopAt:
      return e
    e += 1 
    if e >= len(s) :
      return -1
  return e
  
def _parseAttributes(s) :
  vals = []
  eat = 0
  while s[0] != ']' :
    if s[0] == ',' :
      s = s[1:]
      eat += 1

    nameEnd = _findIndex(s, '=', ",]\"{}")
    if s[nameEnd] != '=' :
      raise "error"
    name = s[:nameEnd]
    s = s[nameEnd+1:]
    eat += nameEnd+1
    
    if s[0] == '"' :
        e = _getStuff(s[1:],'"')
        v = s[1:e+1]
        s = s[e+2:]
        eat += e+2
    elif s[0] == '{' :
        e = _getStuff(s[1:], '}')
        v = s[1:e+1]
        s = s[e+2:]
        eat += e+2
    else :
        e = _findIndex(s, ',', "]")
        if e == -1 :
          raise "error"
        v = s[:e]
        s = s[e:]
        eat += e

    vals.append((name.strip(),v.strip()))
  return eat,vals


def _skipSpaces(txt) :
  i = 0
  while i < len(txt) and txt[i].isspace() :
    i += 1
  return i

def _readSubTree(txt, nodesList) :
  n = _skipSpaces(txt)
  txt = txt[n:]
  if txt[0] == '(' :
    subs = []
    while True:
      n1 = _readSubTree(txt[1:], nodesList)
      n += 1 + n1
      txt = txt[1+n1:]
      subs.append(len(nodesList)-1)

      n1 = _skipSpaces(txt)
      n += n1
      txt = txt[n1:]
      if txt[0] == ',' :
        continue
      if txt[0] == ')' :
        nodesList.append([None, None, subs, None])
        n += 1
        txt = txt[1:]
        break
      raise "error"
  else :
    # a terminal
    n1 = 0
    while not txt[n1].isspace() and txt[n1] not in ":[,()]":
      n1 += 1
    nodesList.append([txt[:n1], None, None, None])
    n += n1
    txt = txt[n1:]

  n1 = _skipSpaces(txt)
  txt = txt[n1:]
  n += n1
  
  if len(txt) and txt[0] == '[' and txt[1] == '&':
      n1, vals = _parseAttributes(txt[2:])
      n1 += 3
      n1 += _skipSpaces(txt[n1:])
      n += n1
      txt = txt[n1:]
      nodesList[-1][3] = vals
       
  if len(txt) and txt[0] == ':' :
      n += 1
      n1 = _skipSpaces(txt[1:])
      n += n1
      txt = txt[1+n1:]
      n1 = 0
      while n1 < len(txt) and txt[n1] in ".0123456789+-Ee" :
        n1 += 1
      b = float(txt[:n1])
      txt = txt[n1:]
      n += n1
      nodesList[-1][1] = b
      
  return n



def _build(nodes, weight=1.0, rooted=True, name='') :
  tb = TreeBuilder(weight=weight, rooted = rooted, name=name)
  t = [None]*len(nodes)
  for k, x in enumerate(nodes):
    if x[0] is not None :
      t[k] = tb.createLeaf(x[0])
    else :
      t[k] = tb.mergeNodes([ [t[l], nodes[l][1]] for l in x[2]])
    if x[3] is not None:
      t[k].data.attributes = dict(x[3])
  return tb.finalize(t[-1])

def parseNewick(txt, weight=1.0, rooted=True, name='') :
#  if 0 :
#    nodes = []
#    _readSubTree(txt, nodes)
#  else :
  nodes = parsetree(txt)
  return _build(nodes, weight=weight, rooted = rooted, name=name)
