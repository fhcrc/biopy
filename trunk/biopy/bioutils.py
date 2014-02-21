## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

__all__ = ["readFasta"]

def readFasta(fasFile, stripGaps = False, comments = False) :
  assert (hasattr(fasFile, 'read') and hasattr(fasFile, 'write'))
  head = None
  if comments:
    cc = []
  for l in fasFile :
    if l[0] == '#' or l[0] == ';' :
      if comments :
        cc.append(l.strip())
      continue
    if l[0] == '>':
      if head :
        if stripGaps:
          body = ''.join([x for x in body if x != '-'])
        if body :
          # skip over empty sequences
          if comments :
            yield (head,body,comment)
          else :
            yield (head,body)

      head = l.strip()
      body = ""
      if comments :
        comment = cc
        cc = []
    else :
      body = body + l.strip()
      
  if head :
    if stripGaps:
      body = ''.join([x for x in body if x != '-'])
    if comments :
      yield (head,body,comment)
    else :
      yield (head,body)
