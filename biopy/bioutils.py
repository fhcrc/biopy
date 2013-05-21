## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

__all__ = ["readFasta"]

def readFasta(fasFile) :
  assert (hasattr(fasFile, 'read') and hasattr(fasFile, 'write'))
  head = None
  for l in fasFile :
    if l[0] == '#' :
      continue
    if l[0] == '>':
      if head :
        yield (head,body)
      head = l.strip()
      body = ""
    else :
      body = body + l.strip()
      
  if head :
    yield (head,body)
