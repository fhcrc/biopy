## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

"""
Generic helpers
===============
"""

__all__ = ["fileFromName", "sumViaLog"]

import os.path, gzip, bz2

def fileFromName(fname) :
  """A Python file object from (possibly compressed) disk file.

  If file has a common suffix (.gz,.bz2) use that as a guide. If fname does
  not exist, look for a compressed file with the same stem.
  
  :param fname: file name 
  """
  
  if os.path.exists(fname) :
    if fname.endswith(".gz") :
      return gzip.open(fname)
    if fname.endswith(".bz2") :
      return bz2.BZ2File(fname)
    return file(fname)
  if os.path.exists(fname + ".gz") :
    return gzip.open(fname + ".gz")
  if os.path.exists(fname + ".bz2") :
    return bz2.BZ2File(fname + ".bz2")
  
  raise IOError("no such file " + fname)

from numpy import log1p
import math

def sumViaLog(x) :
  """calculate log of the sum of numbers in x, each number in x represented by
  its log"""

  x = sorted(x)
  v = x[0]
  for z in x[1:] :
    # log-add v and z into v
    v = z + log1p(math.exp(v-z))
  return v
