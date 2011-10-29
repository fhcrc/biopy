## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

""" Generic helpers """

import os.path, gzip, bz2

def fileFromName(fname) :
  """A Python file object from (possibly compressed) disk file.

  If file has a common suffix (.gz,.bz2) use that. If 'fname' does not exist,
  look for a compressed file with the same stem.
  
  @param fname: File name 
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
