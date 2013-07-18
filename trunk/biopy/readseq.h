// This file is part of biopy.
// Copyright (C) 2013 Joseph Heled
// Author: Joseph Heled <jheled@gmail.com>
// See the files gpl.txt and lgpl.txt for copying conditions.

#if !defined(READSEQ_H)
#define READSEQ_H

// This is so ugly but,
// Tell me what is the right way to share code between 2 python extensions?

typedef unsigned char byte;
static byte const anynuc = 4;
static byte const gap = 5;

static byte*
readSequence(PyObject* seq, uint& nseq, bool const strip = false)
{
  if( ! PySequence_Check(seq) ) {
    PyErr_SetString(PyExc_ValueError, "not a sequence") ;
    return 0;
  }

  nseq = PySequence_Size(seq);
  byte* s = new byte[nseq];

  int nsite = 0;
  
  if( PyString_Check(seq) ) {
    const char* const c = PyString_AsString(seq);
    for(uint i = 0; i < nseq; ++i) {
      switch( c[i] ) {
	case 'a': case 'A': {
	  s[nsite] = 0;
	  break;
	}
	case 'g': case 'G': {
	  s[nsite] = 1;
	  break;
	}
	case 'c': case 'C': {
	  s[nsite] = 2;
	  break;
	}
	case 't': case 'T': {
	  s[nsite] = 3;
	  break;
	}
	case '-' : {
	  if( ! strip ) {
	    s[nsite] = gap;
	  } else {
	    nsite -= 1;
	  }
	  break;
	}
	default : {
	  s[nsite] = anynuc;
	  break;
	}
      }
      nsite += 1;
    }
  } else {
    for(uint i = 0; i < nseq; ++i) {
      PyObject* const o = PySequence_Fast_GET_ITEM(seq, i);
      
      s[nsite] = PyInt_AsLong(o);
      
      if( !( 0 <= s[nsite] && s[nsite] <= gap ) ) {
	delete [] s;
	PyErr_SetString(PyExc_ValueError, "Invalid sequence") ;
	return 0;
      }
      
      if( !strip || s[nsite] != gap ) {
	nsite += 1;
      }
    }
  }
  if( strip ) {
    nseq = nsite;
  }
  return s;
}

#endif
