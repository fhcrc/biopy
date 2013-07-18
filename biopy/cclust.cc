// This file is part of biopy.
// Copyright (C) 2010 Joseph Heled
// Author: Joseph Heled <jheled@gmail.com>
// See the files gpl.txt and lgpl.txt for copying conditions.


#include <Python.h>
// keeps asserts
#undef NDEBUG
#include <cassert>

#include <algorithm>
#include <unordered_map>
using std::unordered_map;

#include <vector>
using std::vector;


/**
     cans = [0]*len(seqs)
     for i in range(len(seq)-11) :
       p = seq[i:i+11]
       ms = matches.get(p)
       if ms :
         for m in ms:
           cans[m] += 1
     cans[iSeq] = 0
   
 **/

PyObject*
getCounts(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"seq", "matches", "result",
				 static_cast<const char*>(0)};
  char* seq = 0;
  int lseq = 0;
  PyObject* matches = 0;
  PyObject* result = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "s#OO", const_cast<char**>(kwlist),
				    &seq, &lseq, &matches, &result)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! PyDict_Check(matches) || ! PyList_Check(result) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args type");
    return 0;
  }

  //  Validity checks at a minimum to gain speed
  
  int const size = PyList_Size(result);
  
  long* cans = new long[size];
  std::fill(cans, cans+size, 0L);

  uint const nFragment = 11;
  
  for(char* s = seq; s < seq+lseq-nFragment; ++s) {
    // dirty, write EOS in string (temporarily)
    char const c = s[nFragment];
    s[nFragment] = 0;
    
    if( PyObject* const ms = PyDict_GetItemString(matches, s) ) {
      // in large sizes there is a penalty for using more abstract types
      // use explicit code for the list or tuple, which we actually use.
      
      if( PyList_Check(ms) ) {
	int const sz = PyList_GET_SIZE(ms);
	for(int j = 0; j < sz; ++j) {
	  PyObject* const o = PyList_GET_ITEM(ms, j);
	  int const m = PyInt_AS_LONG(o);                        assert ( 0 <= m && m < size );
	  cans[m] += 1;
	}
      } else if( PyTuple_Check(ms) ) {
	int const sz = PyTuple_GET_SIZE(ms);
	for(int j = 0; j < sz; ++j) {
	  PyObject* const o = PyTuple_GET_ITEM(ms, j);
	  int const m = PyInt_AS_LONG(o);                        assert ( 0 <= m && m < size );	  
	  cans[m] += 1;
	}
      } else { 
	PyObject* const msf = PySequence_Fast(ms, "error");
	int const sz = PySequence_Fast_GET_SIZE(msf);
	for(int j = 0; j < sz; ++j) {
	  PyObject* const o = PySequence_Fast_GET_ITEM(msf, j);
	  int const m = PyInt_AS_LONG(o);                        assert ( 0 <= m && m < size );
	  cans[m] += 1;
	}
	Py_DECREF(msf);
      }
    }
    s[nFragment] = c;
  }

  for(int k = 0; k < size; ++ k) {
    PyList_SET_ITEM(result, k, PyInt_FromLong(cans[k]));
  }

  delete [] cans;
  
  Py_INCREF(Py_None);
  return Py_None;
}

#include "readseq.h"

/**
def _buildLookup(seqs, asSet=False) :
  bl = dict()
  
  for j,s in enumerate(seqs) :
    for i in range(len(s)-11) :
      p = s[i:i+11]
      m = bl.get(p)
      if m is not None :
        m.append(j)
      else :
        bl[p] = [j]
  matches = dict([(x,set(y) if asSet else y) for x,y in bl.iteritems() if len(y) > 1])
        
  return matches
**/

typedef unsigned long ulong;

PyObject*
buildLookup(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"seqs", "fragmentSize", "removeSingles",
				 static_cast<const char*>(0)};
  PyObject* pSeqs = 0;
  int fragmentSize = 11;
  PyObject* pRemoveSingles = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "O|iO", const_cast<char**>(kwlist),
				    &pSeqs, &fragmentSize, &pRemoveSingles)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PySequence_Check(pSeqs) || fragmentSize <= 0 || fragmentSize > 21 ) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args");
    return 0;
  }

  bool const removeSingles = (pRemoveSingles == 0 || PyObject_IsTrue(pRemoveSingles));

  PyObject* d = PyDict_New();

  // code duplication unless I figure how to be clever with a template
  if( fragmentSize > 13 ) {
    typedef unsigned long long keytype;

    typedef unordered_map<keytype, vector<int> > mm;
    mm matches;

    typedef mm::value_type valtype;
    //typedef matches::value_type valtype;
  
    uint const nSeqs = PySequence_Size(pSeqs);
    PyObject* const sqs = PySequence_Fast(pSeqs, "error");

    vector<int> const emptyv;
    keytype const mask = ~(keytype(07) << 3*(fragmentSize-1));
  
    for(uint ns = 0; ns < nSeqs; ++ns) {
      PyObject* ps = PySequence_Fast_GET_ITEM(sqs, ns);
      uint lseq = 0;
      byte* s = readSequence(ps, lseq, true);
      if( ! s ) {
	PyErr_SetString(PyExc_ValueError, "wrong sequences");
	Py_DECREF(d);
	return 0;
      }
      
      keytype key = s[0];
      for(int l = 1; l < fragmentSize-1; ++l) {
	key <<= 3;
	key += s[l];
      }
      for(uint i = fragmentSize; i < lseq; ++i) {
	// i - i+fragmentSize
	key = (key << 3) + s[i-1];

	{
	  keytype k = s[i-fragmentSize];
	  for(uint l = i-fragmentSize+1; l < i; ++l) {
	    k = (k<<3) + s[l];
	  }
	  assert( k == key );
	}
      
	auto const& m = matches.find(key);
	if( m == matches.end() ) {
	  auto const& a = matches.insert( valtype(key, emptyv) );
	  a.first->second.push_back(ns);
	} else {
	  m->second.push_back(ns);
	}
	key &= mask;
      }
      delete [] s;
    }
    Py_DECREF(sqs);

    char key[fragmentSize+1];
    key[fragmentSize] = 0;

    for(auto i = matches.begin(); i != matches.end(); ++i) {
      vector<int>& v = i->second;
      if( removeSingles && v.size() == 1 ) {
	continue;
      }

      PyObject* ls = PyList_New(v.size());
      if( ! ls ) {
	return PyErr_NoMemory();
      }
      
      for(uint k = 0; k < v.size(); ++k) {
	PyList_SET_ITEM(ls, k, PyInt_FromLong(v[k]));
      }
      v.clear();
      {
	ulong k = i->first;
	for(int l = fragmentSize-1; l >= 0; --l) {
	  key[l] = "AGCTN"[k & 0x7];
	  k = k>>3;
	}
      }
      PyDict_SetItemString(d, key, ls);
      Py_DECREF(ls);
    }
  } else {
    assert(fragmentSize <= 13);

    typedef unordered_map<ulong, vector<int> > mm;
    mm matches;

    typedef mm::value_type valtype;
  
    uint const nSeqs = PySequence_Size(pSeqs);
    PyObject* const sqs = PySequence_Fast(pSeqs, "error");

    ulong ml = 1L;
    for(int l = 0; l < fragmentSize-1; ++l) {
      ml *= 5;
    }

    vector<int> const emptyv;
  
    for(uint ns = 0; ns < nSeqs; ++ns) {
      PyObject* ps = PySequence_Fast_GET_ITEM(sqs, ns);
      uint lseq = 0;
      byte* s = readSequence(ps, lseq, true);
      if( ! s ) {
	PyErr_SetString(PyExc_ValueError, "wrong sequences");
	Py_DECREF(d);
	return 0;
      }
      ulong key = s[0];
      for(int l = 1; l < fragmentSize-1; ++l) {
	key *= 5;
	key += s[l];
      }
      for(uint i = fragmentSize; i < lseq; ++i) {
	// i-fragmentSize to i
	key = key*5 + s[i-1];
      
	auto const& m = matches.find(key);
	if( m == matches.end() ) {
	  auto const& a = matches.insert( valtype(key, emptyv) );
	  a.first->second.push_back(ns);
	} else {
	  m->second.push_back(ns);
	}
	key -= s[i-fragmentSize] * ml;
      }
      delete [] s;
    }
    Py_DECREF(sqs);

    char key[fragmentSize+1];
    key[fragmentSize] = 0;

    for(auto i = matches.begin(); i != matches.end(); ++i) {
      vector<int>& v = i->second;
      if( removeSingles && v.size() == 1 ) {
	continue;
      }

      PyObject* const ls = PyList_New(v.size());
      if( ! ls ) {
	return PyErr_NoMemory();
      }
      
      for(uint k = 0; k < v.size(); ++k) {
	PyList_SET_ITEM(ls, k, PyInt_FromLong(v[k]));
      }
      v.clear();
      {
	ulong k = i->first;
	for(int l = fragmentSize-1; l >= 0; --l) {
	  ulong const r = k / 5;
	  key[l] = "AGCTN"[k - r*5];
	  k = r;
	}
      }
      PyDict_SetItemString(d, key, ls);
      Py_DECREF(ls);  // gets me every time. new ref passed again are
		      // not borrowed and should be released
      matches.erase(i);
    }
  }
  
  return d;
}


static PyMethodDef clustMethods[] = {
  {"counts",		(PyCFunction)getCounts, METH_VARARGS|METH_KEYWORDS,
   ""},
  {"lookupTable",	(PyCFunction)buildLookup, METH_VARARGS|METH_KEYWORDS,
   ""},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initcclust(void)
{
  //PyObject* const m =
  Py_InitModule("cclust", clustMethods);
}


	  
