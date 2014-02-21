// This file is part of biopy.
// Copyright (C) 2010 Joseph Heled
// Author: Joseph Heled <jheled@gmail.com>
// See the files gpl.txt and lgpl.txt for copying conditions.


#include <Python.h>
// keeps asserts
#undef NDEBUG
#include <cassert>

#include <algorithm>
#include <memory>
#include <unordered_map>
using std::unordered_map;

#include <vector>
using std::vector;

struct Element {
  long	value;
  uint	index;

  static bool compare(Element const& lhs, Element const& rhs)
  {
    return (rhs.value != lhs.value ? rhs.value < lhs.value : rhs.index < lhs.index);
  }

  PyObject* asTuple(void) const {
    PyObject* t = PyTuple_New(2);
    PyTuple_SET_ITEM(t, 0, PyInt_FromLong(value));
    PyTuple_SET_ITEM(t, 1, PyInt_FromLong(index));
    return t;
  }
};

class SimpleQ {
public:
  SimpleQ(long* vals, uint _size) :
    size(0),
    q(new Element[_size])
    {
      for(uint k = 0; k < _size; ++k) {
	long v = vals[k];
	if( v != 0 ) {
	  q[size].value = v;
	  q[size].index = k;
	  size += 1;
	}
      }
      std::make_heap(q,q+size,Element::compare);
    }

  ~SimpleQ() {
    delete [] q;
  }

  Element* pop(void) {
    if( size == 0 ) {
      return 0;
    }
    
    std::pop_heap(q, q+size,Element::compare);
    size -= 1;
    return q + size;
  }
  
private:
  uint		size;
  Element*	q;
};

static const char* const qecapName = "QELEMENTS";

static void
elementsqDestructor(PyObject* const o)
{
  if( PyCapsule_IsValid(o, qecapName) ) {
    void* c = PyCapsule_GetPointer(o, qecapName);
    delete reinterpret_cast<SimpleQ*>(c);
  }
}

PyObject*
popElement(PyObject* const o)
{
  if( PyCapsule_IsValid(o, qecapName) ) {
    void* c = PyCapsule_GetPointer(o, qecapName);
    SimpleQ& q = *reinterpret_cast<SimpleQ*>(c);
    if( Element* e = q.pop() ) {
      return e->asTuple();
    }
    Py_INCREF(Py_None);
    return Py_None;
  }
  PyErr_SetString(PyExc_ValueError, "wrong args.") ;
  return 0;
}

PyObject*
popQelement(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"queue",
				 static_cast<const char*>(0)};
  PyObject* q = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "O", const_cast<char**>(kwlist),
				    &q) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  return popElement(q);
}

typedef unordered_map<ulong, vector<int> > mm;
static const char* const capName = "MATCHTABLE";

static inline int
ntoi(char c) {
  switch( c ) {
    case 'a': case 'A': {
      return 0;
    }
    case 'g': case 'G': {
      return 1;
    }
    case 'c': case 'C': {
      return 2;
    }
    case 't': case 'T': {
      return 3;
    }
  }
  return 4;
}

static inline const vector<int>*
findFragment(mm const& matches, char* s, uint const fragmentSize)
{
  ulong key = ntoi(s[0]);
  for(uint l = 1; l < fragmentSize; ++l) {
    key *= 5;
    key += ntoi(s[l]);
  }
  auto const& m = matches.find(key);
  if( m == matches.end() ) {
    return 0;
  } else {
    return &m->second;
  }
}


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
  static const char *kwlist[] = {"seq", "matches", "result", "seqslens", "fragmentSize",
				 static_cast<const char*>(0)};
  char* seq = 0;
  int lseq = 0;
  PyObject* matches = 0;
  PyObject* result = 0;
  PyObject* lseqs = 0;
  int fragmentSize = 11;

  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "s#OO|Oi", const_cast<char**>(kwlist),
				    &seq, &lseq, &matches, &result, &lseqs, &fragmentSize)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  bool const isNative = PyCapsule_IsValid(matches, capName);
  // if true, must have lseqs
  bool const returnQueue = PyInt_Check(result);
  
  if( ! ( (isNative || PyDict_Check(matches)) && (returnQueue || PyList_Check(result))
	  &&  ( !returnQueue || lseqs) ) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args type");
    return 0;
  }

  //  Validity checks at a minimum to gain speed
  
  int const size = returnQueue ? PyInt_AS_LONG(result) : PyList_Size(result);

  if( lseqs ) {
    if( ! PyList_Check(lseqs) || PyList_Size(lseqs) != size ) {
      PyErr_SetString(PyExc_ValueError, "incompatible args");
      return 0;
    }
  }

  const mm* const pmatches = isNative ?
    reinterpret_cast<mm*>(PyCapsule_GetPointer(matches, capName)) : 0;
  
  long* const cans = new long[size];
  std::fill(cans, cans+size, 0L);

  uint const nFragment = fragmentSize;

  for(char* s = seq; s < seq+lseq-nFragment; ++s) {
    // dirty, write EOS in string (temporarily)
    char const c = s[nFragment];
    s[nFragment] = 0;
    
    if( isNative ) {
      if( const vector<int>* i = findFragment(*pmatches, s, nFragment) ) {
	for(auto m = i->begin(); m != i->end(); ++m) {
	  cans[*m] += 1;
	}
      }
    } else {
      if( PyObject* const ms = PyDict_GetItemString(matches, s) ) {
	// in large sizes there is a penalty for using more abstract types
	// use explicit code for a list or tuple (instead of sequence),
	// which is the typical use case.
      
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
    }
    s[nFragment] = c;
  }

  if( lseqs ) {
    PyObject* const lsqs = PySequence_Fast(lseqs, "error");
    long const m = std::numeric_limits<long>::min() + 1;
    for(int k = 0; k < size; ++ k) {
      long n = cans[k];
      if( n > 0 ) {
	int const l = PyInt_AS_LONG(PySequence_Fast_GET_ITEM(lsqs, k));
	n = long(m * (double(n)/ ((l < lseq) ? l : lseq)));
      }
      cans[k] = n;
    }
    Py_DECREF(lsqs);
  }

  PyObject* r;
  
  if( returnQueue ) {
    SimpleQ* qe = new SimpleQ(cans, size);
    r = PyCapsule_New(qe, qecapName, elementsqDestructor);
  } else {
    for(int k = 0; k < size; ++ k) {
      PyList_SET_ITEM(result, k, PyInt_FromLong(cans[k]));
    }
    Py_INCREF(Py_None);
    r = Py_None;
  }

  delete [] cans;
  return r;
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


static void
matchesTableDestructor(PyObject* const o)
{
  if( PyCapsule_IsValid(o, capName) ) {
    void* c = PyCapsule_GetPointer(o, capName);
    delete reinterpret_cast<mm*>(c);
  }
}



PyObject*
buildLookup(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"seqs", "fragmentSize", "removeSingles", "native",
				 static_cast<const char*>(0)};
  PyObject* pSeqs = 0;
  int fragmentSize = 11;
  PyObject* pRemoveSingles = 0;
  PyObject* pNative = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "O|iOO", const_cast<char**>(kwlist),
				    &pSeqs, &fragmentSize, &pRemoveSingles, &pNative)) {
    PyErr_SetString(PyExc_ValueError, "wrong args (1).") ;
    return 0;
  }

  if( ! (PySequence_Check(pSeqs) && fragmentSize > 0 && fragmentSize <= 21) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args (2)");
    return 0;
  }

  bool const removeSingles = (pRemoveSingles == 0 || PyObject_IsTrue(pRemoveSingles));
  bool const native = (pNative != 0 && PyObject_IsTrue(pNative));

  PyObject* d = !native ? PyDict_New() : 0;

  // code duplication unless I figure how to be clever with a template
  if( fragmentSize > 13 ) {
    assert( ! native ) ; // (FIXME)
    
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

    std::unique_ptr<mm> pmatches(new mm);
    mm& matches = *pmatches.get();

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
	if( !native ) {
	  Py_DECREF(d);
	}
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

    if( native ) {
      if( removeSingles ) {
	for(auto i = matches.begin(); i != matches.end(); /** **/) {
	  if( i->second.size() == 1 ) {
	    i = matches.erase(i);
	  } else {
	    ++i;
	  }
	}
      }
      d = PyCapsule_New(pmatches.release(), capName, matchesTableDestructor);
    } else {
      char key[fragmentSize+1];
      key[fragmentSize] = 0;

      for(auto i = matches.begin(); i != matches.end(); i = matches.erase(i) /** ++i **/) {
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
	// matches.erase(i);
      }
    }
  }
  
  return d;
}


static PyMethodDef clustMethods[] = {
  {"popq",		(PyCFunction)popQelement, METH_VARARGS|METH_KEYWORDS,
   "Pop top element from the (encapsulated) priority queue returned from 'counts'."},
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


// PyMODINIT_FUNC
// initcclustnew(void)
// {
//   //PyObject* const m =
//   Py_InitModule("cclustnew", clustMethods);
// }
	  
