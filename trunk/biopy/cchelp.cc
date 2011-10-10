// This file is part of biopy.
// Copyright (C) 2010 Joseph Heled
// Author: Joseph Heled <jheled@gmail.com>
// See the files gpl.txt and lgpl.txt for copying conditions.

#include <Python.h>
#include "numpy/arrayobject.h"
#include <string>
#include <vector>
using std::vector;

static inline bool
is1Darray(PyArrayObject* a) {
  return PyArray_Check(a) && PyArray_NDIM(a) == 1;
}

static PyObject*
nonEmptyIntersection(PyObject*, PyObject* args)
{
  PyArrayObject* al;
  PyArrayObject* ar;
  PyArrayObject* sSet;

  if( !PyArg_ParseTuple(args, "OOO", &al, &ar, &sSet) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  
  if( !(is1Darray(al) && is1Darray(ar) && is1Darray(sSet) ) ) {
    PyErr_SetString(PyExc_ValueError, "not arrays.") ;
    return 0;
  }

  int const aLen = PyArray_DIM(al, 0);
  if( ! (aLen == PyArray_DIM(ar, 0) && aLen == PyArray_DIM(sSet, 0)) ) {
    PyErr_SetString(PyExc_ValueError, "length mismatch.") ;
    return 0;
  }

  const long* const alb = reinterpret_cast<long*>(PyArray_BYTES(al));
  const long* const arb = reinterpret_cast<long*>(PyArray_BYTES(ar));
  const long* const aSetb = reinterpret_cast<long*>(PyArray_BYTES(sSet));
  
  bool inone = false;

  //std::cout << aLen << std::endl;
  //  << " " <<  alb[k] << " " <<  aSetb[k] << std::endl;
  for(int k = 0; k < aLen; ++k) {

    if( alb[k] && aSetb[k] ) {
      inone = true;
      break;
    }
  }
  
  PyObject* result =  Py_False;
  if( inone ) {
    for(int k = 0; k < aLen; ++k) {
      if( arb[k] && aSetb[k] ) {
	result = Py_True;
	break;
      }
    }
  }

  Py_INCREF(result);
  return result;
}

static inline double
darrayElement(PyObject* a, int k)
{
  return PyFloat_AsDouble(PySequence_Fast_GET_ITEM(a, k));
}
  
static PyObject*
demoLPpopulation(PyObject*, PyObject* args)
{
  PyObject* xvals;
  PyObject* vals;
  double t;

  if( !PyArg_ParseTuple(args, "OOd", &vals, &xvals, &t) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( !(PySequence_Check(vals) && PySequence_Check(xvals) ) ) {
    PyErr_SetString(PyExc_ValueError, "not sequences.") ;
    return 0;
  }

  unsigned int const ll =  PySequence_Size(xvals);

  unsigned int k = 0;
  while ( k < ll && darrayElement(xvals, k) < t ) {
    k += 1;
  }
  double pop;
  if( k == ll ) {
    pop = darrayElement(vals, k);
  } else {
    // make t dt relative to start
    double width = darrayElement(xvals, k);
    if( k > 0 ) {
      double const x = darrayElement(xvals, k-1);
      t -= x;
      width -= x;
    }
    double vk = darrayElement(vals, k);
    double vkp1 =  darrayElement(vals, k+1);;
    pop = vk + (t/width) * (vkp1 - vk);
  }
  
  return PyFloat_FromDouble(pop);
}

static PyObject*
demoLPintegrate(PyObject*, PyObject* args)
{
  PyObject* xvals;
  PyObject* vals;
  double xHigh;

  if( !PyArg_ParseTuple(args, "OOd", &vals, &xvals, &xHigh) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( !(PySequence_Check(vals) && PySequence_Check(xvals) ) ) {
    PyErr_SetString(PyExc_ValueError, "not sequences.") ;
    return 0;
  }

  double x = 0.0;
  unsigned int k = 0;
  double v = 0.0;
  unsigned int const ll = PySequence_Size(xvals);

  while( x < xHigh ) {
    PyObject* const a = PySequence_Fast_GET_ITEM(vals, k);
    double const pop0 = PyFloat_AsDouble(a);
    
    if( k == ll ) {
      v += (xHigh - x) / pop0;
      break;
    } else {
      PyObject* const a = PySequence_Fast_GET_ITEM(vals, k+1);
      double pop1 = PyFloat_AsDouble(a);
      
      PyObject* const b = PySequence_Fast_GET_ITEM(xvals, k);
      double const x1 =  PyFloat_AsDouble(b);
      
      double dx = x1 - x;
        
      if (xHigh < x1) {
	double ndx = xHigh-x;
	pop1 = pop0 + (ndx/dx)*(pop1-pop0);
	dx = ndx;
      }
          
      if( pop0 == pop1 ) {
	v += dx/pop0;
      } else {
	double m = dx / (pop1 - pop0);
	v += m * log(pop1/pop0);
      }
      x = x1;
      ++k;
    }
  }
  
  return PyFloat_FromDouble(v);
}

static PyObject*
seqEvolve(PyObject*, PyObject* args)
{
  PyObject* pmat;
  PyObject* seq;
  //  double xHigh;

  if( !PyArg_ParseTuple(args, "OO", &pmat, &seq) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PyArray_Check(pmat) && PyArray_NDIM(pmat) == 2) ) {
     PyErr_SetString(PyExc_ValueError, "wrong args: not a 2d matrix") ;
     return 0;
  }
  
  if( ! PyList_Check(seq) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a list");
    return 0;
  }

  PyObject* nn[4];
  for(int k = 0; k < 4; ++k) {
    nn[k] = PyLong_FromLong(k);
  }
  
  const double* v = reinterpret_cast<const double*>(PyArray_GETPTR1(pmat, 0));
  double const p0[4] = {v[0],v[0]+v[1],v[0]+v[1]+v[2],1.0};
  double const p1[4] = {v[4],v[4]+v[5],v[4]+v[5]+v[6],1.0};
  double const p2[4] = {v[8],v[8]+v[9],v[8]+v[9]+v[10],1.0};
  double const p3[4] = {v[12],v[12]+v[13],v[12]+v[13]+v[14],1.0};
  const double* const p[4] = {p0,p1,p2,p3};
  
  int const nseq = PyList_Size(seq);

  PyObject** s = &PyList_GET_ITEM(seq, 0);
  for(int k = 0; k < nseq; ++k) {
    PyObject* o = s[k];
    long const nuc = PyInt_AsLong(o);
    const double* const pn = p[nuc];
    const double r = random()/double(RAND_MAX);
    int v;
    if( r < pn[1] ) {
      v = (r < pn[0]) ? 0 : 1;
    } else {
      v = (r < pn[2]) ? 2 : 3;
    }
    Py_DECREF(o);
    Py_INCREF(nn[v]);
    s[k] = nn[v];
  }

  Py_INCREF(seq);
  return seq;
}

//#include <iostream>
//using std::cout;
//using std::endl;

static PyObject*
seqsMinDiff(PyObject*, PyObject* args)
{
  PyObject* seqs1;
  PyObject* seqs2;

  if( !PyArg_ParseTuple(args, "OO", &seqs1, &seqs2) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PySequence_Check(seqs1) && PySequence_Check(seqs2)) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not sequences");
    return 0;
  }

  int const n1 = PySequence_Size(seqs1);
  int const n2 = PySequence_Size(seqs2);

  unsigned int gMin = PyString_Size(PySequence_Fast_GET_ITEM(seqs1,0))+1;
  
  for(int i1 = 0; i1 < n1; ++i1) {
     PyObject* const s1 = PySequence_Fast_GET_ITEM(seqs1, i1);
     const char* const s1c = PyString_AS_STRING(s1);

     for(int i2 = 0; i2 < n2; ++i2) {
       PyObject* const s2 = PySequence_Fast_GET_ITEM(seqs2, i2);
       const char* s2c = PyString_AS_STRING(s2);
       unsigned int count = 0;

       // cout << s2c << " vs. " << s1c << endl;
       
       for(const char* s = s1c;  *s ; ++s, ++s2c ) {
	 count += (*s != *s2c);
       }
       
       //cout << count << endl;
       if( count < gMin ) {
	 gMin = count;
       }
     }
  }

  //cout << gMin << endl;
  
  return PyLong_FromLong(gMin);
}



static inline bool
has(char const ch, const char* any) {
  for(/**/ ; *any; ++any) {
    if( *any == ch ) {
      return true;
    }
  }
  return false;
}
    
static inline int
_getStuff(const char* s, char const sep) {
  int e = 0;
  while( s[e] != sep || s[e-1] == '\\' ) {
    e += 1;
  }
  return e;
}

static int
findIndex(const char* s, char const ch, const char* stopAt)
{
  const char* s0 = s;

  while( *s != ch ) {
    if( has(*s, stopAt) ) {
      return s - s0;
    }
    ++s;
    if( ! *s ) {
      return -1;
    }
  }
  return s - s0;
}

static int
parseAttributes(const char* s, vector<PyObject*>& vals)
{
  int eat = 0;
  while( *s != ']' ) {
    if( *s == ',' ) {
      s += 1;
      eat += 1;
    }

    int const nameEnd = findIndex(s, '=', ",]\"{}");
    if( s[nameEnd] != '=' ) {
      return -1;
    }
    std::string name(s, nameEnd);
    s += nameEnd+1;
    eat += nameEnd+1;

    std::string v;
    if( *s == '"' ) {
      int const e = _getStuff(s+1, '"');
      v = std::string(s+1, e+1);
      s += e+2;
      eat += e+2;
    } else if( *s == '{' ) {
      int const e = _getStuff(s+1, '}');
      v = std::string(s+1, e+1);
      s += e+2;
      eat += e+2;
    } else {
      int const e = findIndex(s, ',', "]");
      if( e == -1 ) {
	return -1;
      }
      v = std::string(s, e);
      s += e;
      eat += e;
    }

    PyObject* o = PyTuple_New(2);
    PyTuple_SET_ITEM(o, 0, PyString_FromString(name.c_str()));
    PyTuple_SET_ITEM(o, 1, PyString_FromString(v.c_str()));
    vals.push_back(o);
  }
  return eat;
}


static inline int
skipSpaces(const char* txt)
{
  const char* s = txt;
  while( *s && isspace(*s) ) {
    ++s;
  }
  return s - txt;
}

static int
readSubTree(const char* txt, vector<PyObject*>& nodes)
{
  int n = skipSpaces(txt);
  txt += n;
  
  PyObject* vals = NULL;
  PyObject* nodeData;
  
  if( *txt == '(' ) {
    vector<int> subs;
    while( true ) {
      int n1 = readSubTree(txt+1, nodes);
      n += 1+n1;
      txt += 1+n1;
      subs.push_back(nodes.size()-1);

      n1 = skipSpaces(txt);
      n += n1;
      txt += n1;
      if( *txt == ',' ) {
        continue;
      }
      if( *txt == ')' ) {
	int const ns = subs.size();
	PyObject* psubs = PyList_New(ns);
	for(int k = 0; k < ns; ++k) { 
	  PyList_SET_ITEM(psubs, k, PyLong_FromLong(subs[k]));
	}
	
	nodeData = PyList_New(4);
	PyList_SET_ITEM(nodeData, 2, psubs);

	Py_INCREF(Py_None);
	PyList_SET_ITEM(nodeData, 0, Py_None);

	n += 1;
        txt += 1;
        break;
      }
      return -1; // error
    }
  } else {
    // a terminal
    const char* s = txt;
    while( ! isspace(*s) && *s != ':' && *s != '[' && *s != ','
	   && *s != '(' && *s != ')' && *s != ']' ) {
      ++s;
    }
    int const n1 = s - txt;

    nodeData = PyList_New(4);
    PyList_SET_ITEM(nodeData, 0, PyString_FromStringAndSize(txt, n1));

    Py_INCREF(Py_None);
    PyList_SET_ITEM(nodeData, 2, Py_None);

    n += n1;
    txt += n1;
  }

  {
    int n1 = skipSpaces(txt);
    txt += n1;
    n += n1;
  }
  
  if( *txt && *txt == '[' && txt[1] == '&' ) {
    vector<PyObject*> vs;
    int n1 = parseAttributes(txt+2, vs);
    if( n1 < 0 ) {
      return -1;
    }
    n1 += 3;
    n1 += skipSpaces(txt+n1);
    n += n1;
    txt += n1;

    if( vs.size() > 0 ) {
      vals = PyTuple_New(vs.size());
      for(int k = 0; k < vs.size(); ++k) {
	PyTuple_SET_ITEM(vals, k, vs[k]);
      }
    }
  }
  
  if( ! vals ) {
    Py_INCREF(Py_None);
    vals = Py_None;
  }
  PyList_SET_ITEM(nodeData, 3, vals);

  PyObject* branch = NULL;
  if( *txt && *txt == ':' ) {
    int n1 = skipSpaces(txt+1);
    n += n1 + 1;
    txt += 1+n1;
    
    char* endp;
    double b = strtod(txt, &endp);
    n1 = endp - txt;
    if( n1 == 0 ) {
      return -1;
    }
    
    branch = PyFloat_FromDouble(b);
    txt += n1;
    n += n1;
  } else {
    Py_INCREF(Py_None);
    branch = Py_None;
  }
  
  PyList_SET_ITEM(nodeData, 1, branch);
  nodes.push_back(nodeData);

  return n;
}

static PyObject*
parseSubTree(PyObject*, PyObject* args)
{
  const char* treeTxt;
  //PyObject* nodes;

  if( !PyArg_ParseTuple(args, "s", &treeTxt) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  vector<PyObject*> nodes;
  
  // if( ! PyList_Check(nodes) && PyList_Size(nodes) == 0 ) {
  //   PyErr_SetString(PyExc_ValueError, "wrong args: nodes should be an empty list") ;
  if( readSubTree(treeTxt, nodes) < 0 ) {
    // clean !!!
    for(int k = 0; k < nodes.size(); ++k) {
      Py_DECREF(nodes[k]);
    }
    //Py_DECREF(nodes);
    
    PyErr_SetString(PyExc_ValueError, "failed parsing.") ;
    return 0;
  }
  
  PyObject* n = PyTuple_New(nodes.size());
  for(int k = 0; k < nodes.size(); ++k) {
    PyTuple_SET_ITEM(n, k, nodes[k]);
  }
  
  return n;
}

static PyObject*
test(PyObject*, PyObject* args)
{
  const char* s;
  //PyObject* nodes;

  if( !PyArg_ParseTuple(args, "s", &s) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  // causes problems
  PyObject* nodes = PyList_New(0);
  PyList_Append(nodes, PyString_FromString(s));
  PyList_Append(nodes, PyString_FromString(s));

  //for k in range(10000) :   a = cchelp.test('a'*1000000)
// ---------------------------------------------------------------------------
// SystemError                               Traceback (most recent call last)

// /home/joseph/research/Projects/migration/data/case21s/04nl_b0p4_d0/m0_s0/<ipython console> in <module>()

// SystemError: ../Objects/listobject.c:290: bad argument to internal function

  // //  PyList_Append(nodes, PyString_FromStringAndSize(s, strlen(s)));
  // //  PyList_Append(nodes, PyString_FromStringAndSize(s, strlen(s)));


  // PyObject* nodes = PyTuple_New(2);
  // PyTuple_SET_ITEM(nodes, 0, PyString_FromString(s));
  // PyTuple_SET_ITEM(nodes, 1, PyString_FromString(s));

  // PyObject* nodes = PyList_New(2);
  // PyList_SET_ITEM(nodes, 0, PyString_FromString(s));
  // PyList_SET_ITEM(nodes, 1, PyString_FromString(s));

  return nodes;
}

static PyMethodDef cchelpMethods[] = {
  {"nonEmptyIntersection",  nonEmptyIntersection, METH_VARARGS,
   ""},

  {"demoLPintegrate",  demoLPintegrate, METH_VARARGS,
   ""},

  {"demoLPpopulation",  demoLPpopulation, METH_VARARGS,
   ""},

  {"seqevolve",  seqEvolve, METH_VARARGS,
   ""},

  {"seqsmindiff",  seqsMinDiff, METH_VARARGS,
   "Mimimum distance between all pairs of sequences from S1 x S2." 
   " Distance is total sum of mismatched characters."},

  {"parsetree",  parseSubTree, METH_VARARGS,
   ""},

  {"test",  test, METH_VARARGS,
   ""},
  
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initcchelp(void)
{
  import_array();
  /*PyObject* m = */ Py_InitModule("cchelp", cchelpMethods);
}