// This file is part of biopy.
// Copyright (C) 2010 Joseph Heled
// Author: Joseph Heled <jheled@gmail.com>
// See the files gpl.txt and lgpl.txt for copying conditions.


#include <Python.h>
#undef NDEBUG
#include <cassert>

#include <numpy/arrayobject.h>


#include <string>
using std::string;
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

  unsigned int const ll = PySequence_Size(xvals);

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
	
  // speed, not beuty 
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

  for(int k = 0; k < 4; ++k) {
    // remove creation reference count 
    Py_DECREF(nn[k]);
  }

  Py_INCREF(seq);
  return seq;
}

// #include <iostream>
// using std::cout;
// using std::endl;

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
  while( s[e] && (s[e] != sep || s[e-1] == '\\') ) {
    e += 1;
  }
  return s[e] ? e : -1;
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

static inline int
skipSpaces(const char* txt)
{
  const char* s = txt;
  while( *s && isspace(*s) ) {
    ++s;
  }
  return s - txt;
}

// trim spaces from both ends of string, return python string
static inline PyObject*
trimString(string const& s)
{
  const char* txt = s.c_str();
  int const n0 = skipSpaces(txt);
  int n1 = s.length()-1;
  while( n1 >= 0 && isspace(txt[n1]) ) {
    --n1;
  }
  return PyString_FromStringAndSize(txt + n0, n1+1-n0);
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
    string name(s, nameEnd);
    s += nameEnd+1;
    eat += nameEnd+1;

    string v;
    if( *s == '"' ) {
      int const e = _getStuff(s+1, '"');
      if( e < 0 ) {
	return -1;
      } else {
        v = string(s+1, e);
        s += e+2;
        eat += e+2;
      }
    } else if( *s == '{' ) {
      int const e = _getStuff(s+1, '}');
      if( e < 0 ) {
	return -1;
      } else {
	v = string(s+1, e);
	s += e+2;
	eat += e+2;
      }
    } else {
      int const e = findIndex(s, ',', "]");
      if( e == -1 ) {
	return -1;
      }
      v = string(s, e);
      s += e;
      eat += e;
    }

    PyObject* o = PyTuple_New(2);
    PyTuple_SET_ITEM(o, 0, trimString(name));
    PyTuple_SET_ITEM(o, 1, trimString(v));
    vals.push_back(o);
  }
  return eat;
}


static int
readSubTree(const char* txt, vector<PyObject*>& nodes)
{
  int eat = skipSpaces(txt);
  txt += eat;
  
  PyObject* vals = NULL;
  PyObject* const nodeData = PyList_New(4);
  
  if( *txt == '(' ) {
    vector<int> subs;
    while( true ) {
      int n1 = readSubTree(txt+1, nodes);
      eat += 1+n1;
      txt += 1+n1;
      subs.push_back(nodes.size()-1);

      n1 = skipSpaces(txt);
      eat += n1;
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
	
	PyList_SET_ITEM(nodeData, 2, psubs);

	Py_INCREF(Py_None);
	PyList_SET_ITEM(nodeData, 0, Py_None);

	eat += 1;
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

    PyList_SET_ITEM(nodeData, 0, PyString_FromStringAndSize(txt, n1));

    Py_INCREF(Py_None);
    PyList_SET_ITEM(nodeData, 2, Py_None);

    eat += n1;
    txt += n1;
  }

  {
    int const n1 = skipSpaces(txt);
    txt += n1;
    eat += n1;
  }

  string nodeTxt;
  while( *txt ) {
    if( *txt == '[' ) {
      if( txt[1] == '&' ) {
	vector<PyObject*> vs;
	int n1 = parseAttributes(txt+2, vs);
	if( n1 < 0 ) {
	  return -1;
	}
	n1 += 3;
	n1 += skipSpaces(txt+n1);
	eat += n1;
	txt += n1;

	if( vs.size() > 0 ) {
	  int b;
	  if( vals == NULL ) {
	    vals = PyTuple_New(vs.size());
	    b = 0;
	  } else {
	    b = PyTuple_Size(vals);
	    _PyTuple_Resize(&vals, b + vs.size());
	  }
	  for(unsigned int k = 0; k < vs.size(); ++k) {
	    PyTuple_SET_ITEM(vals, b+k, vs[k]);
	  }
	}
      } else {
	// skip comment
	int const e = _getStuff(txt+1, ']');
	if( e < 0 ) {
	  return -1;
	} else {
	  txt += e+2;
	  eat += e+2;
	}
      }
    } else {
      if( isspace(*txt) || isdigit(*txt) || has(*txt, ":.+-Ee") ) {
	nodeTxt.append(txt, 1);
	txt += 1;
	eat += 1;
      } else {
	break;
      }
    }
  }

  if( ! vals ) {
    Py_INCREF(Py_None);
    vals = Py_None;
  }
  PyList_SET_ITEM(nodeData, 3, vals);

  PyObject* branch = NULL;

  //std::cout << nodeTxt << std::endl;  
  const char* nTxt = nodeTxt.c_str();
  nTxt += skipSpaces(nTxt);
  if( *nTxt && *nTxt == ':' ) {
    int n1 = skipSpaces(nTxt+1);
    nTxt += 1+n1;
    
    //std::cout << nTxt << std::endl;  
    char* endp;
    double b = strtod(nTxt, &endp);
    n1 = endp - nTxt;
    if( n1 == 0 ) {
      return -1;
    }
    
    branch = PyFloat_FromDouble(b);
  } else {
    Py_INCREF(Py_None);
    branch = Py_None;
  }

  PyList_SET_ITEM(nodeData, 1, branch);
  nodes.push_back(nodeData);

  return eat;
  
#if 0
  if( *txt && *txt == '[' && txt[1] == '&' ) {
    vector<PyObject*> vs;
    int n1 = parseAttributes(txt+2, vs);
    if( n1 < 0 ) {
      return -1;
    }
    n1 += 3;
    n1 += skipSpaces(txt+n1);
    eat += n1;
    txt += n1;

    if( vs.size() > 0 ) {
      vals = PyTuple_New(vs.size());
      for(unsigned int k = 0; k < vs.size(); ++k) {
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
    eat += n1 + 1;
    txt += 1+n1;
    
    char* endp;
    double b = strtod(txt, &endp);
    n1 = endp - txt;
    if( n1 == 0 ) {
      return -1;
    }
    
    branch = PyFloat_FromDouble(b);
    txt += n1;
    eat += n1;
  } else {
    Py_INCREF(Py_None);
    branch = Py_None;
  }
  
  PyList_SET_ITEM(nodeData, 1, branch);
  nodes.push_back(nodeData);

  return eat;
#endif
}

static PyObject*
parseSubTree(PyObject*, PyObject* args)
{
  const char* treeTxt;

  if( !PyArg_ParseTuple(args, "s", &treeTxt) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  vector<PyObject*> nodes;
  
  if( readSubTree(treeTxt, nodes) < 0 ) {
    // clean !!!
    for(unsigned int k = 0; k < nodes.size(); ++k) {
      for(int i = 0; i < 4; ++i) {
	Py_DECREF(PySequence_GetItem(nodes[k], i));
      }
      Py_DECREF(nodes[k]);
    }
    
    PyErr_SetString(PyExc_ValueError, "failed parsing.") ;
    return 0;
  }
  
  PyObject* n = PyTuple_New(nodes.size());
  for(unsigned int k = 0; k < nodes.size(); ++k) {
    PyTuple_SET_ITEM(n, k, nodes[k]);
  }
  
  return n;
}


static PyObject*
effectiveSampleStep(PyObject*, PyObject* args)
{
  PyObject* data;

  if( !PyArg_ParseTuple(args, "O", &data) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! PySequence_Check(data) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a sequence");
    return 0;
  }

  int const nSamples = PySequence_Size(data);

  if( nSamples <= 3 ) {
    PyErr_SetString(PyExc_ValueError, "less that 4 samples");
    return 0;
  }
  
  int const maxLag = std::max(nSamples/3, 3);
  
  double* normalizedData = new double [nSamples];
  double sm = 0.0;
  
  for(int k = 0; k < nSamples; ++k) {
    PyObject* const dk = PySequence_Fast_GET_ITEM(data,k);
    double const g = PyFloat_AsDouble(dk);
    normalizedData[k] = g;
    sm += g;
  }

  sm /= nSamples;
  
  for(int k = 0; k < nSamples; ++k) {
    normalizedData[k] -= sm;
  }

  vector<double> gammaStats;
  
  double gammaStat0, gammaPrev;
  double varStat = 0.0;
  int lag;
  for(lag = 0; lag < maxLag; ++lag) {
    sm = 0.0;
    int const n = nSamples-lag;
    for(int k = 0; k < n; ++k) {
      sm += normalizedData[k] * normalizedData[k+lag];
    }
    double const gammaStat = sm / n;
    gammaStats.push_back(gammaStat);
      
    if( lag == 0 ) {
      gammaStat0 = gammaStat;
      varStat = gammaStat;
    } else {
      if( (lag & 0x1) == 0 ) {
	double const s = gammaPrev + gammaStat;
	if( s > 0 ) {
	  varStat += 2*s;
	} else {
	  break;
	}
      }
      gammaPrev = gammaStat;
    }
  }
  
  double const act = varStat / gammaStat0;

  // effective sample size
  // double const ess = nSamples / act;

  lag -= 1;

  //double const ess2 = lag > 0 ? nSamples/double(lag) : nSamples;
  double const act2 = lag > 0 ? lag : 1;

  double act1;

  double cvarStat = varStat + gammaStats[lag];
  // ppp = False
  // if ppp : print "clag", clag, cvarStat
  
  double const back = cvarStat * 0.05;
  //if ppp : print "back", back
  double totBack = 0;
  int const nn = gammaStats.size();
  
  assert(nn >= 3);
  while( lag > 0 ) {
    totBack += gammaStats[lag] + gammaStats[lag-1];
    //if ppp : print "totBack", totBack, clag
    
    lag -= 1;
    if( totBack >= back ) {
      break;
    }
  }
  
  if (lag == 0 ) {
    act1 = 1;
  } else {
    double cutBack = totBack - back;
    assert(lag+1 < nn);
    //double tot = gammaStats[lag] + gammaStats[lag+1];
    //assert cutBack <= tot, (cutBack, tot)
    //if ppp : print clag, cutBack, tot

    double const dif = (gammaStats[lag]-gammaStats[lag+1]);
    double const d2 = gammaStats[lag]*gammaStats[lag] - dif * cutBack;
    assert(d2 >= 0);
    double const sol = (gammaStats[lag] - sqrt(d2)) / dif;

    act1 = lag + sol;
  }
  
  delete [] normalizedData;

  PyObject* t = PyTuple_New(3);
  PyTuple_SET_ITEM(t, 0, PyFloat_FromDouble(act));
  PyTuple_SET_ITEM(t, 1, PyFloat_FromDouble(act1));
  PyTuple_SET_ITEM(t, 2, PyFloat_FromDouble(act2));
  
  return t;

  //  return PyFloat_FromDouble(ess);
}

static PyMethodDef cchelpMethods[] = {
  {"nonEmptyIntersection",  nonEmptyIntersection, METH_VARARGS,
   ""},

  {"demoLPintegrate",  demoLPintegrate, METH_VARARGS,
   ""},

  {"demoLPpopulation",  demoLPpopulation, METH_VARARGS,
   ""},

  {"seqevolve",  seqEvolve, METH_VARARGS,
   "(P-matrix, seqeunce). Evolve sequence according to P matrix."
   " Result changes sequence in place."},

  {"seqsmindiff",  seqsMinDiff, METH_VARARGS,
   "Mimimum distance between all pairs of sequences from S1 x S2." 
   " Distance is total sum of mismatched characters."},

  {"parsetree",  parseSubTree, METH_VARARGS,
   ""},

  {"effectiveSampleStep",  effectiveSampleStep, METH_VARARGS,
   ""},

  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initcchelp(void)
{
  import_array();
  /*PyObject* m = */ Py_InitModule("cchelp", cchelpMethods);
}
