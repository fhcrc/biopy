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
    double const vk = darrayElement(vals, k);
    double const vkp1 =  darrayElement(vals, k+1);;
    if( width == 0.0 ) {
      pop = vk;
    } else {
      pop = vk + (t/width) * (vkp1 - vk);
    }
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

#include <iostream>
using std::cout;
using std::endl;

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
}



static PyObject*
sumNonIntersect(PyObject*, PyObject* args)
{
  PyObject* lu2;
  double l1,u1;
  
  if( !PyArg_ParseTuple(args, "ddO", &l1, &u1, &lu2) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! PySequence_Check(lu2) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a sequence");
    return 0;
  }

  double tot = 0.0;
  double const b1 = u1-l1;
  
  int const nSeq = PySequence_Size(lu2);
  
  for(int k = 0; k < nSeq; ++k) {
    PyObject* const l2u2 = PySequence_Fast_GET_ITEM(lu2, k);

    double const l2 = PyFloat_AsDouble( PySequence_Fast_GET_ITEM(l2u2,0) );
    double const u2 = PyFloat_AsDouble( PySequence_Fast_GET_ITEM(l2u2,1) );
    double const b1b2 = b1 + (u2-l2);
    tot += std::min(b1b2 + 2*(std::max(l1,l2) - std::min(u1,u2)), b1b2);
  }

  return PyFloat_FromDouble(tot);
}

static PyObject*
sumNonIntersectDer(PyObject*, PyObject* args)
{
  PyObject* lu2;
  double l1,u1;
  
  if( !PyArg_ParseTuple(args, "ddO", &l1, &u1, &lu2) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! PySequence_Check(lu2) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a sequence");
    return 0;
  }

  double tot = 0.0;
  double const b1 = u1-l1;
  int dvdl1 = 0, dvdu1 = 0;
    
  int const nSeq = PySequence_Size(lu2);
  
  for(int k = 0; k < nSeq; ++k) {
    PyObject* const l2u2 = PySequence_Fast_GET_ITEM(lu2, k);

    double const l2 = PyFloat_AsDouble( PySequence_Fast_GET_ITEM(l2u2,0) );
    double const u2 = PyFloat_AsDouble( PySequence_Fast_GET_ITEM(l2u2,1) );
    double const b1b2 = b1 + (u2-l2);
    tot += std::min(b1b2 + 2*(std::max(l1,l2) - std::min(u1,u2)), b1b2);

    dvdl1 += ( l2 <= l1 && l1 <= u2 ) ? 1 : -1;
    dvdu1 += ( l2 <= u1 && u1 <= u2 ) ? -1 : 1;
  }

  PyObject* o = PyTuple_New(3);
  PyTuple_SET_ITEM(o, 0, PyFloat_FromDouble(tot));
  PyTuple_SET_ITEM(o, 1, PyFloat_FromDouble(dvdl1));
  PyTuple_SET_ITEM(o, 2, PyFloat_FromDouble(dvdu1));

  return o;
}

static PyObject*
varianceAndDerive(PyObject*, PyObject* args)
{
  PyObject* bsnr;
  PyObject* bsr;
  double a;
  
  if( !PyArg_ParseTuple(args, "OOd", &bsnr, &bsr, &a) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PySequence_Check(bsnr) && PySequence_Check(bsr)
	 && PySequence_Size(bsr) == 2 )) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a sequence");
    return 0;
  }

  uint const n1 = PySequence_Size(bsnr);
  uint const n = n1 + 2;
  
  PyObject* bsr0 = PySequence_Fast_GET_ITEM(bsr,0);
  PyObject* bsr1 = PySequence_Fast_GET_ITEM(bsr,1);
  double const bsr00 = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr0,0));
  double const bsr10 = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr1,0));

  double b2 = (PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr0,1)) +
	       PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr1,1)));
  double rl = a/bsr00;
  double rr = (b2 - a)/bsr10;
  
  //uint ne = PySequence_Size(PySequence_Fast_GET_ITEM(bs,0));

  double sum = 0.0, sum2 = 0.0;
  double r[n];
  double b[n];
  
  for(uint k = 0; k < n1; ++k) {
    PyObject* const x = PySequence_Fast_GET_ITEM(bsnr,k);
    double const br = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(x, 0));
    double const s = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(x, 1));
    double const v = s/br;
    r[k] = v;
    b[k] = br;
    sum += v;
    sum2 += v*v;
  }
  sum += rl + rr;
  sum2 += rl*rl + rr*rr;
  
  double const avg = sum / n;

  PyObject* ret = PyFloat_FromDouble(sum2/n - avg*avg);
  
  if( PySequence_Size(bsr0) == 3 ) {
    r[n1] = rl;
    r[n1+1] = rr;
    b[n1] = bsr00;
    b[n1+1] = bsr10;
    
    double dfdir[n];
    double const f2 = 2./(n*n);
    for(uint k = 0; k < n; ++k) {
      dfdir[k] = f2 * (n*r[k] - sum) * -(r[k]/b[k]);
    }

    PyObject* a[n];
    for(uint k = 0; k < n1; ++k) {
      PyObject* const x = PySequence_Fast_GET_ITEM(bsnr,k);
      a[k] = PySequence_Fast_GET_ITEM(x, 2);
    }
    a[n1] = PySequence_Fast_GET_ITEM(bsr0,2);
    a[n1+1] = PySequence_Fast_GET_ITEM(bsr1,2);
       
    uint const m = PySequence_Size(a[0]);
    PyObject* p = PyTuple_New(m+1);
    for(uint j = 0; j < m; ++j) {
      double v(0);
      for(uint k = 0; k < n; ++k) {
	PyObject* const akj = PySequence_Fast_GET_ITEM(a[k], j);
	//assert ( PyFloat_Check(akj) ||  PyInt_Check(akj) ) ;
	v += dfdir[k] * PyFloat_AsDouble(akj);
      }
      PyTuple_SET_ITEM(p, j, PyFloat_FromDouble(v));
    }

    double const pa = f2 * ( (n*r[n1] - sum)/bsr00 -  (n*r[n1+1] - sum)/bsr10 );
    PyTuple_SET_ITEM(p, m, PyFloat_FromDouble(pa));
    
    PyObject* r = PyTuple_New(2);
    PyTuple_SET_ITEM(r, 0 , ret);
    PyTuple_SET_ITEM(r, 1 , p);
    ret = r;
  }
  return ret;
}



static PyObject*
normedVarianceAndDeriv(PyObject*, PyObject* args)
{
  PyObject* bsnr;
  PyObject* bsr;
  double a;
  
  if( !PyArg_ParseTuple(args, "OOd", &bsnr, &bsr, &a) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PySequence_Check(bsnr) && PySequence_Check(bsr) &&  PySequence_Size(bsr) == 2 )) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a sequence");
    return 0;
  }

  uint const n1 = PySequence_Size(bsnr);
  uint const n = n1 + 2;
  
  PyObject* const bsr0 = PySequence_Fast_GET_ITEM(bsr,0);
  PyObject* const bsr1 = PySequence_Fast_GET_ITEM(bsr,1);
  double const bsr00 = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr0,0));
  double const bsr10 = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr1,0));

  double const b2 = (PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr0,1)) +
	       PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr1,1)));
  double const rl = a/bsr00;
  double const rr = (b2 - a)/bsr10;
  
  double sum = 0.0, sum2 = 0.0;
  double r[n];
  double b[n];
  
  for(uint k = 0; k < n1; ++k) {
    PyObject* const x = PySequence_Fast_GET_ITEM(bsnr,k);
    double const br = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(x, 0));
    double const s = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(x, 1));
    double const v = s/br;
    r[k] = v;
    b[k] = br;
    sum += v;
    sum2 += v*v;
  }
  sum += rl + rr;
  sum2 += rl*rl + rr*rr;
  
  double const avg = sum / n;

  double const iavg2 = 1.0/(avg*avg);
  double const val = (sum2/n) * iavg2 - 1;
  
  PyObject* ret = PyFloat_FromDouble(val);
  
  if( PySequence_Size(bsr0) == 3 ) {
    r[n1] = rl;
    r[n1+1] = rr;
    b[n1] = bsr00;
    b[n1+1] = bsr10;

    double const coeff = 2/(sum*sum);
    double const a1 = n * coeff;
    double const b1 = coeff * sum * (1 + val);
    
    double dfdr[n];
    for(uint k = 0; k < n; ++k) {
      dfdr[k] = a1 * r[k] - b1;
    }
    double const dfdr_l = dfdr[n-2];
    double const dfdr_r = dfdr[n-1];

    for(uint k = 0; k < n; ++k) {
      dfdr[k] *= -(r[k]/b[k]);
    }
    
    PyObject* a[n];
    for(uint k = 0; k < n1; ++k) {
      PyObject* const x = PySequence_Fast_GET_ITEM(bsnr,k);
      a[k] = PySequence_Fast_GET_ITEM(x, 2);
    }
    a[n1] = PySequence_Fast_GET_ITEM(bsr0,2);
    a[n1+1] = PySequence_Fast_GET_ITEM(bsr1,2);
       
    uint const m = PySequence_Size(a[0]);
    PyObject* p = PyTuple_New(m+1);
    for(uint j = 0; j < m; ++j) {
      double v(0);
      for(uint k = 0; k < n; ++k) {
	PyObject* const akj = PySequence_Fast_GET_ITEM(a[k], j);
	//assert ( PyFloat_Check(akj) ||  PyInt_Check(akj) ) ;
	v += dfdr[k] * PyFloat_AsDouble(akj);
      }
      PyTuple_SET_ITEM(p, j, PyFloat_FromDouble(v));
    }

    double const pa = dfdr_l/bsr00 - dfdr_r/bsr10;
    PyTuple_SET_ITEM(p, m, PyFloat_FromDouble(pa));
    
    PyObject* r = PyTuple_New(2);
    PyTuple_SET_ITEM(r, 0 , ret);
    PyTuple_SET_ITEM(r, 1 , p);
    ret = r;
  }
  return ret;
}

static PyObject*
normedVarianceAndDerivNew(PyObject*, PyObject* args)
{
  PyObject* bsnr;
  PyObject* bsr;
  double a;
  int variant;
  
  if( !PyArg_ParseTuple(args, "OOdi", &bsnr, &bsr, &a, &variant) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PySequence_Check(bsnr) && PySequence_Check(bsr) &&  PySequence_Size(bsr) == 2 )
      && (0 <= variant && variant <= 8) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a sequence");
    return 0;
  }

  uint const n1 = PySequence_Size(bsnr);
  uint const n = n1 + 2;
  
  PyObject* const bsr0 = PySequence_Fast_GET_ITEM(bsr,0);
  PyObject* const bsr1 = PySequence_Fast_GET_ITEM(bsr,1);
  double const bsr00 = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr0,0));
  double const bsr10 = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr1,0));

  double const b2 = (PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr0,1)) +
	       PyFloat_AsDouble(PySequence_Fast_GET_ITEM(bsr1,1)));
  double const rl = a/bsr00;
  double const rr = (b2 - a)/bsr10;
  
  double sum = 0.0, sum2 = 0.0;
  double r[n];
  double b[n];
  
  for(uint k = 0; k < n1; ++k) {
    PyObject* const x = PySequence_Fast_GET_ITEM(bsnr,k);
    double const br = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(x, 0));
    double const s = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(x, 1));
    double const v = s/br;
    r[k] = v;
    b[k] = br;
    sum += v;
    sum2 += v*v;
  }
  r[n1] = rl;
  r[n1+1] = rr;
  b[n1] = bsr00;
  b[n1+1] = bsr10;

  sum += rl + rr;
  sum2 += rl*rl + rr*rr;
  
  double const avg = sum / n;
  double const avg2 = avg*avg;
  double const iavg2 = 1.0/avg2;
  
  double val;
  switch( variant ) {
    case 0 :
    {
      val = (sum2/n) * iavg2 - 1;
      break;
    }
    case 1 :
    {
      double const a = (avg-1);
      val = sqrt(sum2 * iavg2 - n) + a*a;
      break;
    }
    case 2 :
    {
      double const a = (avg-1);
      val = (sum2 * iavg2 - n) + a*a;
      break;
    }
    case 3 :
    {
      double const a = (avg-1);
      val = (sum2 * iavg2 - n) + fabs(a);
      break;
    }
    case 4 :
    {
      double const a = (avg-1);
      val = sqrt(sum2 * iavg2 - n) + fabs(a);
      break;
    }
    case 5 :
    {
      double const a = (avg-1);
      val = (sum2 * iavg2 - n) + n * a * a;
      break;
    }
    case 6 :
    {
      double const a = (avg-1);
      val = (sum2 * iavg2 - n) + n * n * a * a;
      break;
    }
    case 7 : case 8 :
    {
      double lx[n];
      double slx = 0;
      for(uint k = 0; k < n; ++k) {
	double x = r[k];
	if( x <= 0 ) {
	  x = 1e-10;
	}
	lx[k] = log(x);
	slx += lx[k];
      }
      double const av = slx/n;
      double const a = (avg-1);
      val = n * a * a;
      if( variant == 8 ) {
	val *= n;
      }
      for(uint k = 0; k < n; ++k) {
	double const v = lx[k] - av;
	val += v*v;
      }
      break;
    }
  }
    
  PyObject* ret = PyFloat_FromDouble(val);
  
  if( PySequence_Size(bsr0) == 3 ) {

    double dfdr[n];
    if( variant == 0 ) {
      double const coeff = 2/(sum*sum);
      double const a1 = n * coeff;
      double const b1 = coeff * sum * (1 + val);
    
      for(uint k = 0; k < n; ++k) {
	dfdr[k] = a1 * r[k] - b1;
      }
    } else if( 7 <= variant && variant <= 8 ) {
      double lx[n];
      double slx = 0;
      for(uint k = 0; k < n; ++k) {
	double x = r[k];
	if( x <= 0 ) {
	  x = 1e-10;
	}
	lx[k] = log(x);
	slx += lx[k];
      }
      double const av = slx/n;
      double const c1 = (variant == 7) ? 2 * (avg - 1) : 2 * (avg - 1) * n;
      double const a1 = (2.0 * (n-2))/n;
      for(uint k = 0; k < n; ++k) {
	dfdr[k] = ((lx[k] - av) * a1)/r[k] + c1;
      }
    } else {
      double const m = sum * avg - sum2;
      double const sumu = (2*m)/(n * avg2*avg);

      double a1 = 2 * iavg2;
      double b1;
      
      switch( variant ) {
	case 1:
	{
	  double const p = 1/(2*sqrt(sum2 * iavg2 - n));
	  a1 *= p;
	  double const c1 = (2.0/n) * (avg - 1);
	  b1 = p * sumu - a1 * avg + c1;
	  break;
	}
	case 2:
	{
	  double const c1 = (2.0/n) * (avg - 1);
	  
	  b1 = sumu - a1 * avg + c1;
	  break;
	}
	case 3:
	{
	  double const ep = 1e-6;
	  double const c1 = fabs(avg-1) < ep ? 0 : (avg > 1 ? (1./n) : (-1./n));
	  
	  b1 = sumu - a1 * avg + c1;
	  break;
	}
	case 4:
	{
	  double const p = 1/(2*sqrt(sum2 * iavg2 - n));
	  a1 *= p;

	  double const ep = 1e-6;
	  double const c1 = fabs(avg-1) < ep ? 0 : (avg > 1 ? (1./n) : (-1./n));
	  
	  b1 = p * sumu - a1 * avg + c1;
	  break;
	}
	case 5:
	{
	  double const c1 = (2.0) * (avg - 1);
	  
	  b1 = sumu - a1 * avg + c1;
    	  break;
	}
	case 6:
	{
	  double const c1 = (2.0*n) * (avg - 1);
	  
	  b1 = sumu - a1 * avg + c1;
	  break;
	}
      }
      
      for(uint k = 0; k < n; ++k) {
	dfdr[k] = a1 * r[k] + b1;
      }
    }
    
    double const dfdr_l = dfdr[n-2];
    double const dfdr_r = dfdr[n-1];

    for(uint k = 0; k < n; ++k) {
      dfdr[k] *= -(r[k]/b[k]);
    }
    
    PyObject* a[n];
    for(uint k = 0; k < n1; ++k) {
      PyObject* const x = PySequence_Fast_GET_ITEM(bsnr,k);
      a[k] = PySequence_Fast_GET_ITEM(x, 2);
    }
    a[n1] = PySequence_Fast_GET_ITEM(bsr0,2);
    a[n1+1] = PySequence_Fast_GET_ITEM(bsr1,2);
       
    uint const m = PySequence_Size(a[0]);
    PyObject* p = PyTuple_New(m+1);
    for(uint j = 0; j < m; ++j) {
      double v(0);
      for(uint k = 0; k < n; ++k) {
	PyObject* const akj = PySequence_Fast_GET_ITEM(a[k], j);
	//assert ( PyFloat_Check(akj) ||  PyInt_Check(akj) ) ;
	v += dfdr[k] * PyFloat_AsDouble(akj);
      }
      PyTuple_SET_ITEM(p, j, PyFloat_FromDouble(v));
    }

    double const pa = dfdr_l/bsr00 - dfdr_r/bsr10;
    PyTuple_SET_ITEM(p, m, PyFloat_FromDouble(pa));
    
    PyObject* r = PyTuple_New(2);
    PyTuple_SET_ITEM(r, 0 , ret);
    PyTuple_SET_ITEM(r, 1 , p);
    ret = r;
  }
  return ret;
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

  {"sumNonIntersect",  sumNonIntersect, METH_VARARGS,
   ""},

  {"sumNonIntersectDer",  sumNonIntersectDer, METH_VARARGS,
   ""},

  {"effectiveSampleStep",  effectiveSampleStep, METH_VARARGS,
   ""},

  {"varianceAndDerive", varianceAndDerive, METH_VARARGS,
   ""},
  
  {"normedVarianceAndDeriv", normedVarianceAndDeriv, METH_VARARGS,
   ""},

  {"normedVarianceAndDerivNew", normedVarianceAndDerivNew, METH_VARARGS,
   ""},
  
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initcchelp(void)
{
  import_array();
  /*PyObject* m = */ Py_InitModule("cchelp", cchelpMethods);
}
