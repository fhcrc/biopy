#include <Python.h>
#include "numpy/arrayobject.h"
//#include <iostream>

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
  unsigned int const ll =  PySequence_Size(xvals);

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


static PyMethodDef cchelpMethods[] = {
  {"nonEmptyIntersection",  nonEmptyIntersection, METH_VARARGS,
   ""},

  {"demoLPintegrate",  demoLPintegrate, METH_VARARGS,
   ""},

  {"demoLPpopulation",  demoLPpopulation, METH_VARARGS,
   ""},

  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initcchelp(void)
{
  import_array();
  /*PyObject* m = */ Py_InitModule("cchelp", cchelpMethods);
}
