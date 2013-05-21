// This file is part of biopy.
// Copyright (C) 2010 Joseph Heled
// Author: Joseph Heled <jheled@gmail.com>
// See the files gpl.txt and lgpl.txt for copying conditions.


#include <Python.h>
// keeps asserts
#undef NDEBUG
#include <cassert>

#include <random>

struct PatchTrace {
  PatchTrace(double t, uint pl, uint pt, uint s, PatchTrace* prev);
  
  double pTime;
  uint  fromPlot : 16;
  uint  fromPatch : 16;
  uint	speciesIndex;
  PatchTrace*	prev;
  int	refCount;
  PyObject* asp;
  
  PyObject* asPyObject(void);
};


PatchTrace::PatchTrace(double t, uint pl, uint pt, uint s, PatchTrace* _prev) :
  pTime(t),
  fromPlot(pl),
  fromPatch(pt),
  speciesIndex(s),
  prev(_prev),
  refCount(1),
  asp(0)
{}

PyObject*
PatchTrace::asPyObject(void)
{
  if( asp ) {
    Py_INCREF(asp);
    return asp;
  }
  asp = PyTuple_New(4);
  PyTuple_SET_ITEM(asp, 0, PyFloat_FromDouble(pTime));
  PyObject* u2 = PyTuple_New(2);
  PyTuple_SET_ITEM(u2, 0, PyInt_FromLong(fromPlot));
  PyTuple_SET_ITEM(u2, 1, PyInt_FromLong(fromPatch));
  PyTuple_SET_ITEM(asp, 1, u2);
  PyTuple_SET_ITEM(asp, 2, PyInt_FromLong(speciesIndex));
  PyObject* p;
  if( prev ) {
    p = prev->asPyObject();
  } else {
    Py_INCREF(Py_None);
    p = Py_None;
  }
  PyTuple_SET_ITEM(asp, 3, p);
  return asp;
}

struct Patch {
  uint	speciesIndex;
};
  
class MetaCommunity {
public:
  MetaCommunity(uint nPlots, uint plotSize, double timeStamp = 0);
  ~MetaCommunity();

  uint const nPlots;
  uint const plotSize;

  double getTime(void) const { return timeStamp; }
  void setTime(double t) { timeStamp = t; }
  unsigned long nIndividuals(void) const { return nPlots * plotSize; }

  Patch& get(uint np, uint ni) { return community[np][ni]; }
  Patch const& get(uint np, uint ni) const { return community[np][ni]; }

  PyObject* asPyObject(void) const;

  void setSpecies(uint np, uint ni, uint s);

  uint getLastSpecies(void) const { return lastSpecies; }
  
private:
  Patch** community;
  double  timeStamp;
  uint    lastSpecies;
};

MetaCommunity::MetaCommunity(uint _nPlots, uint _plotSize, double _timeStamp) :
  nPlots(_nPlots),
  plotSize(_plotSize),
  timeStamp(_timeStamp),
  lastSpecies(0)
{
  community = new Patch* [nPlots];
  for(uint np = 0; np < nPlots; ++np) {
    community[np] = new Patch [plotSize];
    for(uint ni = 0; ni < plotSize; ++ni) {
      Patch& p = community[np][ni];
      p.speciesIndex = 0;
    }
  }
}

MetaCommunity::~MetaCommunity()
{
   for(uint np = 0; np < nPlots; ++np) {
     delete [] community[np];
   }
   delete [] community;
}

inline void
MetaCommunity::setSpecies(uint np, uint ni, uint s)
{
  get(np,ni).speciesIndex = s;
  if( s > lastSpecies ) {
    lastSpecies = s;
  }
}

PyObject*
MetaCommunity::asPyObject(void) const
{
  PyObject* o = PyTuple_New(nPlots);
  for(uint np = 0; np < nPlots; ++np) {
    PyObject* p = PyTuple_New(plotSize);
    for(uint ni = 0; ni < plotSize; ++ni) {
      PyTuple_SET_ITEM(p, ni, PyInt_FromLong(get(np,ni).speciesIndex));
    }
    PyTuple_SET_ITEM(o, np, p);
  }
  return o;
}

class TracedMetaCommunity : public MetaCommunity {
public:
  TracedMetaCommunity(uint nPlots, uint plotSize, double timeStamp = 0);
  ~TracedMetaCommunity();
  
  void replace(uint np, uint nt, uint np1, uint nt1, double t, uint sx);

  void setSpecies(uint np, uint ni, uint s);

  const PatchTrace*	cleanUp(void);

  PyObject* asPyObject(void);

  const PatchTrace*	ca(uint p0, uint i0, uint p1, uint i1) const;

  void 	optTraces(void);
  
private:
  PatchTrace*** trace;
  mutable PatchTrace founder;

  PatchTrace*		ca(void) const;
  const PatchTrace*	ca(uint n) const;
};

TracedMetaCommunity::TracedMetaCommunity(uint nPlots, uint plotSize, double timeStamp) :
  MetaCommunity(nPlots, plotSize, timeStamp),
  founder(-1, -1, -1, -1, 0)
{
  trace = new PatchTrace** [nPlots];
  for(uint np = 0; np < nPlots; ++np) {
    trace[np] = new PatchTrace* [plotSize];
    for(uint ni = 0; ni < plotSize; ++ni) {
      PatchTrace* pt = new PatchTrace(timeStamp, np, ni, 0, &founder);
      ++founder.refCount;
      trace[np][ni] = pt;
    }
  }
}

TracedMetaCommunity::~TracedMetaCommunity()
{
  for(uint np = 0; np < nPlots; ++np) {
    for(uint ni = 0; ni < plotSize; ++ni) {
      PatchTrace* p = trace[np][ni];
      PatchTrace* pp = p->prev;
      while( pp && pp != &founder ) {
	if( pp->refCount > 1 ) {
	  -- pp->refCount;
	  break;
	}
	PatchTrace* px = pp->prev;
	delete pp;
	pp = px;
      }
      delete p;
    }
    delete [] trace[np];
  }
  delete [] trace;
}

void
TracedMetaCommunity::setSpecies(uint np, uint ni, uint s)
{
  MetaCommunity::setSpecies(np, ni, s);
  trace[np][ni]->speciesIndex = s;
}

void
TracedMetaCommunity::replace(uint np, uint nt, uint np1, uint nt1, double t, uint sx)
{
  Patch& p = get(np, nt);
  p.speciesIndex = sx;
  PatchTrace* const prev = trace[np1][nt1];
  ++ prev->refCount;
  
  PatchTrace*& cur = trace[np][nt];
  -- cur->refCount;
  {
    PatchTrace* p = cur;
    while( p->refCount == 0 ) {
      PatchTrace* p1 = p->prev;
      delete p;
      p = p1;
      assert( p && p->refCount > 0 );
      --p->refCount;
      if( p == &founder ) {
	break;
      }
      // if( ! p || p == &founder ) {
      // 	// founder
      // 	break;
      // }
      // assert( p->refCount > 0 );
      // --p->refCount;
    }
  }
    
  PatchTrace* n = new PatchTrace(t, np1, nt1, sx, prev);

  cur = n;
}

#if 0 

static const PatchTrace*
commonAnc(const PatchTrace* x, const PatchTrace* y)
{
  while( !(x->pTime == y->pTime &&
	   x->fromPlot == y->fromPlot &&
	   x->fromPatch == y->fromPatch) )  {
    if( x->pTime > y->pTime ) {
      x = x->prev;
    } else {
      y = y->prev;
    }
  }
  assert( x == y ) ;
  return x;
}

#else

static const PatchTrace*
commonAnc(const PatchTrace* x, const PatchTrace* y)
{
  while( x != y ) {
    assert( x && y && !(x->pTime == y->pTime &&
			x->fromPlot == y->fromPlot &&
			x->fromPatch == y->fromPatch) );
    if( x->pTime > y->pTime ) {
      x = x->prev;
    } else {
      y = y->prev;
    }
  }
  return x;
}

#endif


const PatchTrace*
TracedMetaCommunity::ca(uint n) const
{
  PatchTrace** pn = trace[n];
  const PatchTrace* x = pn[0];
  for(uint k = 1; k < plotSize; ++k) {
    x = commonAnc(x, pn[k]);
  }
  return x;
}

PatchTrace*
TracedMetaCommunity::ca(void) const
{
  const PatchTrace* x = ca(0);
  for(uint n = 1; n < nPlots; ++n) {
    if( x == &founder ) {
      break;
    }
    x = commonAnc(x, ca(n));
  }
  return const_cast<PatchTrace*>(x);
}

const PatchTrace*
TracedMetaCommunity::ca(uint p0, uint i0, uint p1, uint i1) const
{
  return commonAnc(trace[p0][i0], trace[p1][i1]);
}

#if !defined(NDEBUG)
static bool
onPath(const PatchTrace* p, const PatchTrace* c)
{
  while( p != c && p->prev ) {
    p = p->prev;
  }
  return p == c;
}
#endif

const PatchTrace*
TracedMetaCommunity::cleanUp(void)
{
  PatchTrace* const can = ca();
#if !defined(NDEBUG)
  {
    for(uint np = 0; np < nPlots; ++np) {
      for(uint ni = 0; ni < plotSize; ++ni) {
	assert( onPath(trace[np][ni], can) );
      }
    }
  }
#endif
  
  {
    PatchTrace* p = can->prev;
    while( p && p->prev ) {
      assert(p->refCount == 1);
      p = p->prev;
    }
  }
  
  PatchTrace* p = can->prev;
  can->prev = 0;
  while( p && p->prev ) {
    assert(p->refCount == 1);
    PatchTrace* c = p->prev;
    delete p;
    p = c;
  }
  return can;
}


PyObject*
TracedMetaCommunity::asPyObject(void)
{
  const PatchTrace* const ca = cleanUp();
  if( ca == &founder ) {
    founder.asp = PyTuple_New(4);
    PyTuple_SET_ITEM(founder.asp, 0, PyFloat_FromDouble(-1));
    PyObject* u2 = PyTuple_New(2);
    PyTuple_SET_ITEM(u2, 0, PyInt_FromLong(-1));
    PyTuple_SET_ITEM(u2, 1, PyInt_FromLong(-1));
    PyTuple_SET_ITEM(founder.asp, 1, u2);
    PyTuple_SET_ITEM(founder.asp, 2, PyInt_FromLong(-1));
    Py_INCREF(Py_None);
    PyTuple_SET_ITEM(founder.asp, 3, Py_None);
  }
  
  PyObject* p2 = PyTuple_New(2);
  PyTuple_SET_ITEM(p2, 0, MetaCommunity::asPyObject());

  PyObject* o = PyTuple_New(nPlots);
  for(uint np = 0; np < nPlots; ++np) {
    PyObject* p = PyTuple_New(plotSize);
    for(uint ni = 0; ni < plotSize; ++ni) {
      PatchTrace* t = trace[np][ni];
      PyObject* x = t->asPyObject();
      PyTuple_SET_ITEM(p, ni, x);
    }
    PyTuple_SET_ITEM(o, np, p);
  }
  
  PyTuple_SET_ITEM(p2, 1, o);

  const_cast<PatchTrace*>(ca)->asp = 0;
  
  for(uint np = 0; np < nPlots; ++np) {
    for(uint ni = 0; ni < plotSize; ++ni) {
      PatchTrace* t = trace[np][ni];
      while( t->asp ) {
	t->asp = 0;
	t = t->prev;
	assert( t ) ;
      }
    }
  }
  assert( founder.asp == 0 ) ;
  
  return p2;
}

void
TracedMetaCommunity::optTraces(void)
{
  for(uint np = 0; np < nPlots; ++np) {
    for(uint ni = 0; ni < plotSize; ++ni) {
      PatchTrace* p = trace[np][ni];
      while( p != &founder && p->prev ) {
	if( p->prev != &founder && p->prev->refCount == 1 &&
	    p->prev->prev != &founder &&  p->prev->prev->refCount == 1 ) {
	  PatchTrace* x = p->prev;
	  p->prev = x->prev;
	  delete x;
	} else {
	  p = p->prev;
	}
      }
    }
  }
}

class MetaSimulator {
public:
  MetaSimulator(double mu, double lam, double rho);
  ~MetaSimulator() {}

  double advance(double targetTime, MetaCommunity& com);

  double advance(double targetTime, TracedMetaCommunity& com, double minCAtime, double cleanEvery);
  
private:
  double mu;
  double lam;
  double rho;
  
  double pLocal;
  double pLocalOrImi;
};

MetaSimulator::MetaSimulator(double _mu, double _lam, double _rho) :
  mu(_mu),
  lam(_lam),
  rho(_rho)
{
  double p3[3];
  p3[0] = 1-(lam+rho)/mu;
  p3[1] = rho/mu;
  p3[2] = lam/mu;

  pLocal = p3[0];
  pLocalOrImi = p3[0]+p3[1];
}
  
static std::mt19937 randomizer;

uint maxSpecies(MetaCommunity const& com) {
  uint s = 0;
  for(uint np = 0; np < com.nPlots; ++np) {
    for(uint ni = 0; ni < com.plotSize; ++ni) {
      s = std::max(com.get(np,ni).speciesIndex, s);
    }
  }
  return s;
}

double
MetaSimulator::advance(double targetTime, MetaCommunity& com)
{
#if !defined(NDEBUG)
  uint const nPlots = com.nPlots;
#endif
  uint const plotSize = com.plotSize;
  uint const nIndividuals = com.nIndividuals();

  uint lastSpecies = maxSpecies(com);
  
  double const dr = nIndividuals * mu;
  std::exponential_distribution<double> timeToNextDeath(dr);
  std::uniform_int_distribution<int> pickGlobIndividual(0, nIndividuals-1);
  std::uniform_int_distribution<int> pickLocalIndividual(0, plotSize-1);
  std::uniform_real_distribution<double> zeroOne(0.0, 1.0);

  double curTime = com.getTime();

  while( curTime <= targetTime ) {
    double const deltaT = timeToNextDeath(randomizer);
    curTime += deltaT;
    uint const k = pickGlobIndividual(randomizer);
    uint const np = k / plotSize;         assert(0 <= np && np < nPlots);
    uint const nt = k - np * plotSize;    assert(0 <= nt && nt < plotSize );
    double const r = zeroOne(randomizer);

    uint np1,nt1,sx;
    
    if( r < pLocal ) {
      uint const i = pickLocalIndividual(randomizer);
      sx = com.get(np, i).speciesIndex;
      np1 = np;
      nt1 = i;
    } else {
      uint const k1 = pickGlobIndividual(randomizer);
      np1 = k1 / plotSize;         
      nt1 = k1 - np1 * plotSize; 
      if( r < pLocalOrImi ) {
	sx = com.get(np1, nt1).speciesIndex;
      } else {
	lastSpecies += 1;
	sx = lastSpecies;
	// log
      }
    }
    com.setSpecies(np, nt, sx);
    // Patch& p = com.get(np, nt);
    // p.speciesIndex = sx;
  }
  com.setTime(curTime);
  return curTime;
}


double
MetaSimulator::advance(double targetTime, TracedMetaCommunity& com,
		       double const minCAtime, double const xcleanEvery = 1)
{
#if !defined(NDEBUG)
  uint const nPlots = com.nPlots;
#endif
  uint const plotSize = com.plotSize;
  uint const nIndividuals = com.nIndividuals();

  uint lastSpecies = com.getLastSpecies();
  
  double const dr = nIndividuals * mu;
  std::exponential_distribution<double> timeToNextDeath(dr);
  std::uniform_int_distribution<int> pickGlobIndividual(0, nIndividuals-1);
  std::uniform_int_distribution<int> pickLocalIndividual(0, plotSize-1);
  std::uniform_real_distribution<double> zeroOne(0.0, 1.0);

  double curTime = com.getTime();
  // at least one clean in inner loop
  //double const cleanEvery = std::min(1.0, targetTime * 0.99); // ???
  double const cleanEvery = std::max(xcleanEvery, minCAtime); // ???
  
  double cleanAt = curTime + cleanEvery;
  // at least one loop to inital target time
  // double curCAtime = -1;

  // look until target time or CA time > minimum CA time
  while( curTime <= targetTime ) {
    {
      double deltaT = timeToNextDeath(randomizer);
      double cpd = curTime + deltaT;
      while( cpd == curTime ) {
	deltaT = timeToNextDeath(randomizer);
	cpd = curTime + deltaT;
      }
      curTime = cpd;
    }
      
    uint const k = pickGlobIndividual(randomizer);
    uint const np = k / plotSize;         assert(0 <= np && np < nPlots);
    uint const nt = k - np * plotSize;    assert(0 <= nt && nt < plotSize );
    double const r = zeroOne(randomizer);

    uint np1,nt1,sx;
    
    if( r < pLocal ) {
      uint const i = pickLocalIndividual(randomizer);
      sx = com.get(np, i).speciesIndex;
      np1 = np;
      nt1 = i;
    } else {
      uint const k1 = pickGlobIndividual(randomizer);
      np1 = k1 / plotSize;         
      nt1 = k1 - np1 * plotSize; 
      if( r < pLocalOrImi ) {
	sx = com.get(np1, nt1).speciesIndex;
      } else {
	lastSpecies += 1;
	sx = lastSpecies;
      }
    }
    com.replace(np, nt, np1, nt1, curTime, sx);

    if( curTime >= cleanAt ) {
      const PatchTrace* const ca = com.cleanUp();
      double const curCAtime = ca->pTime;
      
      if( minCAtime >= 0 && curCAtime > minCAtime ) {
	break;
      }
      
      cleanAt += cleanEvery;
    }
  }

  com.setTime(curTime);
  return curTime;
}

static void
metaDestructor(PyObject* const o)
{
  if( PyCapsule_IsValid(o, "TMC") ) {
    void* c = PyCapsule_GetPointer(o, "TMC");
    delete reinterpret_cast<TracedMetaCommunity*>(c);
  } else if( PyCapsule_IsValid(o, "MC") ) {
    void* c = PyCapsule_GetPointer(o, "MC");
    delete reinterpret_cast<MetaCommunity*>(c);
  }
}

PyObject*
newCommunity(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"nPlots", "plotSize", "startTime", "trace", "metaCommunity", 
				 static_cast<const char*>(0)};
  PyObject* metaCom = 0;
  double startTime = 0;
  uint nPlots = 0, plotSize = 0;
  PyObject* pyTrace = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "|iidOO", const_cast<char**>(kwlist),
				    &nPlots,&plotSize,&startTime,&pyTrace,&metaCom)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  
  bool const trace = (pyTrace && PyObject_IsTrue(pyTrace));
  
  MetaCommunity* com = 0;
  TracedMetaCommunity* tcom = 0;
  
  if( metaCom ) {
    if( PySequence_Check(metaCom) ) {
      int const n = PySequence_Size(metaCom);
      if( n > 0 ) {
	nPlots = n;

	bool err = false;
	for(uint k = 0; k < nPlots; ++k) {
	  PyObject* const pk = PySequence_Fast_GET_ITEM(metaCom, k);
	  if( ! PySequence_Check(pk) ) {
	    err = true;
	    break;
	  }
	  if( k == 0 ) {
	    plotSize = PySequence_Size(pk);
	  } else {
	    if( plotSize != static_cast<uint>(PySequence_Size(pk)) ) {
	      err = true;
	      break;
	    }
	  }
	}
	if( ! err ) {
	  if( trace ) {
	    tcom = new TracedMetaCommunity(nPlots, plotSize, startTime);
	  } else {
	    com = new MetaCommunity(nPlots, plotSize, startTime);
	  }
	  
	  //com = new MetaCommunity(nPlots, plotSize, startTime);
	  for(uint k = 0; k < nPlots; ++k) {
	    PyObject* const pk = PySequence_Fast_GET_ITEM(metaCom, k);
	    for(uint i = 0; i < plotSize; ++i) {
	      PyObject* const pki = PySequence_Fast_GET_ITEM(pk, i);
	      int sx = std::max(PyInt_AsLong(pki),0L);
	      if( com) {
		com->setSpecies(k, i, sx);
	      } else {
		tcom->setSpecies(k, i, sx);
	      }
	    }
	  }
	}
      }
    }

    if( (trace && !tcom) || (!trace && !com) ) {
      PyErr_SetString(PyExc_ValueError, "wrong args: not a valid community");
      return 0;
    }
	
  } else {
    if( nPlots == 0 || plotSize == 0 ) {
      PyErr_SetString(PyExc_ValueError, "wrong args: no community or valid sizes");
      return 0;
    }
    if( trace ) {
      tcom = new TracedMetaCommunity(nPlots, plotSize, startTime);
    } else {
      com = new MetaCommunity(nPlots, plotSize, startTime);
    }
  }

  void* c = trace ? tcom : com;
  PyObject* o = PyCapsule_New(c, trace ? "TMC" : "MC", metaDestructor);
  return o;
}

PyObject*
forwardSim(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"metaCommunity","targetTime", "mu", "lam", "rho", 
				 "minCAtime", "seed", "cleanEvery",
				 static_cast<const char*>(0)};
  PyObject* metaCom = 0;
  double targetTime;
  double minCAtime = -1;
  double cleanEvery = 1;
  double mu = 1, lam = 1, rho = 1;
  double seed = -1;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "Odddd|ddd", const_cast<char**>(kwlist),
				    &metaCom,&targetTime,&mu,&lam,&rho,&minCAtime,&seed,
				    &cleanEvery)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( seed >= 0 ) {
    randomizer.seed(seed);
  }
  
  //bool trace; // = (pyTrace && PyObject_IsTrue(pyTrace));
  // bool onTheFlyCom;
  //MetaCommunity* com = 0;
  //TracedMetaCommunity* tcom = 0;

  MetaSimulator s(mu, lam, rho);
  double endTime;
  PyObject* retVal = 0;
  
  if( PyCapsule_IsValid(metaCom, "TMC") ) {
    TracedMetaCommunity& tcom =
      *reinterpret_cast<TracedMetaCommunity*>(PyCapsule_GetPointer(metaCom, "TMC"));
    //trace = true;
    endTime = s.advance(targetTime, tcom, minCAtime, cleanEvery);
    retVal = tcom.asPyObject();
  } else if( PyCapsule_IsValid(metaCom, "MC") ) {
    MetaCommunity& com =
      *reinterpret_cast<MetaCommunity*>(PyCapsule_GetPointer(metaCom, "MC"));
    // trace = false;
    endTime = s.advance(targetTime, com);
    retVal = com.asPyObject();
  } else {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a valid community") ;
    return 0;
  }

  PyObject* o = PyTuple_New(2);
  PyTuple_SET_ITEM(o, 0, PyFloat_FromDouble(endTime));
  PyTuple_SET_ITEM(o, 1, retVal);
  
  return o;
}

PyObject*
CAcounts(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"metaCommunity", "nwithin", "nbetween",
				 static_cast<const char*>(0)};
  PyObject* metaCom = 0;
  int nwithin = 10, nbetween = 10;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "O|ii", const_cast<char**>(kwlist),
				    &metaCom,&nwithin,&nbetween)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! PyCapsule_IsValid(metaCom, "TMC") ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a valid community") ;
    return 0;
  }

  TracedMetaCommunity& tcom =
    *reinterpret_cast<TracedMetaCommunity*>(PyCapsule_GetPointer(metaCom, "TMC"));

  std::uniform_int_distribution<int> pickPlot(0, tcom.nPlots-1);
  std::uniform_int_distribution<int> pickLocalIndividual(0, tcom.plotSize-1);

  PyObject* tup0 = PyTuple_New(nwithin);
  
  for(int k = 0; k < nwithin; ++k) {
    uint const np = pickPlot(randomizer);
    uint const j = pickLocalIndividual(randomizer);
    uint count = 0;
    for(uint i = 0; i < tcom.plotSize; ++i) {
      if( i != j ) {
	const PatchTrace* const p = tcom.ca(np, i, np, j);
	count += p->pTime > 0;
      }
    }
    PyTuple_SET_ITEM(tup0, k, PyInt_FromLong(count));
  }

  PyObject* tup1 = PyTuple_New(nbetween);
  for(int k = 0; k < nbetween; ++k) {
    uint const np0 = pickPlot(randomizer);
    uint np1 = np0;
    while( np1 == np0 ) {
      np1 = pickPlot(randomizer);
    }
    uint const j = pickLocalIndividual(randomizer);
    uint count = 0;
    for(uint i = 0; i < tcom.plotSize; ++i) {
      const PatchTrace* const p = tcom.ca(np0, j, np1, i);
      count += p->pTime > 0;
    }
    PyTuple_SET_ITEM(tup1, k, PyInt_FromLong(count));
  }
  
  return Py_BuildValue("NN", tup0, tup1);
}

PyObject*
optTraces(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"metaCommunity",
				 static_cast<const char*>(0)};
  PyObject* metaCom = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "O", const_cast<char**>(kwlist),
				    &metaCom)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! PyCapsule_IsValid(metaCom, "TMC") ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a valid community") ;
    return 0;
  }

  TracedMetaCommunity& tcom =
    *reinterpret_cast<TracedMetaCommunity*>(PyCapsule_GetPointer(metaCom, "TMC"));

  tcom.optTraces();
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef neutralsimMethods[] = {
  {"forwardSimulation",	(PyCFunction)forwardSim, METH_VARARGS|METH_KEYWORDS,
   ""},
  {"newCommunity",  	(PyCFunction)newCommunity, METH_VARARGS|METH_KEYWORDS,
   ""},
  {"CAcounts",		(PyCFunction)CAcounts, METH_VARARGS|METH_KEYWORDS,
   ""},
  {"optTraces",		(PyCFunction)optTraces, METH_VARARGS|METH_KEYWORDS,
   ""},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initneutralsim(void)
{
  randomizer.seed( time(0) );
  
  //import_array();
  
  /*PyObject* m = */ Py_InitModule("neutralsim", neutralsimMethods);
}
