// This file is part of biopy.
// Copyright (C) 2013 Joseph Heled
// Author: Joseph Heled <jheled@gmail.com>
// See the files gpl.txt and lgpl.txt for copying conditions.


#include <Python.h>
// keeps asserts
#undef NDEBUG
#include <cassert>

#include <limits>
#include<memory>
#include <algorithm>
#include <vector>
using std::vector;

enum ComparisonResult {
  Default = 0,
  JCcorrection = -1,
  IDENTITY = -2,
  DIVERGENCE = -3,
  STATS = -4,
};

#include "readseq.h"
// typedef unsigned char byte;
// byte const anynuc = 4;
// byte const gap = 5;


static inline double 
stats2distance(uint const              matches,
	       uint const              misMatches,
	       uint const              gaps,
	       ComparisonResult const  what)
{
  assert ( what == IDENTITY || what == JCcorrection || what == DIVERGENCE ) ;
  
  int const seqlen = (matches + misMatches + gaps);

  double seqIdent = matches == 0 ? 0.0 : static_cast<double>(matches)/seqlen;
  if( what == IDENTITY ) {
    return seqIdent;
  }
  
  double dis = 1 - seqIdent;
      
  if( what == JCcorrection ) {
    if( seqIdent <= 1/4. ) {
      if( seqlen == 0 ) {
	return 10;
      }
      seqIdent = (int(1+seqlen/4.)/float(seqlen) + .25)/2;
    }
    dis = fabs(-3./4 * log((seqIdent * 4 - 1)/3));
  }
  return dis;
}
      
/**
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
	  s[nsite] = anynuc; break;
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
**/

#if 0
class Choose2 {
private:
  static long* choose2tab;
  static uint ns;
public:
  static void init(uint n) {
    ns = n;
    delete [] choose2tab;
    choose2tab = new long[n];
    for(uint i = 0; i < n; ++i) {
      choose2tab[i] = (i*(i-1)) / 2;
    }
  }

  static inline ulong c2(uint i) {
    return i < ns ? choose2tab[i] : (long(i)*(i-1)) / 2;
  }
};

long* Choose2::choose2tab = 0;
uint Choose2::ns = 0;
#endif

struct ProfileCounts {
public:
  int matches;
  int mis;
  int gaps;

  inline ProfileCounts(void) {}

  inline
  ProfileCounts(const int* p) :
    matches(0),
    mis(0)
    {
    int const tn = (p[0]+p[1]+p[2]+p[3]);
    for(int k = 0; k < 4; ++k) {
      int const pk = p[k];
      matches += pk * (pk-1);
      mis += pk * (tn - pk);
    }
    matches /= 2;
    // every mismatch counted twice
    mis /= 2;
    int const an = p[anynuc];
    matches += tn * an + (an*(an-1))/2;
    gaps = (tn+an) * p[gap];
  }
};

template<typename T>
class MatchScoreValues {
private:
  // initalization hackish trick
  PyObject*	values;
  
public:
  bool valid(void) const {
    return values != reinterpret_cast<PyObject*>(0x1);
  }
  
  MatchScoreValues(T mt, T mis, T gapp, bool free) :
    matchScore(mt),
    misMatchScore(mis),
    gapPenalty(gapp),
    freeEndGaps(free)
    {}

  MatchScoreValues(PyObject* values);

  T const matchScore;
  T const misMatchScore;
  T const gapPenalty;
  bool const freeEndGaps;

  inline T
  scoreMatching(byte const n1, byte const n2) const {
    return ((n1 == n2 || n1 == anynuc || n2 == anynuc) ? matchScore : misMatchScore);
  }

#if defined(ALLASSERTS)
  // complexity slows down .... doubly so since profiles tend to get long
  inline T
  scoreMatching(const int* const p1, uint const n1, const int* const p2, uint const n2) const {
    uint mt = 0;
    for(uint i = 0; i < 4; ++i) {
      uint const x = (p1[i] + p2[i]);
      mt += (x*(x-1));      
    }
    mt /= 2;
    
    int const p1g = p1[gap];
    int const p2g = p2[gap];
    uint const p1NonGap = n1 - p1g;
    uint const p2NonGap = n2 - p2g;
    uint const nx = p1[anynuc] + p2[anynuc];
    mt += (p1NonGap + p2NonGap - nx) * nx +
      (p1[anynuc]*(p1[anynuc]-1))/2 + (p2[anynuc]*(p2[anynuc]-1))/2 +
      p1[anynuc] * p2[anynuc];
    
    uint const gaps = (p1g+p2g) * p2NonGap + (p2g+p1g) * p1NonGap;
    
    uint const totNucs = n1 + n2;
    
    uint const gg = ((p1g*(p1g-1)) + (p2g*(p2g-1)))/2 + p1g*p2g;
    int const ms = (totNucs*(totNucs-1))/2 - (mt + gaps + gg);
    assert (ms >= 0);
    
    return mt * matchScore + ms * misMatchScore + gaps * gapPenalty;
  }

  // cached, so complexity not a big problem
  inline T
  scoreGap(const int* const p, uint n, uint np) const {
#if defined(ALLASSERTS)
    assert (np == uint(p[0] + p[1] + p[2] + p[3] + p[4] + p[5]));
#endif
    
    uint const nNotGap = np - p[gap];
    uint const pNucs = nNotGap - p[anynuc];
    
    uint mt = pNucs * p[anynuc];
    for(uint i = 0; i < 4; ++i) {
      mt += (p[i]*(p[i]-1))/2;
    }
    uint const totNucs = np + n;
    uint const gg = (n*(n-1))/2 + (p[gap]*(p[gap]-1))/2;

    uint const g = nNotGap * n;
    uint const ms = (totNucs*(totNucs-1))/2 - (mt + g + gg);
    return mt * matchScore + ms * misMatchScore + g * gapPenalty;
  }
#endif
  
  // gap penalty across everything
  inline T
  scoreGap(ProfileCounts const& pc, uint n, uint pn) const {
    uint const mt = pc.matches;
    uint const ms = pc.mis;
    uint const gp = pc.gaps + n * pn;
    return mt * matchScore + ms * misMatchScore + gp * gapPenalty;
  }

  // gaps across profiles count as a match. maybe they need another match score?
  inline T
  scoreMatching(const int* const p2j,
		ProfileCounts const& ri, ProfileCounts const& cj,
		const ProfileCounts* const iCounts) const {
    int mt = ri.matches + cj.matches;
    int ms = ri.mis     + cj.mis;
    int gp = ri.gaps    + cj.gaps;

    for(int k = 0; k < 4; ++k) {
      if( int const n = p2j[k] ) {
	auto const& ck = iCounts[k];
	mt += n * ck.matches;
	ms += n * ck.mis;
	//gp += n * ck.gaps;
	mt += n * ck.gaps;
      }
    }
    {
      if( int const n = p2j[gap] ) {
	//gp += n * iCounts[gap].gaps;
	mt += n * iCounts[gap].gaps;	
      }
    }
    {
      if( int const n = p2j[anynuc] ) {
	auto const& an = iCounts[anynuc];
	mt += n * an.matches;
	//gp += n * an.gaps;
	mt += n * an.gaps;	
      }
    }      

    T const sm = matchScore * mt + misMatchScore * ms + gp * gapPenalty;
    return sm;
  }
};


template<typename T>
inline T
getValue(PyObject* o, T dflt)
{
  return (o && o != Py_None) ? PyFloat_AsDouble(o) : dflt;
}

template<int>
inline int
getValue(PyObject* o, int dflt)
{
  return (o && o != Py_None) ? PyInt_AsLong(o) : dflt;
}

template<bool>
inline bool
getValue(PyObject* o, bool dflt)
{
  return (o && o != Py_None) ? PyObject_IsTrue(o) : dflt;
}

template<typename T>
inline T
lgetValue(PyObject* s, int i, T dflt)
{
  if( s ) {
    PyObject* o = PySequence_Fast_GET_ITEM(s,i);
    return getValue<T>(o, dflt);
  }
  return dflt;
}

template<typename T>
MatchScoreValues<T>::MatchScoreValues(PyObject* _values) :
  values( (_values != Py_None &&
	   _values && (PySequence_Check(_values) && PySequence_Size(_values) == 4) ) ?
	  _values : 0 ),
  matchScore(lgetValue<T>(values, 0, 10)),
  misMatchScore(lgetValue<T>(values, 1, -5)),
  gapPenalty(lgetValue<T>(values, 2, -6)) ,
  freeEndGaps(lgetValue<bool>(values, 3, true))
{
  // Ugly hack to flag invalid. not C++ strongest point. I could throw and
  // catch, but then the catch need to extend over the lifetime of the variable.
  if( _values && _values != Py_None && ! values ) {
    values = reinterpret_cast<PyObject*>(0x1);
  }
}
  

template<typename T>
class Alignment {
public:
    Alignment(const byte* const  _s1,
	      uint const         _lseq1,
	      const byte* const  _s2,
	      uint const         _lseq2,
	      Alignment*         _prev = 0) :
    s1(_s1),
    lseq1(_lseq1),
    s2(_s2),
    lseq2(_lseq2),
    sz((lseq1+1) * (lseq2+1)),
    score(new float[sz]),
    prev(_prev)
    {}
 
  ~Alignment() {
    delete [] score;
  }
  
  const byte* const s1;
  uint const        lseq1;
  const byte* const s2;
  uint const        lseq2;

  void	getStats(MatchScoreValues<T> const&  scores,
		 int&                        matches,
		 int&                        misMatches,
		 int&                        gaps);

  uint  align(MatchScoreValues<T> const&  scores,
	      byte* const                 alignment);

private:
  uint initScores(T gapPenalty);
  
  void fillScoreTable(MatchScoreValues<T> const& mScores);
  
  void fillScoreTable(MatchScoreValues<T> const&  mScores,
		      Alignment const&            prev);

  void fillScores(MatchScoreValues<T> const& mScores);

  uint const  sz;
  T*   const  score;

  const Alignment*	prev;
};

template<typename T>
inline void
Alignment<T>::fillScores(MatchScoreValues<T> const& mScores) {
  initScores(mScores.freeEndGaps ? 0 : mScores.gapPenalty);
  if( prev ) {
    fillScoreTable(mScores, *prev);
  } else {
    fillScoreTable(mScores);
  }
}

template<typename T>
uint
Alignment<T>::initScores(T const gapPenalty)
{
  std::fill(score, score+sz, 0);

  if( gapPenalty != 0 ) {
    T penalty = gapPenalty;
    for(T* s = score + 1; s < score + lseq2+1; ++s) {
#if defined(ALLASSERTS)
      assert ( s < score + sz );
#endif
      *s = penalty;
      penalty += gapPenalty;
    }
  
    penalty = gapPenalty;
    for(T* s = score + (lseq2+1); s < score + sz; s += lseq2+1) {
#if defined(ALLASSERTS)
      assert ( s < score + sz );
#endif
      *s = penalty;
      penalty += gapPenalty;
    }
  }
  return sz;
}

template<typename T>
void
Alignment<T>::fillScoreTable(MatchScoreValues<T> const& mScores)
{
  float* m1m1 = score;
  for(uint i = 1; i <= lseq1; ++i) {
    byte const s1i = s1[i-1];
    for(const byte* s2j = s2; s2j < s2 + lseq2; ++s2j) {
      T const match = *m1m1 + mScores.scoreMatching(s1i, *s2j);      
      T const del = *(m1m1+1) + mScores.gapPenalty;
      T const ins = *(m1m1+1+lseq2) + mScores.gapPenalty;
#if defined(ALLASSERTS)
      assert ( uint((m1m1+lseq2+2) - score) < ((lseq1+1) * (lseq2+1)) );
#endif
      
      *(m1m1+lseq2+2) = std::max(std::max(match, del), ins);
      m1m1 += 1;
    }
    m1m1 += 1;
  }
}

template<typename T>
void
Alignment<T>::fillScoreTable(MatchScoreValues<T> const&  mScores,
			     Alignment const&            prev)
{
  // assume same match scores (caller responsibility)
  int s1MatchLen = 0;
  int const s1mx = std::min(lseq1,prev.lseq1);
  for(; s1MatchLen < s1mx; ++s1MatchLen) {
    if( s1[s1MatchLen] != prev.s1[s1MatchLen] ) {
      break;
    }
  }

  if( s1MatchLen ) {
    int s2MatchLen = 0;
    int const s2mx = std::min(lseq2,prev.lseq2);
    for(; s2MatchLen < s2mx; ++s2MatchLen) {
      if( s2[s2MatchLen] != prev.s2[s2MatchLen] ) {
	break;
      }
    }
  
    for(int i = 1; i < s1MatchLen+1; ++i) {
      const T* p = prev.score + i*(prev.lseq2+1) + 1;
      //float* m = score + off;
      uint const off = i*(lseq2+1) + 1;
      T* m = score + off;
      T* end = m + s2MatchLen;
      for(/**/; m < end; ++m, ++p) {
	*m = *p;
      }
      if( uint(s2MatchLen) < lseq2 ) {
	m -= (lseq2+2);
	byte const s1i = s1[i-1];
	for(const byte* s2j = s2+s2MatchLen; s2j < s2 + lseq2; ++s2j) {
	  T const match = *m + mScores.scoreMatching(s1i, *s2j);      
	  T const del = *(m+1) + mScores.gapPenalty;
	  T const ins = *(m+1+lseq2) + mScores.gapPenalty;
#if defined(ALLASSERTS)
	  assert ( uint((m+lseq2+2) - score) < ((lseq1+1) * (lseq2+1)) );
#endif
      
	  *(m+lseq2+2) = std::max(std::max(match, del), ins);
	  m += 1;
	}
      }
    }
  }
  
  T* m1m1 = score + s1MatchLen*(lseq2+1);
  for(uint i = s1MatchLen+1; i <= lseq1; ++i) {
    byte const s1i = s1[i-1];
    for(const byte* s2j = s2; s2j < s2 + lseq2; ++s2j) {
      T const match = *m1m1 + mScores.scoreMatching(s1i, *s2j);      
      T const del = *(m1m1+1) + mScores.gapPenalty;
      T const ins = *(m1m1+1+lseq2) + mScores.gapPenalty;
#if defined(ALLASSERTS)
      assert ( uint((m1m1+lseq2+2) - score) < ((lseq1+1) * (lseq2+1)) );
#endif
      
      *(m1m1+lseq2+2) = std::max(std::max(match, del), ins);
      m1m1 += 1;
    }
    m1m1 += 1;
  }
}

// alignment: at least [2*(lseq1 + lseq2)], returned reversed

template<typename T>
uint
Alignment<T>::align(MatchScoreValues<T> const&  scores,
		    byte* const                 alignment)
{
  byte* al0 = alignment;
  byte* al1 = alignment + (lseq1 + lseq2);

  int iRow = lseq1;
  int jCol = lseq2;

  fillScores(scores);
  if( scores.freeEndGaps ) {
    uint iMaxRow = sz-lseq2-1;
    T mxLastRow = score[iMaxRow];
    for(uint l = iMaxRow+1; l < sz; ++l) {
      if( score[l] >= mxLastRow ) {
	mxLastRow = score[l];
	iMaxRow = l;
      }
    }

    uint iMaxCol = lseq2;
    T mxLastCol = score[iMaxCol];
    for(uint l = iMaxCol+lseq2+1; l < sz; l += lseq2+1) {
      if( score[l] >= mxLastCol ) {
	mxLastCol = score[l];
	iMaxCol = l;
      }
    }

    if( mxLastCol > mxLastRow ) {
      int const im = (iMaxCol - lseq2)/(lseq2+1);
      while( iRow > im ) {
	*al0 = s1[iRow-1];
	*al1 = gap;
	al0 += 1;
	al1 += 1;
	iRow -= 1;
      }
    } else {
      int const jm = iMaxRow - (sz-lseq2-1);
      while( jCol > jm ) {
	*al0 = gap;
	*al1 = s2[jCol-1];
	al0 += 1;
	al1 += 1;
	jCol -= 1;
      }
    }
  }
  
  while( iRow > 0 and jCol > 0 ) {
    uint const cur = iRow*(lseq2+1) + jCol;
    T const score_current = score[cur];
    T const score_diagonal = score[cur - lseq2 - 2];

    if( score_current == score_diagonal +
	scores.scoreMatching(s1[iRow-1], s2[jCol-1]) ) {
      *al0 = s1[iRow-1]; 
      *al1 = s2[jCol-1];
      iRow -= 1;
      jCol -= 1;
    } else {
      T const score_up = score[cur - lseq2 - 1];

      if( score_current == score_up + scores.gapPenalty ) {
	*al0 = s1[iRow-1];
	*al1 = gap;
	iRow -= 1;
      } else {
#if !defined(NDEBUG)
	T const score_left = score[cur - 1];
#endif
	assert ( score_current == score_left + scores.gapPenalty ) ;
	
	*al0 = gap;
	*al1 = s2[jCol-1];
	jCol -= 1;
      }
    }
    al0 += 1;
    al1 += 1;
  }
  
  while( iRow > 0 ) {
    *al0 = s1[iRow-1];
    *al1 = gap;
    al0 += 1;
    al1 += 1;
    iRow -= 1;
  }
  
  while( jCol > 0 ) {
    *al0 = gap;
    *al1 = s2[jCol-1];
    al0 += 1;
    al1 += 1;
    jCol -= 1;
  }
  
  return al0 - alignment;
}

template<typename T>
void
Alignment<T>::getStats(MatchScoreValues<T> const&  scores,
		       int&                        matches,
		       int&                        misMatches,
		       int&                        gaps)
{
  matches = misMatches = gaps = 0;

  int iRow = lseq1;
  int jCol = lseq2;
  
  fillScores(scores);
  if( scores.freeEndGaps ) {
    uint iMaxRow = sz-lseq2-1;
    T mxLastRow = score[iMaxRow];
    for(uint l = iMaxRow+1; l < sz; ++l) {
      if( score[l] >= mxLastRow ) {
	mxLastRow = score[l];
	iMaxRow = l;
      }
    }

    uint iMaxCol = lseq2;
    T mxLastCol = score[iMaxCol];
    for(uint l = iMaxCol+lseq2+1; l < sz; l += lseq2+1) {
      if( score[l] >= mxLastCol ) {
	mxLastCol = score[l];
	iMaxCol = l;
      }
    }

    if( mxLastCol > mxLastRow ) {
      iRow = (iMaxCol - lseq2)/(lseq2+1);
    } else {
      jCol = iMaxRow - (sz-lseq2-1);
    }
  } 
  
  while( iRow > 0 and jCol > 0 ) { 
    uint const cur = iRow*(lseq2+1) + jCol;
    T const score_current = score[cur];
    T const score_diagonal = score[cur - lseq2 - 2];

    if( score_current == score_diagonal + scores.scoreMatching(s1[iRow-1], s2[jCol-1]) ) {
      if( s1[iRow-1] == s2[jCol-1] ) {
	matches += 1;
      } else {
	misMatches += 1;
      }
      iRow -= 1;
      jCol -= 1;
    } else {
      gaps += 1;
      T const score_left = score[cur - 1];
      
      if ( score_current == score_left + scores.gapPenalty ) {
	jCol -= 1;
      } else {
#if !defined(NDEBUG)
	T score_up = score[cur - lseq2 - 1];
#endif
	assert (score_current == score_up + scores.gapPenalty);
	
	iRow -= 1;
      }
    }
  }

  if( !scores.freeEndGaps ) {
    gaps += iRow + jCol;
  }
}

PyObject*
globAlign(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char* kwlist[] = {"seq0", "seq1", "report", "strip", "scores",
				 static_cast<const char*>(0)};
  PyObject* pseq1 = 0;
  PyObject* pseq2 = 0;

  ComparisonResult resultType = Default;
  PyObject* pStrip = 0;
  PyObject* mScores = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "OO|iOO", const_cast<char**>(kwlist),
				    &pseq1,&pseq2,&resultType,&pStrip,&mScores)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PySequence_Check(pseq1) && PySequence_Check(pseq2)) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not sequences") ;
    return 0;
  }

  MatchScoreValues<float> const scores(mScores);
  if( ! scores.valid() ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: invalid scores") ;
    return 0;
  }
  
  uint lseq1;
  uint lseq2;
  // always strip for alignment
  byte* s1 = readSequence(pseq1, lseq1, true);
  byte* s2 = readSequence(pseq2, lseq2, true);
  if( ! (s1 && s2 ) ) {
    return 0;
  }

  Alignment<float> al(s1, lseq1, s2, lseq2);
  
  PyObject* ret = 0;
  
  uint const ltot = lseq1 + lseq2;
  if( resultType != Default ) {
    int matches, misMatches, gaps;
    al.getStats(scores, matches, misMatches, gaps);
    
    if( resultType == STATS ) {
	ret = PyTuple_New(3);
	PyTuple_SET_ITEM(ret, 0, PyLong_FromLong(matches));
	PyTuple_SET_ITEM(ret, 1, PyLong_FromLong(misMatches));
	PyTuple_SET_ITEM(ret, 2, PyLong_FromLong(gaps));
    } else {
      double const dis = stats2distance(matches, misMatches, gaps, resultType);
      ret = PyFloat_FromDouble(dis);
    }
  } else {
    byte* const alignment = new byte [2*ltot];
    
    uint const alen = al.align(scores, alignment);

    PyObject* al = PyTuple_New(2);
    PyObject* ps0 = PyTuple_New(alen);
    const byte* ps = alignment + alen-1;
    for(uint i = 0; i < alen; ++i,--ps) {
      PyTuple_SET_ITEM(ps0, i, PyInt_FromLong(*ps));
    }
    PyTuple_SET_ITEM(al, 0, ps0);
  
    PyObject* ps1 = PyTuple_New(alen);
    ps = alignment + ltot + alen-1;
    for(uint i = 0; i < alen; ++i,--ps) {
      PyTuple_SET_ITEM(ps1, i, PyInt_FromLong(*ps));
    }
    PyTuple_SET_ITEM(al, 1, ps1);

    delete [] alignment;

    ret = al;
  }
  
  delete [] s1;
  delete [] s2;

  return ret;
}


class SeqsList {
public:
  SeqsList(uint n, const uint* sl, const byte** sq) :
    nSeqs(n),
    seqslen(sl),
    seqs(sq)
    {}
  
  ~SeqsList() {
    for(uint j = 0; j < nSeqs; ++j) {
      delete [] seqs[j];
    }
    
    delete [] seqslen;
    delete [] seqs;
  }

  int aligned(void) const {
    if( nSeqs == 0 ) {
      return -1;
    }
    uint n = seqslen[0];
    for(uint j = 1; j < nSeqs; ++j) {
      if( seqslen[j] != n ) {
	return -1;
      }
    }
    return n;
  }
  
  uint const nSeqs;
  const uint* const seqslen;
  const byte** seqs ;
};

static bool inline
isSequence(PyObject* const s)
{
  // only check the type of the first element. can be fooled.
  return (PyString_Check(s) ||
	  ( PySequence_Check(s) && PySequence_Size(s) > 0 &&
	    PyInt_Check(PySequence_Fast_GET_ITEM(s,0))) );
}

// Read in a sequence of one or more DNA sequences. With 'strip', remove
// gaps. Without, check that all sequences have equal length.
// should be a constructor!!

SeqsList*
readSeqsIn(PyObject* const pseqs, bool const strip)
{
  bool const single = isSequence(pseqs);
  
  int const nseqs = single ? 1 : PySequence_Size(pseqs);
  uint* seqslen = new uint [nseqs];
  const byte** seqs = new const byte* [nseqs];

  if( single ) {
    uint nseq;
    const byte* s = readSequence(pseqs, nseq, strip);
    
    if( ! s ) {
      PyErr_Format(PyExc_ValueError, "wrong args: not a sequence") ;
      delete [] seqslen;
      delete [] seqs;
      return 0;
    }
    seqslen[0] = nseq;
    seqs[0] = s;
  } else {
    for(int j = 0; j < nseqs; ++j) {
      PyObject* seq = PySequence_Fast_GET_ITEM(pseqs, j);

      uint nseq;
      const byte* s = readSequence(seq, nseq, strip);
    
      if( ! s ) {
	PyErr_Format(PyExc_ValueError, "wrong args: seqs[%d] not a sequence", j) ;
	delete [] seqslen;
	delete [] seqs;
	return 0;
      }

      seqslen[j] = nseq;
    
      if( !strip && j > 0 && nseq != seqslen[j-1] ) {
	delete [] seqslen;
	delete [] seqs;
	PyErr_SetString(PyExc_ValueError, "Sequences not aligned") ;
	return 0;
      }
      
      seqs[j] = s;
    }
  }
  return new SeqsList(nseqs, seqslen, seqs);
}

template<typename T>
PyObject*
distmat(PyObject*         pseqs,
	PyObject*         palign,
	ComparisonResult  resultType,
	PyObject*         pReorder,
	PyObject*         mScores,
	T*                retSpace)
{
  MatchScoreValues<float> const scores(mScores);
  if( ! scores.valid() ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: invalid scores") ;
    return 0;
  }

  bool const align = (palign == 0 || PyObject_IsTrue(palign));
  bool const strip = align;

  std::unique_ptr<const SeqsList> sq1(readSeqsIn(pseqs, strip));
  if( ! sq1 || sq1->nSeqs == 0 ) {
    return 0;
  }
  
  uint const nseqs = sq1->nSeqs;
  int* const order = (pReorder && PySequence_Check(pReorder)) ? new int[nseqs] : 0;  

  if( order ) {
    uint const sz = PySequence_Size(pReorder);
    if( sz != nseqs ) {
      delete [] order;
      return 0;
    }
    vector<int> check(sz, -2); //  = new order(n);

    for(uint n = 0; n < sz; ++n) {
      PyObject* pn = PySequence_Fast_GET_ITEM(pReorder, n);
      uint v = PyInt_AsLong(pn);
      if( ! (0 <= v && v < sz && check[v] == -2) ) {
	delete [] order;
	return 0;
      }
      check[v] = n;
      order[n] = v;
    }
  }

  uint const nResults = (nseqs * (nseqs-1)) / 2;
  PyObject* res = !retSpace ? PyTuple_New(nResults) : 0;
  int nr = 0;
 
  int matches, misMatches, gaps;

  Alignment<float>* prev = 0;
  
  for(uint j = 0; j < nseqs-1; ++j) {
    for(uint i = j+1; i < nseqs; ++i) {
      if( align ) {
	if( order ) {
	  Alignment<float>* al =
	    new Alignment<float>(sq1->seqs[i], sq1->seqslen[i],
				 sq1->seqs[j], sq1->seqslen[j], prev);
	
	  al->getStats(scores, matches, misMatches, gaps);
	  delete prev; prev = al;
	} else {
	  Alignment<float>  al(sq1->seqs[i], sq1->seqslen[i],
			       sq1->seqs[j], sq1->seqslen[j]);
	
	  al.getStats(scores, matches, misMatches, gaps);
	}
      } else {
	const byte* s1 = sq1->seqs[i];
	const byte* s2 = sq1->seqs[j];
	matches = misMatches = gaps = 0;

	int k = 0;
	int send = sq1->seqslen[i];
	if( scores.freeEndGaps ) {
	  while( (s1[k] == gap || s2[k] == gap) && k < send ) {
	    ++k;
	  }
	  while( k < send && (s1[send-1] == gap || s2[send-1] == gap) ) {
	    --send;
	  }
	}
	for(/**/; k < send; ++k) {
	  byte c1 = s1[k];
	  byte c2 = s2[k];
	  if( c1 != gap && c2 != gap ) {
	    if( c1 == c2 ) {
	      matches += 1;
	    } else {
	      misMatches += 1;
	    }
	  } else if( ! (c1 == gap && c2 == gap) ) {
	    gaps += 1;
	  }
	}
      }

      double const dis = stats2distance(matches, misMatches, gaps, resultType);

      if( order ) {
	int const o1 = order[i];
	int const o2 = order[j];
	
	int const mij = std::min(o1,o2);
	int const xij = std::max(o1,o2);
	uint const pos = mij*(2*nseqs - 1 - mij)/2 + (xij-mij-1);
	assert( 0 <= pos && pos < nResults );
	
	if( retSpace ) {
	  retSpace[pos] = dis;
	} else {
	  PyTuple_SET_ITEM(res, pos, PyFloat_FromDouble(dis));
	}
      } else {
	if( retSpace ) {
	  *retSpace = dis;
	  retSpace += 1;
	} else {
	  PyTuple_SET_ITEM(res, nr, PyFloat_FromDouble(dis));
	  nr += 1;
	}
      }
    }
  }
  delete prev;

  delete [] order;
  return retSpace ? Py_None : res;
}

PyObject*
distMat(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"seqs", "align", "report", "reorder", "scores",
				 static_cast<const char*>(0)};
  PyObject* pseqs = 0;
  PyObject* palign = 0;
  PyObject* mScores = 0;
  ComparisonResult resultType = DIVERGENCE;
  PyObject* pReorder = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "O|OiOO", const_cast<char**>(kwlist),
				    &pseqs, &palign, &resultType, &pReorder,&mScores)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  
  if( ! (PyTuple_Check(pseqs) || PyList_Check(pseqs)) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a list of sequences") ;
    return 0;
  }

  if( ! (STATS < resultType && resultType <= JCcorrection) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args (return type)");
    return 0;
  }
    
  PyObject* res = distmat(pseqs, palign, resultType, pReorder, mScores,
			  static_cast<float*>(0));
  return res;
}


PyObject*
distPairs(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"seqs1", "seqs2", "align", "report", "order","scores",
				 static_cast<const char*>(0)};
  PyObject* pseqs1 = 0;
  PyObject* pseqs2 = 0;
  PyObject* palign = 0;
  ComparisonResult resultType = DIVERGENCE;
  PyObject* mScores = 0;
  PyObject* order = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "OO|OiOO", const_cast<char**>(kwlist),
				    &pseqs1, &pseqs2, &palign, &resultType,
				    &order,&mScores)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  
  if( ! ((PyTuple_Check(pseqs1) || PyList_Check(pseqs1) || PyString_Check(pseqs1)) &&
	 (PyTuple_Check(pseqs2) || PyList_Check(pseqs2) || PyString_Check(pseqs2))) )  {
    PyErr_SetString(PyExc_ValueError, "wrong args: not lists of sequences") ;
    return 0;
  }

  if( ! (STATS < resultType && resultType <= JCcorrection) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args (return type)");
    return 0;
  }

  MatchScoreValues<float> const scores(mScores);
  if( ! scores.valid() ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: invalid scores") ;
    return 0;
  }
  
  bool const align = (palign == 0 || PyObject_IsTrue(palign));

  std::unique_ptr<const SeqsList> sq1(readSeqsIn(pseqs1, align));
  std::unique_ptr<const SeqsList> sq2(readSeqsIn(pseqs2, align));
  
  if( ! sq1 || ! sq2 ) {
    return 0;
  }

  if( ! align ) {
    int const n = sq1->aligned();
    if( n == -1 || n != sq2->aligned() ) {
      PyErr_SetString(PyExc_ValueError, "Sequences not aligned") ;
      return 0;
    }
  }

  uint* orders[2] = {0,0};
  if( align && order ) {
    if( ! (PySequence_Check(order) && PySequence_Size(order) == 2) ) {
      return 0;
    }
    // order should be a permutation (not validated)
    for(int i = 0; i < 2; ++i) {
      PyObject* o1 = PySequence_Fast_GET_ITEM(order,i);
      if( o1 == Py_None ) {
	continue;
      }
      uint const n = i == 0 ? sq1->nSeqs : sq2->nSeqs;
      if( ! (PySequence_Check(o1) && PySequence_Size(o1) == int(n)) ) {
	return 0;
      }
      orders[i] = new uint [n];
      for(uint k = 0; k < n; ++k) {
	uint ok = PyInt_AsLong(PySequence_Fast_GET_ITEM(o1,k));
	if( ! (0 <= ok && ok < n) ) {
	  return 0;
	}
	orders[i][k] = ok;
      }
    }
  }
  
  bool const returnFlat = isSequence(pseqs1);
  
  PyObject* res = PyTuple_New( returnFlat ? sq2->nSeqs : sq1->nSeqs );

  Alignment<float>* prev = 0;
  
  int matches, misMatches, gaps;
  for(uint j = 0; j < sq1->nSeqs; ++j) {
    uint const rj = orders[0] ? orders[0][j] : j;
    
    PyObject* resj = returnFlat ?  0 : PyTuple_New( sq2->nSeqs );
    if( resj ) {
      PyTuple_SET_ITEM(res, rj, resj);
    }
    
    for(uint i = 0; i < sq2->nSeqs; ++i) {
      uint const ri = orders[1] ? orders[1][i] : i;
      if( align ) {
	Alignment<float>* al =
	  new Alignment<float> (sq2->seqs[ri], sq2->seqslen[ri],
				 sq1->seqs[rj], sq1->seqslen[rj], prev);
	al->getStats(scores, matches, misMatches, gaps);
	delete prev; prev = al;
      } else {
	assert(ri == i && rj == j);
	
	const byte* s1 = sq2->seqs[i];
	const byte* s2 = sq1->seqs[j];
	matches = misMatches = gaps = 0;

	int k = 0;
	int send = sq1->seqslen[j];
	if( scores.freeEndGaps ) {
	  while( (s1[k] == gap || s2[k] == gap) && k < send ) {
	    ++k;
	  }
	  while( k < send && (s1[send-1] == gap || s2[send-1] == gap) ) {
	    --send;
	  }
	}
	for(/**/; k < send; ++k) {
	  byte const c1 = s1[k];
	  byte const c2 = s2[k];
	  if( c1 != gap && c2 != gap ) {
	    if( c1 == c2 ) {
	      matches += 1;
	    } else {
	      misMatches += 1;
	    }
	  } else if( ! (c1 == gap && c2 == gap) ) {
	    gaps += 1;
	  }
	}
      }

      double const dis = stats2distance(matches, misMatches, gaps, resultType);
      PyTuple_SET_ITEM(returnFlat ? res : resj, ri, PyFloat_FromDouble(dis));
    }
  }

  delete prev;
  for(int i = 0; i < 2; ++i) {
    delete [] orders[i];
  }

  return res;
}

static inline float
scoreSiteGap(int const                 nSite,
	     const int* const*  const  prof,
	     uint  const               nSeqs,
	     float const               gapPenalty)
{
  return gapPenalty * (nSeqs - prof[nSite-1][gap]);
}

static inline float
scoreSite(int const           nSite,
	  int const           nuc,
	  const int* const*  const prof,
	  uint  const         nSeqs,
	  float const         matchScore,
	  float const         misMatchScore,
	  float const         gapPenalty)
{
  const int* const c = prof[nSite-1];
#if defined(ALLASSERTS)
  assert (nuc != gap);
#endif
  
  return nuc == anynuc ? (gapPenalty * c[gap]) :
    (matchScore * c[nuc] +
     misMatchScore * (c[0]+c[1]+c[2]+c[3] - c[nuc])
     + gapPenalty * c[gap]);
}


// result reversed
template<typename T>
static int
alignToProf(byte*        inseq,
	    uint const 	 seqLen,
	    int**        profile,
	    uint const   ns,
	    T const  matchScore,
	    T const  misMatchScore,
	    T const  gapPenalty)
{
  uint nNucs = 0;
  byte seq[ns];
  for(uint k = 0; k < seqLen; ++k) {
    if( inseq[k] != gap ) {
      seq[nNucs] = inseq[k];
      nNucs += 1;
    }
  }

  if( nNucs > ns ) {
    return -1;
  }
  
  uint const sz = (nNucs+1) * (ns+1);
  T* const rawscr = new T[sz + (ns+1)*6];
  T** const score = new T* [ns+1];
  
  for(uint i = 0; i <= ns; ++i) {
    score[i] = rawscr + i * (nNucs+1);
  }

  T const bla = -123465.887;
  for(uint i = 0; i < sz; ++i) {
    rawscr[i] = bla;
  }

  uint tot = 0;
  for(uint i = 0; i < 6; ++i) {
    tot += profile[0][i];
  }

  // caching site scores is a worthwhile speedup (x3)
  T* const gapScores = rawscr + sz;
  
  T* siteScores[5];
  siteScores[0] = gapScores + ns+1;
  for(uint i = 1; i < 5; ++i) {
    siteScores[i] = siteScores[i-1] + ns + 1;
  }
  
  for(uint ni = 1; ni <= ns; ++ni) {
    for(uint i = 0; i < 5; ++i) {
      siteScores[i][ni] = scoreSite(ni, i, profile, tot, matchScore, misMatchScore, gapPenalty);
    }
  }

  score[0][0] = 0;
  for(uint ni = 1; ni <= ns; ++ni) {
    gapScores[ni] = scoreSiteGap(ni, profile, tot, gapPenalty);
    score[ni][0] = 0;
  }

  // nk + ns - nNucs
  for(uint nk = 1; nk <= nNucs; ++nk) {
#if defined(ALLASSERTS)
    assert( siteScores[seq[nk-1]][nk] ==
	    scoreSite(nk, seq[nk-1], profile, tot, matchScore, misMatchScore, gapPenalty) );
#endif
    assert (score[nk-1][nk-1] != bla);
    
    score[nk][nk] = score[nk-1][nk-1] + siteScores[seq[nk-1]][nk];

    for(uint ni = nk+1; ni <= nk + ns - nNucs; ++ni) {
      // Expensive asserts
      // assert (scoreSiteGap(ni, profile, tot, gapPenalty) == gapScores[ni]);
      // assert (siteScores[seq[nk-1]][ni] == scoreSite(ni, seq[nk-1], profile, tot,
      // 					     matchScore, misMatchScore, gapPenalty));
      assert (score[ni-1][nk] != bla && score[ni-1][nk-1] != bla );

      score[ni][nk] = std::max(score[ni-1][nk]   + gapScores[ni],
			       score[ni-1][nk-1] + siteScores[seq[nk-1]][ni]);
    }
  }

  uint j = nNucs;

  for(uint ni = nNucs+1; ni <= ns; ++ni) {
    if( score[ni][nNucs] > score[j][nNucs] ) {
      j = ni;
    }
  }
  
  uint pos = 0;
  for(/**/; pos < ns - j; ++pos) {
    inseq[pos] = gap;
  }
	
  uint i = nNucs;
  while( j > i ) {
    //assert( i >= 0 && j > 0 );
    // Expensive asserts
    // assert(i == 0 || siteScores[seq[i-1]][j] ==
    // 	   scoreSite(j, seq[i-1], profile, tot, matchScore, misMatchScore, gapPenalty));
    
    if( i > 0 && score[j][i] == score[j-1][i-1] + siteScores[seq[i-1]][j] ) {
      inseq[pos] = seq[i-1];
      i -= 1;
      j -= 1;
    } else {
      inseq[pos] = gap;
      j -= 1;
    }
    
    ++pos;
  }
  
  while( i > 0 ) {
    inseq[pos++] = seq[i-1];
    i -= 1;
  }
  assert( pos <= ns );

  delete [] rawscr;
  delete [] score;
  
  return pos;
}

bool
readProfile(PyObject*   pProfile,
	    uint const  nProfile,
	    int*        storage,
	    int**       profile,
	    uint        pad = 0)
{
  for(uint n = 0; n < nProfile; ++n) {
    PyObject* pn = PySequence_Fast_GET_ITEM(pProfile, n);
    if( ! (PySequence_Check(pn) && PySequence_Size(pn) == 6) ) {
      // Leak
      PyErr_SetString(PyExc_ValueError, "wrong args: bad profile") ;
      return false;
    }
    profile[n+pad] = storage + 6*(n+pad);
    for(uint i = 0; i < 6; ++i) {
      long const v = PyInt_AsLong(PySequence_Fast_GET_ITEM(pn, i));
      if( v == -1 ) {
	// Leak
	PyErr_SetString(PyExc_ValueError, "wrong args: bad profile");
	return false;
      }
      profile[n+pad][i] = v;
    }
  }
  return true;
}

// Copies profile cell p[i-1] to al, adds n gaps, decrementing i and
// incrementing al.
inline void
profileFillMatchedGap(const int* const* p, int& i, int*& al, int const n)
{
  i -= 1;
  const int* pr = p[i];
  for(uint k = 0; k < 6; ++k) {
    *al++ = *pr++;
  }
  al[-1] += n;  // assumes gap is last
}

template<typename T>
uint
alignProfToProf(const int* const*          p1,
		uint const                 lp1,
		const int* const*          p2,
		uint const                 lp2,
		MatchScoreValues<T> const  matchScores,
		int* const                 alignment)
{
  int* al0 = alignment;
  
  uint const sz = (lp1+1)*(lp2+1);
  T* score = new T [sz];
  std::fill(score, score+sz, 0);

  int n1 = 0, n2 = 0;
  for(int i = 0; i < 6; ++i) {
    n1 += p1[0][i];
    n2 += p2[0][i];
  }

  ProfileCounts colProfCount[lp2];
  for(uint j = 0; j < lp2; j += 1) {
    colProfCount[j] = ProfileCounts(p2[j]);
  }

  ProfileCounts countPerType[lp1][6];
  ProfileCounts rowProfCount[lp1];
  
  for(uint i = 0; i < lp1; i += 1) {
    const int* p1i = p1[i];
    rowProfCount[i] = ProfileCounts(p1i);

    {
      int const an = p1i[anynuc];
      int const gp  = p1i[gap];
      int const nonGap = n1 - gp;
      int const nc = nonGap - an;    assert (  nc == p1i[0]+p1i[1]+p1i[2]+p1i[3] );
      
      for(int k = 0; k < 4; ++k) {
	auto& ck = countPerType[i][k];
	ck.matches = p1i[k] + an;
	ck.mis = nc - p1i[k];
	ck.gaps = gp;
      }
      
      auto& anc = countPerType[i][anynuc];
      anc.matches = nonGap;
      anc.mis = 0;
      anc.gaps = gp;

      auto& cgp = countPerType[i][gap];
      cgp.matches = cgp.mis = 0;
      cgp.gaps = nonGap;
    }
  }

  T* const colGapScore = new T [lp2+lp1];    std::unique_ptr<T> rel0(colGapScore);
  T* const rowGapScore = colGapScore + lp2;
  
  {
    T penalty = 0;
    T* s = score + 1;
    for(uint jCol = 1; jCol <= lp2; ++jCol, ++s) {
      //colGapScore[jCol-1] = matchScores.scoreGap(p2[jCol-1], n1, n2);
      colGapScore[jCol-1] = matchScores.scoreGap(colProfCount[jCol-1], n1, n2);      
      penalty += colGapScore[jCol-1];
      *s = penalty;
    }
  }
  
  {
    T penalty =  0;
    T* s = score + (lp2+1);
    for(uint iRow = 1; iRow <= lp1; ++iRow) {
      //rowGapScore[iRow-1] = matchScores.scoreGap(p1[iRow-1], n2, n1);
      rowGapScore[iRow-1] = matchScores.scoreGap(rowProfCount[iRow-1], n2, n1);
      penalty += rowGapScore[iRow-1];
      *s = penalty;
      s += lp2+1;
    }
  }

  
  T* m1m1 = score;
  for(uint i = 1; i <= lp1; ++i) {
    
    T const gs = rowGapScore[i-1];
#if defined(ALLASSERTS)
    //                                          assert(rowGapScore[i-1] == matchScores.scoreGap(p1i, n2, n1));
#endif
    for(uint j = 0; j < lp2; j += 1) {

      T const sm = matchScores.scoreMatching(p2[j], rowProfCount[i-1], colProfCount[j], countPerType[i-1]);
#if defined(ALLASSERTS)
      T const sm1 = matchScores.scoreMatching(p1[i-1], n1, p2[j], n2); assert( sm == sm1 );
#endif
      
      T const match = *m1m1 + sm;
      T const del = *(m1m1+1) + gs;
      T const ins = *(m1m1+1+lp2) + colGapScore[j];
#if defined(ALLASSERTS)
      //                                       assert(colGapScore[j] == matchScores.scoreGap(p2[j], n1, n2));
#endif

      *(m1m1+lp2+2) = std::max(std::max(match, del), ins);
      m1m1 += 1;
    }
    m1m1 += 1;
  }
  
  int iRow = lp1;
  int jCol = lp2;

  if( matchScores.freeEndGaps ) {
    uint iMaxRow = sz-lp2-1;
    T mxLastRow = score[iMaxRow];
    for(uint l = iMaxRow+1; l < sz; ++l) {
      if( score[l] >= mxLastRow ) {
	mxLastRow = score[l];
	iMaxRow = l;
      }
    }

    uint iMaxCol = lp2;
    T mxLastCol = score[iMaxCol];
    for(uint l = iMaxCol+lp2+1; l < sz; l += lp2+1) {
      if( score[l] >= mxLastCol ) {
	mxLastCol = score[l];
	iMaxCol = l;
      }
    }

    if( mxLastCol > mxLastRow ) {
      int const im = (iMaxCol - lp2)/(lp2+1);
      while( iRow > im ) {
	profileFillMatchedGap(p1, iRow, al0, n2);
      }
    } else {
      int const jm = iMaxRow - (sz-lp2-1);
      while( jCol > jm ) {
	profileFillMatchedGap(p2, jCol, al0, n1);
      }
    }
  }
  
  while( iRow > 0 and jCol > 0 ) {
    uint const cur = iRow*(lp2+1) + jCol;
    T const score_current = score[cur];
    T const score_diagonal = score[cur - lp2 - 2];

    T const sm = matchScores.scoreMatching(p2[jCol-1], rowProfCount[iRow-1],
					   colProfCount[jCol-1], countPerType[iRow-1]);
#if defined(ALLASSERTS)
    T const sm1 = matchScores.scoreMatching(p1[iRow-1], n1, p2[jCol-1], n2);
    assert ( sm == sm1 );
#endif
    
    if( score_current == score_diagonal + sm ) {
      const int* p1r = p1[iRow-1];
      const int* p2c = p2[jCol-1];
      for(int i = 0; i < 6; ++i) {
	al0[i] = p1r[i] + p2c[i];
      }
      iRow -= 1;
      jCol -= 1;
      al0 += 6;
    } else {
      T const score_up = score[cur - lp2 - 1];

      if( score_current == score_up + rowGapScore[iRow-1] ) {
	profileFillMatchedGap(p1, iRow, al0, n2);
      } else {
#if !defined(NDEBUG)
	T const score_left = score[cur - 1];
#endif
	assert ( score_current == score_left + colGapScore[jCol-1] ) ;

	profileFillMatchedGap(p2, jCol, al0, n1);
      }
    }
  }
  
  while( iRow > 0 ) {
    profileFillMatchedGap(p1, iRow, al0, n2);
  }
  
  while( jCol > 0 ) {
    profileFillMatchedGap(p2, jCol, al0, n1);
  }

  delete [] score;
  
  return (al0 - alignment)/6;
}

PyObject*
alignToProfile(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"seq", "profile", "pad", "chop",
				 "matchScore", "misMatchScore", "gapPenalty",
				 static_cast<const char*>(0)};
  PyObject* pSeq = 0;
  PyObject* pPfofile = 0;
  float matchScore = 10;
  float misMatchScore = -5;
  float gapPenalty = -6;
  int	pad = 0;
  PyObject* pChop = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "OO|iOfff", const_cast<char**>(kwlist),
				    &pSeq, &pPfofile,&pad,&pChop,
				    &matchScore, &misMatchScore, &gapPenalty)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  uint seqLen;
  byte* seq = readSequence(pSeq, seqLen, false);
  if( ! seq ) {
    return 0;
  }

  if( ! PySequence_Check(pPfofile) || pad < 0 ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: bad profile/pad") ;
    return 0;
  }

  int const nProfile = PySequence_Size(pPfofile);
  uint const nSites = nProfile + 2 * pad;
  
  int*  const allp(new int [6*nSites]);      std::unique_ptr<int> rel(allp);
  int** const profile(new int* [nSites]);    std::unique_ptr<int*> rel2(profile);

  if( ! readProfile(pPfofile, nProfile, allp, profile, pad) ) {
    return 0;
  }

  int nSeqs = 0;
  for(uint i = 0; i < 6; ++i) {
    nSeqs += profile[pad][i];
  }
  
  for(int n = 0; n < pad; ++n) {
    profile[n] = allp + 6*n;
    int const n1 = pad + nProfile + n;
    profile[n1] = allp + 6*n1;
    for(uint i = 0; i < gap; ++i) {
      profile[n][i] = profile[n1][i] = 0;
    }
    profile[n][gap] = profile[n1][gap] = nSeqs;
  }
    
  if( seqLen < nSites ) {
    byte* s = new byte[nSites];
    for(uint k = 0; k < seqLen; ++k) {
      s[k] = seq[k];
    }
    delete seq;
    seq = s;
  }
  
  int const newLen = alignToProf<float>(seq, seqLen, profile, nSites,
					matchScore, misMatchScore, gapPenalty);
  
  PyObject* retSeq = 0;

  bool const chop = (pad > 0 && pChop != 0 && PyObject_IsTrue(pChop));
  
  if( newLen < 0 ) {
    PyErr_SetString(PyExc_ValueError, "short profile.");
  } else {
    int chopFront = 0;
    int chopBack = 0;
    if( pad > 0 ) {
      if( chop ) {
	while( seq[newLen - 1 - chopFront] == gap && chopFront < pad ) {
	  chopFront += 1;
	}
      
	while( seq[chopBack] == gap && chopBack < pad ) {
	  chopBack += 1;
	}
      }
    }

    int const totChop = chopBack+chopFront;
    retSeq = PyTuple_New(newLen - totChop);
    for(int i = chopFront; i < newLen-chopBack; ++i) {
      PyTuple_SET_ITEM(retSeq, i-chopFront, PyInt_FromLong(seq[newLen - 1 - i]));
    }
    
    if( chop ) {
      PyObject* t = PyTuple_New(3);
      PyTuple_SET_ITEM(t, 0, retSeq);
      PyTuple_SET_ITEM(t, 1, PyInt_FromLong(pad - chopFront));
      PyTuple_SET_ITEM(t, 2, PyInt_FromLong(pad - chopBack));
      retSeq = t;
    }
  }
  
  delete [] seq;

  return retSeq;
}


PyObject*
alignProfileToProfile(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"p0", "p1", "scores",
				 static_cast<const char*>(0)};
  PyObject* pPfofile0 = 0;
  PyObject* pPfofile1 = 0;
  PyObject* mScores = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "OO|O", const_cast<char**>(kwlist),
				    &pPfofile0, &pPfofile1, &mScores)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PySequence_Check(pPfofile0) && PySequence_Check(pPfofile1)) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: bad profiles") ;
    return 0;
  }

  int const nProfile0 = PySequence_Size(pPfofile0);
  int const nProfile1 = PySequence_Size(pPfofile1);
  
  int*  const allp(new int [12*(nProfile0+nProfile1)]);   std::unique_ptr<int> rel0(allp);
  int** const profile0(new int* [nProfile0]);             std::unique_ptr<int*> rel1(profile0);
  int** const profile1(new int* [nProfile1]);             std::unique_ptr<int*> rel2(profile1);

  if( ! readProfile(pPfofile0, nProfile0, allp, profile0) ) {
    return 0;
  }
  
  if( ! readProfile(pPfofile1, nProfile1, allp + 6*nProfile0, profile1) ) {
    return 0;
  }
  
  MatchScoreValues<double> const matchScores(mScores);

  int* const al = allp + 6*(nProfile0+nProfile1); //  std::unique_ptr<int> rel3(al);

  // Scores in the profile grow large and run down the float accuracy.
  int const alLen = alignProfToProf<double>(profile0, nProfile0, profile1, nProfile1, matchScores, al);
  
  PyObject* retSeq = PyTuple_New(alLen);
  for(int i = 0; i < alLen; ++i) {
    PyObject* t = PyTuple_New(6);
    const int* const a = al + 6*(alLen - 1 - i);
    for(int k = 0; k < 6; ++k) {
      PyTuple_SET_ITEM(t, k, PyInt_FromLong(a[k]));
    }
    PyTuple_SET_ITEM(retSeq, i, t);
  }

  return retSeq;
}

PyObject*
createProfile(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"seqs",
				 static_cast<const char*>(0)};
  PyObject* pSeqs = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "O", const_cast<char**>(kwlist),
				    &pSeqs) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PyTuple_Check(pSeqs) || PyList_Check(pSeqs)) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not a list of sequences") ;
    return 0;
  }

  int const nSeqs = PySequence_Size(pSeqs);

  if( nSeqs == 0 ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: no sequences") ;
    return 0;
  }
  
  uint seqLen = 0;
  long* counts = 0;
  
  for(int j = 0; j < nSeqs; ++j) {
    PyObject* seq = PySequence_Fast_GET_ITEM(pSeqs, j);

    uint lseq;
    const byte* const s = readSequence(seq, lseq, false);
    if( ! s ) {
      return 0;
    }
    
    if( j == 0 ) {
      seqLen = lseq;
      counts = new long [seqLen*6];
      for(long* c = counts; c < counts+seqLen*6; ++c) {
	*c = 0;
      }
    } else {
      if( unsigned(lseq) != seqLen ) {
	delete [] counts;
	delete [] s;
	PyErr_SetString(PyExc_ValueError, "wrong args: not aligned");
	return 0;
      }
    }
    
    long* c = counts;
    for(const byte* seq = s; lseq > 0; --lseq, ++seq) {
      c[*seq] += 1;
      c += 6;
    }
    delete [] s;
  }
  
  PyObject* retSeq = PyTuple_New(seqLen);
  const long* c = counts;
  for(uint i = 0; i < seqLen; ++i, c += 6) {
    PyObject* t = PyList_New(6);
    for(uint j = 0; j < 6; ++j) {
      //PyList_SET_ITEM(t, j, PyLong_FromLong(c[j]));
      PyList_SET_ITEM(t, j, PyInt_FromLong(c[j]));      
    } 
    PyTuple_SET_ITEM(retSeq, i, t);
  }

  delete [] counts;
  return retSeq;
}

template<typename T>
struct ColumnMin {
  ColumnMin(void) :
    value(std::numeric_limits<T>::infinity()),
    index(-1)
    {}
  
  T		value;
  int		index;
};

template<typename T>
PyObject*
upgma(T* const ds, uint const n, const int* const weights)
{
  // keep a pointer for each row. Only O(n)
  
  T** const di = new T* [n];

  if( ds == 0 || di == 0 ) {
    return PyErr_NoMemory();
  }

  {
    int o = 0;
    for(int k = n-1; k > 0; k -= 1) {
      di[n - 1 - k]= ds + o;
      o += k;
    }
  }

  // w[i] the weight (number of elements) in live group i
  vector<int> w(n, 1);
  if( weights ) {
    std::copy(weights, weights+n, w.begin());
  }
  
  // Index/name of live columns (scipy style).
  vector<int> cx(n);
  // Remaining columns mapped to initial upper half matrix
  vector<int> col(n-1);

  for(uint i = 0; i < n; ++i) {
    cx[i] = i;
    if( i+1 < n ) {
      col[i] = i;
    }
  }

  // minimum of column. Global minimum is the minimum of those.
  // We try to keep the column minimum current as much as possible
  // during the update.
  vector< ColumnMin<T> > colMins(n-1);
  for(uint mj = 0; mj < n-1; ++mj) {
    ColumnMin<T>& c = colMins[mj];
    for(uint mi = 0; mi < mj+1; ++mi) {
      if( di[mi][mj-mi] < c.value ) {
	c.value = di[mi][mj-mi];
	c.index = mi;
      }
    }
  }
  
  PyObject* ret = PyTuple_New(n-1);
  for(uint ic = 0; ic < n-1; ++ic) {
    uint lcol = n-1-ic;

    // Find global minimum. Update stale column minimum if needed.
    T mn = std::numeric_limits<T>::infinity();
    int mi = -1,mj = -1;
    
    for(uint xmj = 1; xmj < lcol+1; ++xmj) {
      ColumnMin<T>& c = colMins[xmj-1];
      if( c.index == -1 ) {
	for(uint xmi = 0; xmi < xmj; ++xmi) {
	  int const r = xmi > 0 ? col[xmi-1]+1 : 0;
	  T const v = di[r][col[xmj-1] - r];
	  if( v < c.value ) {
	    c.value = v;
	    c.index = xmi;
	  }
	}
      }
      if( c.value < mn ) {
	mn = c.value;
	mj = xmj;
	mi = c.index;
      }
    }

    // global minimum info
    int const wi = w[mi];
    int const wj = w[mj];
    T const tw = wi+wj;
    T const rwi = w[mi]/tw;
    T const rwj = w[mj]/tw;

    int const cmj = col[mj-1];
    int const cmi = col[mi-1];

    // merge i and j to i. Start with updating cells in i's column,
    // that is i,k with k < i. Keep minimum of updates as we go.
    ColumnMin<T> mm;
    
    for(int ki = 0; ki < mi; ++ki) {
      int const ir = ki > 0 ? col[ki-1]+1 : 0;
      T* const d = di[ir];
      int const i = cmi - ir;
      T const nv = (rwi * d[i] + rwj * d[cmj - ir]);
      if( nv < mm.value ) {
	mm.value = nv;
	mm.index = ki;
      }
      d[i] = nv;
    }

    // 0 never goes. If column minimum was among the updates, and no value is
    // smaller, invalidate the column minimum since it might be in the rest of
    // the distances.
    if( mi > 0 ) {
      if( mm.value < colMins[mi-1].value ) {
	colMins[mi-1] = mm;
      } else if( colMins[mi-1].index < mi ) {
	colMins[mi-1] = ColumnMin<T>();
      }
    }

    // Now do the row. Invalidate column minimum if we update the cell.
    // Duplicate code for speed, part one for cells smaller than j, the second
    // for larger cells.

    int const ir = mi > 0 ? cmi+1 : 0;
    T* const d = di[ir];

    for(int ki = mi; ki < mj-1; ++ki) {
      //assert 0 < ki < mj

      int const i = col[ki] - ir;
      int const ir1 = col[ki]+1;
      T const nv = (rwi * d[i] + rwj * di[ir1][cmj - ir1] );
	
      if( nv <= colMins[ki].value ) {
	colMins[ki].index = mi;
	colMins[ki].value = nv;
      } else if( d[i] == colMins[ki].value ) {
	colMins[ki] = ColumnMin<T>();
      }
      d[i] = nv;
    }
    
    int const ir1 =  mj > 0 ? cmj+1 : 0;
    const T* d2 = uint(mj) < lcol ? di[ir1] : 0;
    for(uint ki = mj; ki < lcol; ++ ki) {
      //assert 0 <= mj < ki

      int const c = col[ki];
      int const i = c - ir;

      T const nv = (rwi * d[i] + rwj * d2[c - ir1]);
      
      if( nv <= colMins[ki].value ) {
	colMins[ki].index = mi;
	colMins[ki].value = nv;
      } else if( d[i] == colMins[ki].value ) {
	colMins[ki] = ColumnMin<T>();
      }

      d[i] = nv;
    }

    // Merge now
    
    w[mi] += w[mj];
    w.erase(w.begin() + mj);

    PyObject* r = PyTuple_New(4);
    PyTuple_SET_ITEM(r, 0, PyInt_FromLong(std::min(cx[mi],cx[mj])) );
    PyTuple_SET_ITEM(r, 1, PyInt_FromLong(std::max(cx[mi],cx[mj])) );
    PyTuple_SET_ITEM(r, 2, PyFloat_FromDouble(mn) );
    PyTuple_SET_ITEM(r, 3, PyInt_FromLong(wi+wj));

    PyTuple_SET_ITEM(ret, ic, r);
    
    cx[mi] = ic + n;

    cx.erase(cx.begin() + mj);
    col.erase(col.begin() + mj-1);

    // Minimums in the deleted row (of j) invalidate the row minimum
    for(uint xmj = 1; xmj < lcol+1; ++xmj) {
      ColumnMin<T>& c = colMins[xmj-1];
      if( c.index == mj ) {
        c = ColumnMin<T>();
      }
      if( c.index > mj ) {
        c.index -= 1;
      }
    }
    colMins.erase(colMins.begin() + mj-1);
  }
  delete [] di;
  
  return ret;
}

PyObject*
UPGMA(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"distances", "seqs", "align", "report",
				 "saveto", "reorder", "scores", "weights",
				 static_cast<const char*>(0)};
  PyObject* dists = 0;
  PyObject* pseqs = 0;
  PyObject* palign = 0;
  ComparisonResult resultType = DIVERGENCE;
  PyObject* saveDistancesTo = 0;
  PyObject* pReorder = 0;
  PyObject* mScores = 0;
  PyObject* pWeights = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "|OOOiOOOO", const_cast<char**>(kwlist),
				    &dists, &pseqs, &palign,
				    &resultType, &saveDistancesTo,&pReorder,mScores,
				    &pWeights)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( (!dists) == (!pseqs) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: expecting either distances or sequences");
    return 0;
  }

  float* ds = 0;
  uint n;

  std::unique_ptr<float> dsReleaser;
  
  if( pseqs ) {
    if( ! (PyTuple_Check(pseqs) || PyList_Check(pseqs)) ) {
      PyErr_SetString(PyExc_ValueError, "wrong args: not a list of sequences") ;
      return 0;
    }

    if( ! (STATS < resultType && resultType <= JCcorrection) ) {
      PyErr_SetString(PyExc_ValueError, "wrong args (return type)");
      return 0;
    }

    if( saveDistancesTo ) {
      if( ! (PyFile_Check(saveDistancesTo) || 
	     PyObject_HasAttrString(saveDistancesTo, "append")) ) {
	PyErr_SetString(PyExc_ValueError, "can't save");
	return 0;
      }
    }
    
    n = PySequence_Size(pseqs);
    uint const nDists = (n*(n-1))/2;
    ds = new float [nDists];             dsReleaser.reset(ds);
	
    if( ! distmat(pseqs, palign, resultType, pReorder, mScores, ds) ) {
      return 0;
    }

    if( saveDistancesTo ) {
      if( PyFile_Check(saveDistancesTo) ) {
	for(uint i = 0; i < nDists; ++i) {
	  char* const b = PyOS_double_to_string(ds[i], 'r', 0, Py_DTSF_ADD_DOT_0, 0);
	  if( !b || PyFile_WriteString(b, saveDistancesTo) != 0 ) {
	    return 0;
	  }
	  PyFile_WriteString("\n", saveDistancesTo);
	  PyMem_Free(b);
	}
      } else {
	PyObject* const ap = PyObject_GetAttrString(saveDistancesTo, "append");
	assert( ap );
	  
	for(uint i = 0; i < nDists; ++i) {
	  PyObject* pt = PyTuple_New(1);
	  //PyObject* f = PyFloat_FromDouble(ds[i]);
	  // value is used by the tuple, so we lose the reference
	  PyTuple_SET_ITEM(pt, 0, PyFloat_FromDouble(ds[i]));
	  PyObject* r = PyObject_CallObject(ap, pt);
	  Py_DECREF(r);
	  Py_DECREF(pt);
	}
	Py_DECREF(ap);
      }
    }
  } else {
    if( ! PySequence_Check(dists) ) {
      PyErr_SetString(PyExc_ValueError, "wrong args: not a list of distances ") ;
      return 0;
    }

    PyObject* rel = 0;
    if( ! (PyTuple_Check(dists) || PyList_Check(dists)) ) {
      dists = PySequence_Fast(dists, "strange sequence");
      rel = dists;
    }
    
    unsigned long const nDists = PySequence_Size(dists);
    n = static_cast<uint>(.5 + (sqrt(1+8*nDists)+1)/2);

    if( nDists == 0 || ((long(n)*(n-1))/2 != nDists) ) {
      PyErr_SetString(PyExc_ValueError, "wrong args: incorrect dists") ;
      return 0;
    }
  
    ds = new float [nDists];             dsReleaser.reset(ds);
    
    for(uint i = 0; i < nDists; ++i) {
      PyObject* d = PySequence_Fast_GET_ITEM(dists, i);
      ds[i] = PyFloat_AsDouble(d);
      if( ds[i] < 0 ) {
	PyErr_SetString(PyExc_ValueError, "wrong args: incorrect diststances") ;
	return 0;
      }
    }
    
    Py_XDECREF(rel);
  }
 
  int* weights = 0;
  if( pWeights and pWeights != Py_None ) {
    if( ! (PySequence_Check(pWeights) && uint(PySequence_Size(pWeights)) == n) ) {
      PyErr_SetString(PyExc_ValueError, "Error with weights");
      return 0;
    }
    weights = new int[n];
      
    for(uint i = 0; i < n; ++i) {
      PyObject* o = PySequence_Fast_GET_ITEM(pWeights, i);
      weights[i] = PyInt_AsLong(o);
	
      if( weights[i] <= 0 ) {
	PyErr_SetString(PyExc_ValueError, "Error with weights (should be positive)");
	delete [] weights;
	return 0;
      }
    }
  }

  PyObject* ret = upgma<float>(ds, n, weights);

  delete [] weights;
  
  return ret;
}

PyDoc_STRVAR(calign__doc__,
"Sequence alignment and UPGMA");

PyDoc_STRVAR(upgma__doc__,
	     "UPGMA tree from distances or sequences. Return a scipy compatible list." 
	     " 'saveto' can be either an open file or an object supporting an append.");

static PyMethodDef calignMethods[] = {
  {"globalAlign",	(PyCFunction)globAlign, METH_VARARGS|METH_KEYWORDS,
   "Global alignment of two DNA sequences. Full Needleman-Wunch with a free flanking gaps option."},
  
  {"profileAlign",	(PyCFunction)alignToProfile, METH_VARARGS|METH_KEYWORDS,
   ""},
  {"prof2profAlign",	(PyCFunction)alignProfileToProfile, METH_VARARGS|METH_KEYWORDS,
   ""},
  {"createProfile",	(PyCFunction)createProfile, METH_VARARGS|METH_KEYWORDS,
   "Profile from alignment."},
  
  {"distances",		(PyCFunction)distMat, METH_VARARGS|METH_KEYWORDS,
   "Distances for all 'n choose 2' pairs (via alignment)."},
  {"allpairs",		(PyCFunction)distPairs, METH_VARARGS|METH_KEYWORDS,
   "Distances for all NxM pairs (via alignment)."},
  
  {"upgma",		(PyCFunction)UPGMA, METH_VARARGS|METH_KEYWORDS,upgma__doc__},  
  {NULL, NULL, 0, NULL}        /* Sentinel */
};


static inline void
setInt(PyObject* m, const char* const attr, int v)
{
  PyObject* const o = PyInt_FromLong(v);
  PyObject_SetAttrString(m, const_cast<char*>(attr), o);
  Py_DECREF(o);
}

PyMODINIT_FUNC
initcalign(void)
{
  PyObject* const m = Py_InitModule3("calign", calignMethods, calign__doc__);

  setInt(m, "GAP", gap);
  setInt(m, "N", anynuc);
  setInt(m, "A", 0);
  setInt(m, "G", 1);
  setInt(m, "C", 2);
  setInt(m, "T", 3);

  setInt(m, "JCcorrection", JCcorrection);
  setInt(m, "IDENTITY", IDENTITY);
  setInt(m, "DIVERGENCE", DIVERGENCE);
  setInt(m, "STATS", STATS);

  //Choose2::init(50000);
}


// ('ACATGTCTTCGGTGGACTTCAGTTCCGCCCGAGAAACTTAGACAA', 'AATGTAATAAATGG')
// order with free gaps can matter

