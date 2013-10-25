// This file is part of biopy.
// Copyright (C) 2013 Joseph Heled
// Author: Joseph Heled <jheled@gmail.com>
// Author: Alexei Drummond <alexei@cs.auckland.ac.nz>
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

#include "readseq.h"

// last 3 are ambiguous codes
static const char aminoAcidsOrder[24] = "ACDEFGHIKLMNPQRSTVWYBXZ";
static uint const nAA = sizeof(aminoAcidsOrder)-1;

static byte*
readAASequence(PyObject* seq, uint& nseq)
{
  if( ! PySequence_Check(seq) ) {
    PyErr_SetString(PyExc_ValueError, "not a sequence") ;
    return 0;
  }

  nseq = PySequence_Size(seq);
  byte* s = new byte[nseq];

  if( PyString_Check(seq) ) {
    const char* const aa = PyString_AsString(seq);
    for(uint i = 0; i < nseq; ++i) {
      char const c = aa[i];
      uint k = 0;
      for(/**/; k < nAA; ++k) {
	if( c == aminoAcidsOrder[k] ) {
	  break;
	}
      }
      if( k < nAA ) {
	s[i] = k;
      } else {
	PyErr_SetString(PyExc_ValueError, "Invalid AA sequence") ;
	delete [] s;
	return 0;
      }
    }
  } else {
    for(uint i = 0; i < nseq; ++i) {
      PyObject* const o = PySequence_Fast_GET_ITEM(seq, i);
      
      s[i] = PyInt_AsLong(o);
      
      if( !( 0 <= s[i] && s[i] < nAA ) ) {
	delete [] s;
	PyErr_SetString(PyExc_ValueError, "Invalid AA sequence") ;
	return 0;
      }
    }
  }
  return s;
};

template<typename T>
class CorrectionScores {
public:
  T const indelPenalty;
  T const correctionPenalty;
  T const stopCodonPenalty;

  CorrectionScores(T d, T c, T s) :
    indelPenalty(d),
    correctionPenalty(c),
    stopCodonPenalty(s)
    {}
};

enum FType {
  match = 0,
  ins_read,
  match_duplicate,
  match_delete,
  ins_read_duplicate,
  ins_read_delete,
  ins_ref
};

enum CodonLocations {
  match_ins = 0,
  dup_0,
  dup_1,
  dup_2,
  del_0,
  del_1,
  del_2,
  del_3,
};

static int const
cd[del_3+1][3] = {
  {-2, -1, 0},
  
  {-1, 0, 0},
  {-1, -1, 0},
  {-2, -1, 0},

  {-3, -2, -1},
  {-3, -2, 0},
  {-3, -1, 0},
  {-2, -1, 0}
};

struct Cell {
  FType          what : 16;
  CodonLocations loc : 16;
  float          score;
};

typedef int Codon[3];

/**
 * Provides scoring for alignment + correction
 */
class FScore {
  const int* const geneticCode;
  const int* const scores;

public:
  CorrectionScores<float> const& correctionScores;

  FScore(const int*                      _geneticCode,
	 CorrectionScores<float> const&  _correctionScores,
	 const int*                      _scores);

  int translate(uint const v) const;
  
  uint getAAIndicesOfCodon(Codon const codon, int (&aas)[16]) const;
  uint getAAofFuzzyCodon(Codon const codon, int (&aas)[16]) const;
    
  double matchScore(Codon const codon, int j) const;
  double insertCodon(Codon const codon) const;
  double matchDelete(Codon const codon, int j) const;
  double insertCodonDelete(Codon const codon) const;
  double matchDuplicate(Codon const codon, int j) const;
};

FScore::FScore(const int*                           _geneticCode,
	       CorrectionScores<float> const&  _correctionScores,
	       const int*                           _scores) :
  geneticCode(_geneticCode),
  scores(_scores),
  correctionScores(_correctionScores)
{}

inline int
FScore::translate(uint const v) const
{
  return geneticCode[v];
}

inline uint
FScore::getAAofFuzzyCodon(Codon const codon, int (&aas)[16]) const
{
  uint n = 0;
  int const c0 = codon[0];
  int const c1 = codon[1];
  int const c2 = codon[2];
  bool const fz0 = ! ( 0 <= c0 && c0 < 4 );
  bool const fz1 = ! ( 0 <= c1 && c1 < 4 );
  bool const fz2 = ! ( 0 <= c2 && c2 < 4 );

  if( fz0 && fz1 && fz2 ) {
    return 0;
  }
  
  for(int k0 = 0; k0 < 4; ++k0) {
    for(int k1 = 0; k1 < 4; ++k1) {
      for(int k2 = 0; k2 < 4; ++k2) {
	if( (fz0 || k0 == c0) && (fz1 || k1 == c1)  && (fz2 || k2 == c2) ) {
	  aas[n] = 16*k0+4*k1+k2;
	  n += 1;
	}
      }
    }
  }
  return n;
}

inline uint
FScore::getAAIndicesOfCodon(Codon const codon, int (&aas)[16]) const
{
  uint n;
  int const c0 = codon[0];
  if( c0 > 3 || codon[1] > 3 || codon[2] > 3 ) {
     n = getAAofFuzzyCodon(codon, aas);
  } else {
    // no N's
    if( c0 >= 0 ) {
      int const c1 = codon[1];
      int const c2 = codon[2];           assert( c1 >= 0 && c2 >= 0 );
      aas[0] = 16*c0 + 4*c1 + c2;
      n = 1;
    } else {
      assert(codon[2] >= 0);
      if( codon[1] >= 0 ) {
	int const b = 4*codon[1] + codon[2];
	for(int i = 0; i < 4; ++i) {
	  aas[i] = 16*i + b;
	}
	n = 4;
      } else {
	for(int i = 0; i < 4; ++i) {
	  for(int j = 0; j < 4; ++j) {
	    aas[4*i+j] = 16*i +4*j + codon[2];
	  }
	}
	n = 16;
      }
    }
  }
  return n;
}
    
inline double
FScore::matchScore(Codon const codon, int j) const
{
  double score = 0;
  int aas[16];
  uint const n = getAAIndicesOfCodon(codon, aas);
  if( n == 1 ) {
    int const a = geneticCode[aas[0]];
    score = (a >= 0) ? scores[nAA * a + j] : correctionScores.stopCodonPenalty;
  } else {
    uint nonStop = 0;
    for(uint i = 0; i < n; ++i) {
      int const a = geneticCode[aas[i]];
      if( a >= 0 ) {
	score += scores[nAA * a + j];
	nonStop += 1;
      }
    }
    if( nonStop > 0 ) {
      score /= nonStop;
    } else {
      if( n > 0 ) {
	score = correctionScores.stopCodonPenalty;
      } else {
	// NNN
	for(uint i = 0; i < 64; ++i) {
	  int const a = geneticCode[i];
	  if( a >= 0 ) {
	    score += scores[nAA * a + j];
	    nonStop += 1;
	  }
	}
	score /= nonStop;
      }
    }
  }
  return score;
}

inline double
FScore::insertCodon(Codon const codon) const
{
  int aas[16];
  uint const n = getAAIndicesOfCodon(codon, aas);
  if( n == 0 ) {
    // all N
    return correctionScores.indelPenalty;
  }
  
  for(uint i = 0; i < n; ++i) {
    if( geneticCode[aas[i]] > 0 ) {
      return correctionScores.indelPenalty;
    }
  }
  // stop Codon
  return correctionScores.stopCodonPenalty + correctionScores.indelPenalty;
}

inline double
FScore::matchDelete(Codon const codon, int j) const
{
  return matchScore(codon, j) + correctionScores.correctionPenalty;
}

inline double
FScore::insertCodonDelete(Codon const codon) const
{
  return insertCodon(codon) + correctionScores.correctionPenalty;
}

inline double
FScore::matchDuplicate(Codon const codon, int j) const
{
  return matchScore(codon, j) + correctionScores.correctionPenalty;
}


class AlignAndCorrect {
public:
  FScore const scores;

  uint nRead;
  const byte* read;
  uint naa;
  const byte* aa;

  // [i,j] is the score of the best alignment of i nucs [0.. i-1] and j AA [0 .. j-1]
  Cell** nodes;

  struct Result {
    uint matches;
    uint mismatches;
    uint gaps;
    uint correctionDeletions;
    uint correctionInsertions;

    uint dnaFreeEnd;
    uint aaFreeEnd;
    uint dnaFreeStart;
    uint aaFreeStart;
    
    vector<int>		alignedFramedRead;
    vector<int>		alignedAA;
    vector<int>		dnaBoundries;
    
    Result(void) :
      matches(0),
      mismatches(0),
      gaps(0),
      correctionDeletions(0),
      correctionInsertions(0)
      {}
  };

  AlignAndCorrect(const int* sub, CorrectionScores<float> const& cscores, int const geneticCode[]);

  inline void
  getCodon(CodonLocations const cloc, int const idna, int (&codon)[3]) const
  {
    auto const cdl = cd[cloc];
    uint const ii = idna - 1;
    for(uint k = 0; k < 3 ; ++k) {
      int const i = ii + cdl[k];
      codon[k] = i >= 0 ? read[i] : -1;
    }
  }

  Result doAlignment(const byte* read, uint nRead, const byte* aa, uint naa);

  void scoreCell(uint idna, uint jaa);
};

AlignAndCorrect::AlignAndCorrect(const int* sub,
				 CorrectionScores<float> const& cscores,
				 int const geneticCode[]) :
  scores(geneticCode, cscores, sub),
  nRead(0),
  naa(0),
  nodes(0)
{}


void
AlignAndCorrect::scoreCell(uint const idna, uint const jaa)
{
  Cell& c = nodes[idna][jaa];
    
  if (idna < 2 || jaa == 0) {
    c.score = 0;
    c.what = /*some kind of none;*/ (FType)200;
  } else {
    // gap in AA
    double const insAA = nodes[idna][jaa - 1].score + scores.correctionScores.indelPenalty;
    c.score = insAA;
    c.what = ins_ref;

    int curaa = aa[jaa-1];
    
    if( idna >= 3 ) {
      // no insertions or deleteions to codon
      int codon[3];
      getCodon(match_ins, idna, codon);    
      double const sc = scores.matchScore(codon, curaa);
      // match to AA
      double const matchNuc = nodes[idna-3][jaa - 1].score + sc;
      if( matchNuc > c.score ) {
	c.score = matchNuc;
	c.what = match;
	c.loc = match_ins;
      }
	
      // codon not matching - insert gap in AA
      double const insNuc = nodes[idna-3][jaa].score + scores.insertCodon(codon);
      if( insNuc > c.score ) {
	c.score = insNuc;
	c.what = ins_read;
	c.loc = match_ins;
      }	  
    }

    if( idna >= 4 ) {
      // nuc was inserted to frame. 4 possible ways	
      int codon[3];
      for(uint m = del_0; m <= del_3 ; ++m) {
	getCodon(static_cast<CodonLocations>(m), idna, codon);    
	
	// nuc was inserted to frame and codon matched to AA
	double const delAndMatch = nodes[idna-4][jaa - 1].score + scores.matchDelete(codon, curaa);
	if( delAndMatch > c.score ) {
	  c.score = delAndMatch;
	  c.what = match_delete;
	  c.loc = static_cast<CodonLocations>(m);
	}
	  
	// nuc was inserted to frame and codon deleted
	double const insAndDelete = nodes[idna-4][jaa].score + scores.insertCodonDelete(codon);
	if( insAndDelete > c.score ) {
	  c.score = insAndDelete;
	  c.what = ins_read_delete;
	  c.loc = static_cast<CodonLocations>(m);
	}
      }
    }

    if( idna >= 2 ) {
      // nuc was deleted in frame. 3 possible ways
      int codon[3];
      for(uint m = dup_0; m <= dup_2 ; ++m) {
	getCodon(static_cast<CodonLocations>(m), idna, codon);  
	
	//  nuc was deleted in frame and codon matched to AA
	double const delAndMatch = nodes[idna-2][jaa - 1].score + scores.matchDuplicate(codon, curaa);
	if( delAndMatch > c.score ) {
	  c.score = delAndMatch;
	  c.what = match_duplicate;
	  c.loc = static_cast<CodonLocations>(m);
	}
	  
	// nuc was inserted to frame and codon deleted
	double const insAndDelete = nodes[idna-2][jaa].score + scores.insertCodonDelete(codon);
	if( insAndDelete > c.score ) {
	  c.score = insAndDelete;
	  c.what = ins_read_duplicate;
	  c.loc = static_cast<CodonLocations>(m);
	}
      }
    }
  }
}

AlignAndCorrect::Result
AlignAndCorrect::doAlignment(const byte* read, uint nRead, const byte* aa, uint naa)
{
  this->nRead = nRead;
  this->naa = naa;
  this->aa = aa;
  this->read = read;
      
  nodes = new Cell* [(nRead+1)];
  nodes[0] = new Cell[(nRead+1)*(naa+1)];
  for(uint k = 1; k <= nRead; ++k) {
    nodes[k] = nodes[k-1] + (naa+1);
  }

  for(uint i = 0; i <= nRead; i++) {
    for(uint j = 0; j <= naa; j++) {
      scoreCell(i, j);
    }
  }

  int icur = -1, jcur = -1;
  double maxval = std::numeric_limits<double>::lowest();
  for( uint i = 0; i <= nRead; i++ ) {
    if( maxval < nodes[i][naa].score ) {
      icur = i;
      maxval = nodes[i][naa].score;
    }
  }
  for( uint j = 0; j <= naa; j++ ) {
    if( maxval < nodes[nRead][j].score ) {
      jcur = j;
      maxval = nodes[nRead][j].score;
    }	
  }

  if( jcur != -1 ) {
    icur = nRead;
  } else {          
    jcur = naa;
  }
    
  Result res;

  res.dnaFreeEnd = nRead - icur;
  res.aaFreeEnd = naa - jcur;
    
  while( icur > 1 && jcur > 0 ) {
    Cell const& cur = nodes[icur][jcur];

    int codon[3];
    if( cur.what != ins_ref ) {
      getCodon(cur.loc, icur, codon);
    } else {
      codon[2] = codon[1] = codon[0] = gap;
      res.dnaBoundries.push_back(0);
    }
    res.alignedFramedRead.push_back(codon[2]);
    res.alignedFramedRead.push_back(codon[1]);
    res.alignedFramedRead.push_back(codon[0]);
    // jcur is decremented adn aa[jcur] needed later
    int const curaa = aa[jcur-1];
      
    bool aagap = false;
    int const icurbefore = icur;
    
    switch( cur.what ) {
      case match: {
	icur -= 3;
	jcur -= 1;
	break;
      }
      case ins_read_delete: {
	icur -= 4;
	aagap = true;
	break;
      }
      case match_delete: {
	icur -= 4;
	jcur -= 1;
	res.correctionDeletions += 1;
	break;
      }
      case ins_read_duplicate: {
	icur -= 2;
	aagap = true;
	break;
      }
      case match_duplicate: {
	icur -= 2;
	jcur -= 1;
	res.correctionInsertions += 1;
	break;
      }
      case ins_read: {
	icur -= 3;
	aagap = true;
	break;
      }
      case ins_ref: {
	res.gaps += 1;
	jcur -= 1;
	break;
      }
    }

    if( icur != icurbefore ) {
      res.dnaBoundries.push_back(icurbefore - icur);
    }
    
    if( aagap ) {
      res.gaps += 1;
      res.alignedAA.push_back( nAA + 1 );
    } else {
      int aas[16];
      uint const n = scores.getAAIndicesOfCodon(codon, aas);
      if( n == 0 ) {
	res.matches += 1;
      } else {
	res.mismatches += 1;
	for(uint k = 0; k < n; ++k) {
	  if( scores.translate(aas[k]) == curaa ) {
	    res.matches += 1;
	    res.mismatches -= 1;
	    break;
	  }
	}
      }
      res.alignedAA.push_back( curaa );
    }
  }
  res.aaFreeStart = jcur;
  res.dnaFreeStart = icur;
	
  delete [] nodes[0];
  delete [] nodes;

  return res;
}

PyObject*
aaCorrect(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char* kwlist[] = {"seq", "aaseq", "scoreMatrix", "geneticCode",
				 "indel", "correction", "stopCodon", 
				 static_cast<const char*>(0)};

  double indelPenalty = -10;
  double correctionPenalty = -10;
  double stopCodonPenalty = -100;

  PyObject* pScoreMatrix;
  PyObject* pGeneticCode;
  
  PyObject* pseq = 0;
  PyObject* paaseq = 0;

  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "OOOO|ddd", const_cast<char**>(kwlist),
				    &pseq,&paaseq,&pScoreMatrix,&pGeneticCode,
				    &indelPenalty,&correctionPenalty,&stopCodonPenalty)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! (PySequence_Check(pseq) && PySequence_Check(paaseq)) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not sequences") ;
    return 0;
  }

  if( ! (PySequence_Check(pScoreMatrix) && PySequence_Size(pScoreMatrix) == nAA * nAA) ) {
     PyErr_SetString(PyExc_ValueError, "wrong args: invalid AA score matrix");
     return 0;
  }

  if( ! (PySequence_Check(pGeneticCode) && PySequence_Size(pGeneticCode) == 64) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: invalid genetic code");
  }
  
  CorrectionScores<float> const scores(indelPenalty, correctionPenalty, stopCodonPenalty);

  int* scoreMatrix = new int[nAA * nAA];
  std::unique_ptr<int> dscoreMatrix(scoreMatrix);
  {
    PyObject* const s = PySequence_Fast(pScoreMatrix, "error");
    for(uint k = 0; k < nAA*nAA; ++k) {
      PyObject* const o = PySequence_Fast_GET_ITEM(s,k);
      int const m = PyInt_AS_LONG(o);
      scoreMatrix[k] = m;
    }
    Py_DECREF(s);
  }

  int geneticCode[64];
  {
    PyObject* const s = PySequence_Fast(pGeneticCode, "error");
    for(uint k = 0; k < 64; ++k) {
      PyObject* const o = PySequence_Fast_GET_ITEM(s,k);
      int const m = PyInt_AS_LONG(o);
      assert( -1 <= m && m < (int)nAA );
      geneticCode[k] = m;
    }
  }

  uint nseq=0, naa=0;
  byte* const dirtyRead = readSequence(pseq, nseq, true);
  if( ! dirtyRead ) {
    return 0;
  }
    
  byte* const peptide = readAASequence(paaseq, naa);

  if( ! peptide ) {
    return 0;
  }
  
  AlignAndCorrect ac(scoreMatrix, scores, geneticCode);

  AlignAndCorrect::Result const& res = ac.doAlignment(dirtyRead, nseq, peptide, naa);

  PyObject* tup = PyTuple_New(4);
  {
    auto const& aread = res.alignedFramedRead;
    int alen = aread.size();
    PyObject* ps0 = PyTuple_New(alen);
    alen -= 1;
    for(uint i = 0; alen >= 0; ++i,--alen) {
      PyTuple_SET_ITEM(ps0, i, PyInt_FromLong(aread[alen]));
    }
    PyTuple_SET_ITEM(tup, 0, ps0);
  }
  {
    auto const& aread = res.alignedAA;
    int alen = aread.size();
    PyObject* ps0 = PyTuple_New(alen);
    alen -= 1;
    for(uint i = 0; alen >= 0; ++i,--alen) {
      PyTuple_SET_ITEM(ps0, i, PyInt_FromLong(aread[alen]));
    }
    PyTuple_SET_ITEM(tup, 1, ps0);
  }
  {
    PyObject* d = PyDict_New();
    {
      PyObject* v = PyInt_FromLong(res.matches);
      PyDict_SetItemString(d, "matches", v);
      Py_DECREF(v);
    }
    {
      PyObject* v = PyInt_FromLong(res.mismatches);
      PyDict_SetItemString(d, "mismatches", v);
      Py_DECREF(v);
    }
    {
      PyObject* v = PyInt_FromLong(res.gaps);
      PyDict_SetItemString(d, "gaps", v);
      Py_DECREF(v);
    }
    {
      PyObject* v = PyInt_FromLong(res.correctionDeletions);
      PyDict_SetItemString(d, "correctionDeletions", v);
      Py_DECREF(v);
    }
    {
      PyObject* v = PyInt_FromLong(res.correctionInsertions);
      PyDict_SetItemString(d, "correctionInsertions", v);
      Py_DECREF(v);
    }
    {
      PyObject* v = PyInt_FromLong(res.dnaFreeEnd);
      PyDict_SetItemString(d, "dnaFreeEnd", v);
      Py_DECREF(v);
    }
    {
      PyObject* v = PyInt_FromLong(res.dnaFreeStart);
      PyDict_SetItemString(d, "dnaFreeStart", v);
      Py_DECREF(v);
    }
    {
      PyObject* v = PyInt_FromLong(res.aaFreeEnd);
      PyDict_SetItemString(d, "aaFreeEnd", v);
      Py_DECREF(v);
    }
    {
      PyObject* v = PyInt_FromLong(res.aaFreeStart);
      PyDict_SetItemString(d, "aaFreeStart", v);
      Py_DECREF(v);
    }
    PyTuple_SET_ITEM(tup, 2, d);
  }
  {
    auto const& b = res.dnaBoundries;
    int alen = b.size();
    PyObject* ps = PyTuple_New(alen);
    alen -= 1;
    for(uint i = 0; alen >= 0; ++i,--alen) {
      PyTuple_SET_ITEM(ps, i, PyInt_FromLong(b[alen]));
    }
    PyTuple_SET_ITEM(tup, 3, ps);
  }

  delete [] dirtyRead;
  delete [] peptide;
  
  return tup;
}

#include "seqslist.cc"

// template <typename T>
// inline T
// maxit(const T* x, uint start, uint end)
// {
//   T v = x[start];
//   start += 1;
  
//   while( start < end ) {
//     v = std::max(v, x[start]);
//     start += 1;
//   }
//   return v;
// }

inline float
maxit(const float* x, uint start, uint end)
{
  float v = x[start];
  start += 1;
  
  while( start < end ) {
    v = std::max(v, x[start]);
    start += 1;
  }
  return v;
}

// template <typename T>
// inline int
// imax(const T* x, uint start, uint end)
// {
//   int i = start;
//   T v = x[i];
//   start += 1;
  
//   while( start < end ) {
//     if( x[start] > v ) {
//       i = start;
//       v = x[i];
//     }
//     start += 1;
//   }
//   return i;
// }

inline int
imax(const float* x, uint start, uint end)
{
  int i = start;
  float v = x[i];
  start += 1;
  
  while( start < end ) {
    if( x[start] > v ) {
      i = start;
      v = x[i];
    }
    start += 1;
  }
  return i;
}

static inline int
decode(uint i, int (&c)[3])
{
  if( i < 64 ) {
    c[0] = i >> 4;
    i -= 16*c[0];
    c[1] = i >> 2;
    c[2] = i - 4*c[1];
    return 3;
  }
  if( i < 64 + 16 ) {
    i -= 64;
    c[0] = i >> 2;
    i -= 4*c[0];
    c[1] = i;
    return 2;
  }
  if( i < 64 + 16 + 4 ) {
    c[0] = i - 64 - 16;
    return 1;
  }
  return 0;
}

// inline byte
// lastnuc(uint i) {
//   return 0;
// }

static inline byte
lastnuc(uint i) {
  if( i > 64+16+4 ) {
    i -= 64+16+4;
  }
  if( i < 64 ) {
    return i & 0x3;
  }
  if( i < 64 + 16 ) {
    return (i-64) & 0x3;
  }
  return i - (64 + 16);
}

static inline bool
compat(uint const i, uint const ip)
{
  int ic[3], ipc[3];
  uint ni = decode(i, ic);
  uint nip = decode(ip, ipc);
  if( nip == 3 && ni == 1 ) {
    return true;
  }
  if( ni == nip+1 ) {
    for(uint k = 0; k < nip; ++k ) {
      if( ic[k] != ipc[k] ) {
	return false;
      }
    }
    return true;
  }
  
  return false;
}

static inline void
sortv(const float* const v, uint const n, int* inds)
{
  for(uint k = 0; k < n; ++k) {
    inds[k] = k;
  }
  auto fcmp = [v] (int i1, int i2) -> bool { return v[i1] > v[i2]; };
  std::sort(inds, inds + n, fcmp);
}

inline float
cost(byte const x, byte const y, const byte* const nb)
{
  // cost '-' -> '-' (deleting) 0
  // cost of '-' -> X (insert) mismatch (should depend on neighbors?)
  // cost of 'X' -> '-' (deleting) depends on neighbors
  // cost X -> X 0 (or positive?)
  // cost X -> Y (change) mismatch (should depend on neighbors?)
  if( x == y ) {
    return 0;
  }
  
  if( y == gap ) {
    uint const n = (x==nb[0]) + (x==nb[1]);
    return ((float[]){-6,-3,-2})[n];
  }
  
  // mismatch score
  return -5;
}

static byte*
getAAcons(SeqsList const& seqs, int const (&geneticCode)[64], uint& fstart)
{
  uint const nSeq = seqs.nSeqs;
  uint const nSites = seqs.seqslen[0];
  // 0 is left, 1 is right
  byte* const nb = new byte [nSites*2];

  uint const dsb = (64 + 16 + 4);
  uint const nv = 2*dsb;

  float** const seqScore = new float* [nSites];
  float** const profileScore = new float* [nSites];
  uint const nCells = nSites * nv;
  seqScore[0] = new float [nCells];
  profileScore[0] = new float [nCells];  // zeros
  std::fill(profileScore[0], profileScore[0]+nCells, 0);
  
  for(uint k = 1; k < nSites; ++k) {
    seqScore[k] = seqScore[k-1] + nv;
    profileScore[k] = profileScore[k-1] + nv;
  }

  for(uint ns = 0; ns < nSeq; ++ns) {
    const byte* seq = seqs.seqs[ns];
    // Establish neighbors
    nb[2*0 + 0] = gap;
    for(uint si = 1; si < nSites; ++si) {
      nb[2*si+0] = (seq[si-1] != gap) ? seq[si-1] : nb[2*(si-1)+0];
    }

    nb[2*(nSites-1) + 1] = gap;
    for(uint si = nSites-1; si > 0; --si) {
      nb[2*(si-1) + 1] = (seq[si] != gap) ? seq[si] : nb[2*si + 1];
    }

    {
      uint const si = 0;
      byte const xi = seq[si];
      const byte* const nbsi = nb + 2*si;
      float const c[4] = {cost(xi, 0, nbsi),cost(xi, 1, nbsi),cost(xi, 2, nbsi),cost(xi, 3, nbsi)};

      for(uint i = 0; i < 4; ++i) {
	seqScore[si][64+16 + i] = c[i]; 
	for(uint j = 0; j < 4; ++j) {
	  seqScore[si][64 + 4*j + i] = c[i]/2;
	  for(uint k = 0; k < 4; ++k) {
	    uint const ii = 16*k + 4*j + i;
	    // assume any wildcard is not all stop codons???
	    seqScore[si][ii] = geneticCode[ii] >= 0 ? c[i]/3 : -50000;
	  }
	}
      }
      float const dc = cost(xi, gap, nbsi);
      for(uint i = 0; i < dsb; ++i) {
	seqScore[si][dsb+i] = dc;
      }
    }
    for(uint si = 1; si < nSites; ++si) {
      byte const xi = seq[si];
      const byte* const nbsi = nb + 2*si;
      float const c[4] = {cost(xi, 0, nbsi),cost(xi, 1, nbsi),cost(xi, 2, nbsi),cost(xi, 3, nbsi)};
      float anyEndU = maxit(seqScore[si-1],0,64);
      float anyEndD = maxit(seqScore[si-1], dsb, dsb+64);
      float anyEnd = std::max(anyEndU,anyEndD);

      for(uint i = 0; i < 4; ++i) {
	seqScore[si][64+16 + i] = anyEnd + c[i]; 
	for(uint j = 0; j < 4; ++j) {
	  float const m0 = std::max(seqScore[si-1][64+16+j],seqScore[si-1][dsb + (64+16+j)]);
	  seqScore[si][64 + 4*j + i] = m0 + c[i]/2;
	  
	  for(uint k = 0; k < 4; ++k) {
	    uint const ii = 16*k + 4*j + i;
	    // assume any wildcard is not all stop codons???
	    if( geneticCode[ii] >= 0 ) {
	      float const m1 = std::max(seqScore[si-1][64 + 4*k + j], seqScore[si-1][dsb + 64 + 4*k + j]);
	      seqScore[si][ii] = m1 + c[i]/3;
	    } else {
	      seqScore[si][ii] = -50000;
	    }
	  }
	}
      }
      float const dc = cost(xi, gap, nbsi);
      for(uint i = 0; i < dsb; ++i) {
	seqScore[si][dsb+i] = std::max(seqScore[si-1][i], seqScore[si-1][dsb + i]) + dc;
      }
    }
    float* s = seqScore[0];
    float* p = profileScore[0];
    for(uint i = 0; i < nCells; ++i) {
      *p++ += *s++;
    }
  }

  byte* states = new byte [nSites];
  int* istates = new int [nSites]();
  std::fill(istates, istates+nSites, -1);
  
  int si = nSites-1;
  while( imax(profileScore[si],0,nv) >= int(dsb) ) {
    states[si] = gap;
    si -= 1;
  }

  uint i = imax(profileScore[si],0,nv);
  states[si] = lastnuc(i);
  istates[si] = i;
  uint isi = i;
  
  while( si > 0 ) {
    si -= 1;
    uint i1 = imax(profileScore[si],0,nv);
    if( i1 < dsb ) {
      if( compat(i, i1) ) {
        states[si] = lastnuc(i1);
	istates[si] = i1;
        
        i = i1 ;
	isi = si;
      } else {
	//print "1",si,
	int sinds[nv];
	sortv(profileScore[si], nv, sinds);
	// i1 might not be top if there are equals to it
	for(uint k = 0; k < nv; ++k) {
	  i1 = sinds[k];
	  if( i1 < dsb ) {
	    if( compat(i, i1) ) {
              states[si] = lastnuc(i1);
	      istates[si] = i1;

              i = i1 ; isi = si;
              break;
	    }
	  } else {
	    if( compat(i, i1-dsb) ) {
              states[si] = gap;
              break;
	    }
	  }
	}
      }
    } else {
      // deleted state for site
      states[si] = gap;
      if( ! compat(i, i1-dsb) ) {
	// deleted state *is not* compatible with current
	int sinds[nv];
	sortv(profileScore[si], nv, sinds);
	// i1 might not be top if there are equals to it
	for(uint k = 0; k < 11; ++k) {
	  i1 = sinds[k];
	  if( i1 >= dsb ) {
	    if( compat(i, i1-dsb) ) {
	      // a deleted state compatible with current, keep deleted status
	      break;
	    }
	  } else {
	    if( compat(i, i1) ) {
	      // a non deleted state compatible with current, cancel gap
	      states[si] = lastnuc(i1)/*[-1]*/ ;
	      istates[si] = i1;

	      i = i1 ; isi = si;
	      break;
	    }
	  }
	}
      }
    }
  }

  i = istates[isi];                assert( i < dsb );
  int t[3];
  int const frame = decode(i,t) - 1;
  fstart = isi + ((3-frame) % 3);
  
  delete [] profileScore[0];
  delete [] seqScore[0];
  delete [] seqScore;
  delete [] profileScore;

  delete [] nb;

  delete [] istates;
  
  return states;
}
	  

PyObject*
aaCons(PyObject*, PyObject* args, PyObject* kwds)
{
  static const char* kwlist[] = {"seqs", "geneticCode",
				 static_cast<const char*>(0)};

  PyObject* pGeneticCode;
  PyObject* pseqs = 0;
  
  if( ! PyArg_ParseTupleAndKeywords(args, kwds, "OO", const_cast<char**>(kwlist),
				    &pseqs,&pGeneticCode) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! PySequence_Check(pseqs) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: not sequences") ;
    return 0;
  }

  if( ! (PySequence_Check(pGeneticCode) && PySequence_Size(pGeneticCode) == 64) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args: invalid genetic code");
  }
  
  int geneticCode[64];
  {
    PyObject* const s = PySequence_Fast(pGeneticCode, "error");
    for(uint k = 0; k < 64; ++k) {
      PyObject* const o = PySequence_Fast_GET_ITEM(s,k);
      int const m = PyInt_AS_LONG(o);
      assert( -1 <= m && m < (int)nAA );
      geneticCode[k] = m;
    }
  }

  std::unique_ptr<const SeqsList> seqs(readSeqsIn(pseqs, false));
  if( ! seqs || seqs->nSeqs == 0 ) {
    return 0;
  }
  uint fstart;
  const byte* s = getAAcons(*seqs, geneticCode, fstart);

  PyObject* tup = PyTuple_New(2);
  
  uint const nSites = seqs->seqslen[0];
  PyObject* ps = PyTuple_New(nSites);
  for(uint i = 0; i < nSites; ++i) {
    PyTuple_SET_ITEM(ps, i, PyInt_FromLong(s[i]));
  }

  PyTuple_SET_ITEM(tup, 0, ps);
  PyTuple_SET_ITEM(tup, 1, PyInt_FromLong(fstart));
  
  delete [] s;
  
  return tup;
  
  // PyObject* ret = Py_None;
  // Py_INCREF(ret);
  
  // return ret;
}

  
PyDoc_STRVAR(aalign__doc__,
"AA Sequence alignment");


static PyMethodDef aalignMethods[] = {
  {"acorrect",	(PyCFunction)aaCorrect, METH_VARARGS|METH_KEYWORDS,
   "Align nucs to AA reference and correct errors."},
  
  {"aacons",	(PyCFunction)aaCons, METH_VARARGS|METH_KEYWORDS,
   "Valid coding amino acid from nuclieotide alignment."},
  
  {NULL, NULL, 0, NULL}        /* Sentinel */
};



PyMODINIT_FUNC
initaalign(void)
{
  PyObject* const m = Py_InitModule3("aalign", aalignMethods, aalign__doc__);

  PyObject* const o = PyString_FromString(aminoAcidsOrder);
  PyObject_SetAttrString(m, const_cast<char*>("AAorder"), o);
  Py_DECREF(o);
}
