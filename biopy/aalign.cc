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

static const char aminoAcidsOrder[24] = "ABCDEFGHIKLMNPQRSTVWXYZ";
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
    c.what = /*some kind of none;*/ (FType)-1;
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
    
  while( icur > 0 && jcur > 0 ) {
    Cell const& cur = nodes[icur][jcur];

    int codon[3];
    if( cur.what != ins_ref ) {
      getCodon(cur.loc, icur, codon);
    } else {
      codon[2] = codon[1] = codon[0] = 5; /* damn */
    }
    res.alignedFramedRead.push_back(codon[2]);
    res.alignedFramedRead.push_back(codon[1]);
    res.alignedFramedRead.push_back(codon[0]);
    // jcur is decremented adn aa[jcur] needed later
    int const curaa = aa[jcur-1];
      
    bool aagap = false;
      
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
  byte* const peptide = readAASequence(paaseq, naa);
  
  AlignAndCorrect ac(scoreMatrix, scores, geneticCode);

  AlignAndCorrect::Result const& res = ac.doAlignment(dirtyRead, nseq, peptide, naa);

  PyObject* tup = PyTuple_New(3);
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

  delete [] dirtyRead;
  delete [] peptide;
  
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
