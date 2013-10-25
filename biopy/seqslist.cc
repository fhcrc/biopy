// Ugly to include this file, but again, tell me wise old man,
// how do I share this code between diffrent c++ modules???

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

static SeqsList*
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
