#include <Python.h>
#include "structmember.h"

#include <cassert>

#include <cmath>
#include <limits>

#include <string>
using std::string;
#include <algorithm>
#include <vector>
using std::vector;
#include <stack>
using std::stack;
#include <list>
using std::list;

// for compilers lacking it
typedef unsigned int uint;

#include <map>
// change to 0 for compilers not-supporting new hash maps
#if 1
#include <unordered_map>
using std::unordered_map;
#else
#define unordered_map map
using std::map;
#endif

template<typename T>
PyObject*
dvector2tuple(vector<T> const& v)
{
  if( &v ) {
    PyObject* t = PyTuple_New(v.size());
    for(uint k = 0; k < v.size(); ++k) {
      PyTuple_SET_ITEM(t, k, PyFloat_FromDouble(v[k]));
    }
    return t;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

template<typename T>
PyObject*
ivector2tuple(vector<T> const& v)
{
  if( &v ) {
    PyObject* t = PyTuple_New(v.size());
    for(uint k = 0; k < v.size(); ++k) {
      PyTuple_SET_ITEM(t, k, PyInt_FromLong(v[k]));
    }
    return t;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static inline bool
areSame(double a, double b)
{
  return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

static inline int
lg2i(int v)
{
  int r = 0;
  v >>= 1;
  while( v ) {
    r += 1;
    v >>= 1;
  }
  return r;
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

#undef PYINTERFACE

#ifdef PYINTERFACE
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
#else
static inline string
trimString(string const& s)
{
  const char* txt = s.c_str();
  int const n0 = skipSpaces(txt);
  int n1 = s.length()-1;
  while( n1 >= 0 && isspace(txt[n1]) ) {
    --n1;
  }
  return string(txt + n0, n1+1-n0);
}
#endif


#ifndef PYINTERFACE
// list
typedef vector< std::pair<string,string> > Attributes;

static PyObject*
attributesAsPyObj(const Attributes* const attributes)
{
  if( attributes ) {
    Attributes const& atr = *attributes;
    PyObject* a = PyDict_New();
    for(uint k = 0; k < atr.size(); ++k) {
      //auto const p = atr[k];
      std::pair<string,string> const p = atr[k];
      PyObject* const val = PyString_FromString(p.second.c_str());
      PyDict_SetItemString(a, p.first.c_str(), val);
    }
    return a;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

class ParsedTreeNode {
public:
  ParsedTreeNode() :
    branch(0),
    attributes(0)
    {}

  ~ParsedTreeNode() {
    delete branch;
    delete attributes;
  }

  ParsedTreeNode(ParsedTreeNode const& o) :
    taxon(o.taxon),
    sons(o.sons)
    {
      // steal memebers from initializer. make sure he can't
      // access them anymore.
      branch = o.branch;
      attributes = o.attributes;
      const_cast<ParsedTreeNode&>(o).branch = 0;
      const_cast<ParsedTreeNode&>(o).attributes = 0;
    }

  PyObject* asPyObject(void) const;

  string 		taxon;
  double* 		branch;
  std::vector<uint> 	sons;
  Attributes* 		attributes;
};

PyObject*
ParsedTreeNode::asPyObject(void) const
{
  PyObject* const nodeData = PyList_New(4);

  if( taxon.size() > 0 ) {
    PyList_SET_ITEM(nodeData, 0, PyString_FromStringAndSize(taxon.c_str(),taxon.size()));
  } else {
    Py_INCREF(Py_None);
    PyList_SET_ITEM(nodeData, 0, Py_None);
  }
    
  PyObject* b;
  if( branch ) {
    b = PyFloat_FromDouble(*branch);
  } else {
    Py_INCREF(Py_None);
    b = Py_None;
  }
  PyList_SET_ITEM(nodeData, 1, b);

  int const ns = sons.size();
  if( ns > 0 ) {
    PyObject* psubs = PyList_New(ns);
    for(int k = 0; k < ns; ++k) { 
      PyList_SET_ITEM(psubs, k, PyLong_FromLong(sons[k]));
    }
    PyList_SET_ITEM(nodeData, 2, psubs);
  } else {
    Py_INCREF(Py_None);
    PyList_SET_ITEM(nodeData, 2, Py_None);
  }
  
  PyList_SET_ITEM(nodeData, 3, attributesAsPyObj(attributes));

  return nodeData;
}

#endif

typedef
#if PYINTERFACE
PyObject*
#else
ParsedTreeNode
#endif
returnType;

static int
parseAttributes(const char* s,
#if PYINTERFACE
vector<PyObject*>&
#else
Attributes&
#endif
		vals)
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

#if PYINTERFACE
    PyObject* o = PyTuple_New(2);
    PyTuple_SET_ITEM(o, 0, trimString(name));
    PyTuple_SET_ITEM(o, 1, trimString(v));
    vals.push_back(o);
#else
    std::pair<string,string> a(trimString(name),trimString(v));
    vals.push_back(a);
#endif
  }
  return eat;
}


static int
readSubTree(const char* txt, vector<returnType>& nodes)
{
  int eat = skipSpaces(txt);
  txt += eat;
  
#if PYINTERFACE
  PyObject* vals = NULL;
  PyObject* const nodeData = PyList_New(4);
#else
  ParsedTreeNode nodeData;
#endif
  
  if( *txt == '(' ) {
#if PYINTERFACE
    vector<int> subs;
#else
    vector<uint>& subs = nodeData.sons;
#endif
    
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
#ifdef PYINTERFACE	
	int const ns = subs.size();
	PyObject* psubs = PyList_New(ns);
	for(int k = 0; k < ns; ++k) { 
	  PyList_SET_ITEM(psubs, k, PyLong_FromLong(subs[k]));
	}
	PyList_SET_ITEM(nodeData, 2, psubs);

	Py_INCREF(Py_None);
	PyList_SET_ITEM(nodeData, 0, Py_None);
#endif
	
	eat += 1;
        txt += 1;
        break;
      }
      return -1; // error
    }
  } else {
    const char* s = txt;
    assert( ! isspace(*s) ) ;
    // a terminal
    if( *s == '\'' || *s == '"' ) {
      int const e = _getStuff(s+1, *s);
      if( e < 0 ) {
	return -1;
      }
      s += e+2;
    } else {
      while( ! isspace(*s) && *s != ':' && *s != '[' && *s != ','
	     && *s != '(' && *s != ')' && *s != ']' ) {
	++s;
      }
    }
    int const n1 = s - txt;
#ifdef PYINTERFACE
    PyList_SET_ITEM(nodeData, 0, PyString_FromStringAndSize(txt, n1));

    Py_INCREF(Py_None);
    PyList_SET_ITEM(nodeData, 2, Py_None);
#else
    nodeData.taxon = string(txt, n1);
#endif
    
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
    if( has(*txt, "(),;") ) {
      break;
    }
    if( *txt == '[' ) {
      if( txt[1] == '&' ) {
#ifdef PYINTERFACE
	vector<PyObject*> vs;
#else
	nodeData.attributes = new Attributes;
	Attributes& vs = *nodeData.attributes;
#endif	
	int n1 = parseAttributes(txt+2, vs);
	if( n1 < 0 ) {
	  return -1;
	}
	n1 += 3;
	n1 += skipSpaces(txt+n1);
	eat += n1;
	txt += n1;

#ifdef PYINTERFACE
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
#endif
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
      //if( isspace(*txt) || isdigit(*txt) || has(*txt, ":.+-Ee") ) {
      nodeTxt.append(txt, 1);
      txt += 1;
      eat += 1;
      // } else {
      // 	break;
      // }
    }
  }

#if PYINTERFACE
  if( ! vals ) {
    Py_INCREF(Py_None);
    vals = Py_None;
  }
  PyList_SET_ITEM(nodeData, 3, vals);

  PyObject* branch = NULL;
#endif

  //std::cout << nodeTxt << std::endl;  
  const char* nTxt = nodeTxt.c_str();
  nTxt += skipSpaces(nTxt);
  if( *nTxt ) {
    int const i = findIndex(nTxt, ':',"");
    if( i != 0 ) {
      int const k = (i>0) ? i : strlen(nTxt);
#if PYINTERFACE
      PyList_SET_ITEM(nodeData, 0, trimString(string(nTxt,k)));
#else
      nodeData.taxon = string(nTxt,k);
#endif
      nTxt += k;
    }
    if( i >= 0 ) {
      int n1 = skipSpaces(nTxt+1);
      nTxt += 1+n1;
    
      //std::cout << nTxt << std::endl;  
      char* endp;
      double const b = strtod(nTxt, &endp);
      n1 = endp - nTxt;
      if( n1 == 0 ) {
	return -1;
      }
#if PYINTERFACE    
      branch = PyFloat_FromDouble(b);
#else
      nodeData.branch = new double(b);
#endif
    }
  }
#if PYINTERFACE
  if( branch == NULL ) {
    Py_INCREF(Py_None);
    branch = Py_None;
  }
  PyList_SET_ITEM(nodeData, 1, branch);
#endif
  nodes.push_back(nodeData);

  return eat;
}

template<typename T> class Packer {
public:
  // Number of stored values. 
  virtual uint size(void) const = 0;
  
  // Stored values as a vector. Valid until next call to unpacked
  // from any instantiation
  virtual vector<T> const&  unpacked(void) const = 0;
  
  // Stored values as a vector. Sets isPermanent to false if vector is transient
  // (becomes invalid on the next call to unpacked from anywhere) and to true if
  // the return value is safe to keep around.
  virtual vector<T> const&  unpacked(bool& isPermanent) const = 0;
};

template<typename T>
class SimplePacker : public Packer<T> {
public:
  SimplePacker<T>(vector<T> _vals) :
    vals(_vals)
    {}

  template<typename U> SimplePacker<T>(vector<U> _vals) {
    for(auto x = _vals.begin(); x != _vals.end(); ++x) {
      vals.push_back(*x);
    }
  }

  virtual uint size(void) const { return vals.size(); }
  
  vector<T> const&  unpacked(void) const { return vals; }
  vector<T> const&  unpacked(bool& isPermanent) const {
    isPermanent = true; return vals;
  }
  
private:
  vector<T> vals;
};


class FixedIntPacker : public Packer<uint> {
public:
  FixedIntPacker(uint nBitsPerValue, vector<uint>::const_iterator from, vector<uint>::const_iterator to);
  virtual ~FixedIntPacker();

  vector<uint> const&  unpacked(void) const;
  vector<uint> const&  unpacked(bool& isPermanent) const;

  virtual uint size(void) const { return len; }

  uint const nBitsPerValue : 8;
  uint const len : 24;
private:
  static uint const usize = 8*sizeof(char);
  static vector<uint> temp;
  
  unsigned char* bits;
};

FixedIntPacker::~FixedIntPacker() {
  delete [] bits;
}

vector<uint> FixedIntPacker::temp;

template<typename T>
inline T lowerNbits(uint n) {
  return (static_cast<T>(1) << n) - 1;
}

template<typename T>
inline T upperNbits(uint n) {
  return ~lowerNbits<T>(8*sizeof(T)-n);
}

FixedIntPacker::FixedIntPacker(uint _nBitsPerValue,
			       vector<uint>::const_iterator from,
			       vector<uint>::const_iterator to) :
  nBitsPerValue(_nBitsPerValue),
  len(to-from)
{
  int const sbits = (nBitsPerValue * len + (usize-1)) / usize;
  bits = new unsigned char [sbits];

  assert( nBitsPerValue <= usize );
  if( nBitsPerValue == usize ) {
    uint k = 0;
    for(auto v = from; v < to; ++v, ++k) {
      bits[k] = *v;
    }
  } else {
    
    std::fill(bits, bits+sbits,0);
  
    unsigned char* cur = bits;
    uint loc = 0;
  
    for(auto v = from; v < to; ++v)  {
      loc += nBitsPerValue;
      if( loc <= usize ) {
	*cur |= (*v << (usize - loc));
      } else {
	loc -= usize; // now loc is number of bits in next block
	*cur |= (*v >> loc);
	++cur;                          assert( cur - bits < sbits );
	uint const s = usize - loc;
	*cur |= (*v << s) & ~((1 << s)-1);
      }
    }
  }
}

vector<uint> const&
FixedIntPacker::unpacked(void) const
{
  temp.clear();
  temp.resize(len, 0);
  unsigned char* cur = bits;

  if( nBitsPerValue == usize ) {
    for(int k = 0; k < len; ++k) {
      temp[k] = bits[k];
    }
    return temp;
  }
  
  int loc = usize;
  unsigned char mask = lowerNbits<unsigned char>(nBitsPerValue);
  
  for(uint k = 0; k < len; ++k) {
    loc -= nBitsPerValue;
    if( loc >= 0 ) {
      temp[k] = ((*cur) >> loc) & mask;
    } else {
      uint const upper = *cur & lowerNbits<unsigned char>(loc + nBitsPerValue);
      uint const left = -loc;

      ++cur;
      loc = usize - left;
      uint const lower = (*cur & upperNbits<unsigned char>(left)) >> loc;
      temp[k] = (upper << left) | lower;
    }
  }
  return temp;
}

vector<uint> const&
FixedIntPacker::unpacked(bool& isPermanent) const
{
  isPermanent = false;
  return unpacked();
}

class TreeRep {
public:
  // steals attributes
  TreeRep(Packer<uint>& t, Packer<uint>* l, vector<const Attributes*>* atrbs);
  virtual ~TreeRep();

  virtual bool isCladogram(void) const = 0;
  
  uint nTaxa(void) const { return ptips.size(); }
  
  vector<uint> const& tips() const;
  vector<uint> const& tips(bool& permanent) const;

  vector<uint>* labels(void) const;
  
  const vector<const Attributes*>* getAttributes(void) const {
    return attributes;
  }
  
protected:
  Packer<uint>&	ptips;
  Packer<uint>*	plabels;
  vector<const Attributes*>* attributes;
};

inline
TreeRep::TreeRep(Packer<uint>& t, Packer<uint>* l, vector<const Attributes*>* atrbs) :
  ptips(t),
  plabels(l),
  attributes(atrbs)
{}

TreeRep::~TreeRep()
{
  delete &ptips;
  if( plabels ) {
    delete plabels;
  }
  if( attributes ) {
    for(auto i = attributes->begin(); i < attributes->end(); ++i) {
      delete *i;
    }
    delete attributes;
  }
}

vector<uint> const&
TreeRep::tips(void) const
{
  return ptips.unpacked();
}

vector<uint> const&
TreeRep::tips(bool& permanent) const
{ 
  return ptips.unpacked(permanent);
}

vector<uint>*
TreeRep::labels(void) const
{
  if( ! plabels ) {
    return 0;
  }
  vector<uint> const& l = plabels->unpacked();
  return new vector<uint>(l);
}

class CladogramRep : public TreeRep {
public:
  // Steals atrbs
  CladogramRep(Packer<uint>& t, Packer<uint>* l,
	       Packer<uint>* h, vector<const Attributes*>* atrbs);
  virtual ~CladogramRep();
 
  virtual bool isCladogram(void) const { return true; }
 
  vector<uint> const& heights() const { return pheights->unpacked(); }

private:
  Packer<uint>*   pheights;
};

inline
CladogramRep::CladogramRep(Packer<uint>& t, Packer<uint>* l, Packer<uint>* h,
			   vector<const Attributes*>* atrbs) :
  TreeRep(t, l, atrbs),
  pheights(h)
{}

CladogramRep::~CladogramRep()
{
  delete pheights;
}

template<typename T>
class PhylogramRep : public TreeRep {
public:
  // Steals atrbs
  PhylogramRep(Packer<uint>& t, Packer<uint>* l, Packer<T>* h,
	       Packer<T>* txh, vector<const Attributes*>* atrbs);
  virtual ~PhylogramRep();

  virtual bool isCladogram(void) const { return false; }
  
  vector<T> const& heights() const { return pheights->unpacked(); } 
  vector<T> const* txheights() const {
    return ptxheights ? &ptxheights->unpacked() : static_cast< vector<T>* >(0);
  } 

private:
  // Packed internal nodes heights 
  Packer<T>*	pheights;
  // Packed taxa heights (if any - null means all zero.)
  Packer<T>* 	ptxheights;
};

template<typename T>
PhylogramRep<T>::PhylogramRep(Packer<uint>& t, Packer<uint>* l, Packer<T>* h,
			      Packer<T>* txh, vector<const Attributes*>* atrbs) :
  TreeRep(t, l, atrbs),
  pheights(h),
  ptxheights(txh)
{}

template<typename T>
PhylogramRep<T>::~PhylogramRep() {
  delete pheights;
  delete ptxheights;
}

// Tree should have been a nested class of Trees set
class TreesSet;

class Tree {
public:
  Tree(TreesSet const& _ts, uint _nt);
  ~Tree();
  
  vector<uint> const& tips() const;
  
  void getTerminals(vector<uint>& terms) const;
  
  uint getRootID(void) const {
    if( ! internals ) setup();
    return internals->size() - 1;
  }

  struct Expanded {
    Expanded(int	itax,
	     uint 	nSons,
	     uint* 	sons,
	     double* 	branch,
	     double 	height,
	     int    	prev,
	     const Attributes* attributes);
    
    int 	itax;
    uint 	nSons;
    uint* 	sons;
    double* 	branch;
    double 	height;
    int    	prev;
    const Attributes* attributes;
  };

  bool isCladogram(void) const;
  
  uint nNodes(void) const;
  
  Expanded const& getNode(uint n) const {
    if( ! internals ) setup();
    return (*internals)[n];
  }

  void toNewick(string& s, int nodeId, bool topoOnly, bool includeStem, bool withAttribute) const;

  TreesSet const& ts;
  uint const 	  nt;
private:
  void setup() const;
  
  void tostr(vector<string>& s, uint nodeId, bool topoOnly, bool includeStem, bool withAttributes) const;
  
  uint rep2treeInternal(vector<Expanded>&                 nodes,
			uint                              low,
			uint                              hi,
			vector<uint> const&               tax,
			vector<double> const&             htax,
			vector<double> const&             hs,
			const vector<const Attributes*>*  atrbs,
			const vector<uint>* const         labels,
			uint*&                            sonsBlock,
			uint*                             curiScratch,
			uint                              bleft) const;

  // only after setup  
  mutable vector<Expanded>* internals;

  // Number of taxa 
  mutable uint              nTaxa;
  
  // one contiguous block for all sons indices
  mutable uint*		    sonsBlockSave;
};

inline
Tree::Tree(TreesSet const& _ts, uint _nt) :
  ts(_ts),
  nt(_nt),
  internals(0),
  sonsBlockSave(0)
{}

Tree::~Tree() {
  delete internals;
  delete [] sonsBlockSave;
}


struct TreesSetObject;
struct TreeNodeObject;

struct TreeObject : PyObject {
  void del(void);

  void init(TreesSetObject* _ts = NULL, Tree* t = NULL);

  PyObject* getTaxa(void) const;

  PyObject* getTerminals(void) const ;
  PyObject* allIds(void) const;
  uint      getRootID(void) { return tr->getRootID(); }
  PyObject* getNode(uint nt) const;

  PyObject* toNewick(int nodeId, bool topoOnly, bool includeStem, bool withAttributes) const;
  void      getInOrder(bool preOrder, vector<int>& ids, int nodeId, bool includeTaxa);
  
  void      setBranch(uint nid, double branch);

  // public to be used in offsetof in type structure
  PyObject*	  dict;

private:
  void      tostr(vector<string>& s, int nodeId,
		  bool topoOnly, bool includeStem, bool withAttributes) const;
  
  Tree* 	  tr;
  TreesSetObject* ts;
  mutable PyObject*       		taxa;
  mutable vector<TreeNodeObject*>*	treeNodes;
};

//#include <iostream>

void
TreeObject::init(TreesSetObject* _ts, Tree* t)
{
  tr = t;
  ts = _ts;
  taxa = 0;
  treeNodes = 0;
}

void
TreeObject::del(void)
{
  delete tr;
  Py_XDECREF(taxa);
  if( treeNodes ) {
    for( auto n = treeNodes->begin(); n != treeNodes->end(); ++n ) {
      Py_XDECREF(*n);
    }
    Py_XDECREF(treeNodes);
  }
  Py_XDECREF(dict);
}

static void
Tree_dealloc(TreeObject* self)
{
  self->del();
  self->ob_type->tp_free((PyObject*)self);
}

static TreeObject*
Tree_new(PyTypeObject* type, PyObject*, PyObject*)
{
  TreeObject* self = (TreeObject *)type->tp_alloc(type, 0);

  if (self != NULL) {
    self->init();
    self->dict = PyDict_New();
  }

  return self;
}

PyObject*
Tree_getattr(TreeObject *self, PyObject* aname)
{
  if( PyString_Check(aname) ) {
    const char* const name = PyString_AS_STRING(aname);
    if( ! strcmp(name, "root") ) {
      return PyInt_FromLong(self->getRootID());
    }
  }
  return PyObject_GenericGetAttr(self, aname);
}

PyObject*
Tree_str(TreeObject *self)
{
  return self->toNewick(-1, false, false, false);
}

PyObject*
tree_getTaxa(TreeObject* self)
{
  return self->getTaxa();
}

PyObject*
tree_getTerminals(TreeObject* self)
{
  return self->getTerminals();
}

PyObject*
tree_allIds(TreeObject* self)
{
  return self->allIds();
}

PyObject*
tree_xorder(TreeObject* self, PyObject* args, PyObject* kwds, bool pre)
{
  static const char *kwlist[] = {"node", "includeTaxa", static_cast<const char*>(0)};
  int nodeId;
  PyObject* incTax = 0;
  if (! PyArg_ParseTupleAndKeywords(args, kwds, "i|O", (char**)kwlist,
				    &nodeId,&incTax)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.");
    return 0;
  }
  bool const incT = (!incTax || PyObject_IsTrue(incTax));
  vector<int> ids;
  self->getInOrder(pre, ids, nodeId, incT);

  PyObject* r = ivector2tuple(ids);
  return r;
}

PyObject*
tree_postorder(TreeObject* self, PyObject* args, PyObject* kwds)
{
  return tree_xorder(self, args, kwds, false);
}

PyObject*
tree_preorder(TreeObject* self, PyObject* args, PyObject* kwds)
{
  return tree_xorder(self, args, kwds, true);
}

PyObject*
tree_node(TreeObject* self, PyObject* args)
{
  int n;
  if( !PyArg_ParseTuple(args, "i", &n) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.");
    return 0;
  }
  return self->getNode(n);
}

PyObject*
tree_setBranch(TreeObject* self, PyObject* args)
{
  int n;
  double br;
  if( !PyArg_ParseTuple(args, "id", &n, &br) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.");
    return 0;
  }
  if( br < 0 ) {
     PyErr_SetString(PyExc_ValueError, "negative branch length.");
     return 0;
  }
  self->setBranch(n,br);
  Py_RETURN_NONE;
}

// PyObject*
// tree_setHeight(TreeObject* self, PyObject* args)
// {
//   int n;
//   double h;
//   if( !PyArg_ParseTuple(args, "id", &n,&h) ) {
//     PyErr_SetString(PyExc_ValueError, "wrong args.");
//     return 0;
//   }
//   self->setHeight(n,h);
//   Py_RETURN_NONE;
// }


PyObject*
tree_2newick(TreeObject* self, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"nodeId", "topologyOnly", "attributes", "includeStem",
				 static_cast<const char*>(0)};
  PyObject* attr = 0;
  PyObject* topo = 0;
  PyObject* incStem = 0;
  int nodeId = -1;
  
  if (! PyArg_ParseTupleAndKeywords(args, kwds, "|iOOO", (char**)kwlist,
				    &nodeId,&topo,&attr,&incStem)) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  // 
  bool const topoOnly = (topo && PyObject_IsTrue(topo));
  bool const incS = (incStem && PyObject_IsTrue(incStem));
  bool const attrs = (attr && PyObject_IsTrue(attr));
  
  return self->toNewick(nodeId, topoOnly, incS, attrs);
}

#if 0
static int
Tree_init(TreeObject* self, PyObject* args, PyObject* kwds)
{
  // static const char *kwlist[] = {"compressed", "precision", static_cast<const char*>(0)};
  // PyObject* comp = 0;
  // int precision = 4;
  
  // if (! PyArg_ParseTupleAndKeywords(args, kwds, "|Oi", (char**)kwlist,
  // 				    &comp,&precision)) {
  //   return -1;
  // }

  // if( ! (precision == 4 || precision == 8) ) {
  //   PyErr_SetString(PyExc_ValueError, "wrong args (precision)");
  //   return -1;
  // }
    
  // bool const compressed = (! comp || PyObject_IsTrue(comp));

  // self->ts = new TreesSet(compressed, precision);
  return 0;
}
#endif

static PyMethodDef tree_methods[] = {
  {"get_taxa", (PyCFunction)tree_getTaxa, METH_NOARGS,
   "get tree taxa."
  },
  {"get_terminals", (PyCFunction)tree_getTerminals, METH_NOARGS,
   "get tree terminals ids."
  },
  {"all_ids", (PyCFunction)tree_allIds, METH_NOARGS,
   "get all tree ids."
  },
  {"in_postorder", (PyCFunction)tree_postorder, METH_VARARGS|METH_KEYWORDS,
   "get sub-tree in post-order."
  },
  {"in_preorder", (PyCFunction)tree_preorder, METH_VARARGS|METH_KEYWORDS,
   "get sub-tree in post-order."
  },
  {"node", (PyCFunction)tree_node, METH_VARARGS,
   "get tree node."
  },
  {"setBranch", (PyCFunction)tree_setBranch, METH_VARARGS,
   "Set branch length of node 'n'"
  },
  {"toNewick", (PyCFunction)tree_2newick, METH_VARARGS|METH_KEYWORDS,
   "tree as string in NEWICK."
  },
  
  {NULL}  /* Sentinel */
};

static PyTypeObject TreeType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "treesset.Tree",        /*tp_name*/
    sizeof(TreeObject),    /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)Tree_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash*/
    0,                         /*tp_call*/
    (reprfunc)Tree_str,        /*tp_str*/
    (getattrofunc)Tree_getattr, /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Tree object",             /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    tree_methods,             /* tp_methods */
    0/*tree_members*/,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    offsetof(TreeObject,dict),  /* tp_dictoffset */
    0/*(initproc)Tree_init*/,      /* tp_init */
    0,                         /* tp_alloc */
    (newfunc)Tree_new,         /* tp_new */
};



struct TreeNodeDataObject : PyObject {
  bool hasBranch(void) const {
    return branchlength != Py_None;
  }
  
  double getBranch(void) const {
    return hasBranch() ? PyFloat_AsDouble(branchlength) : 0;
  }

  // return branch length before setting
  double setBranch(double newLen);
  
  bool hasHeight(void) const {
    return height != Py_None;
  }
  
  double adjustHeight(double dif);
  
  PyObject* taxon;
  PyObject* branchlength;
  PyObject* height;
  PyObject* allData;
};

double
TreeNodeDataObject::setBranch(double newLen)
{
  double const b = getBranch();
  if( hasBranch() ) {
    Py_DECREF(branchlength);
    branchlength = PyFloat_FromDouble(newLen);
  } else {
    branchlength = PyFloat_FromDouble(newLen);
  }
  return b;
}

double
TreeNodeDataObject::adjustHeight(double dif)
{
  if( ! hasHeight() ) {
    height = PyFloat_FromDouble(dif);
    return dif;
  }
    
  double const h = PyFloat_AsDouble(height) + dif;
  Py_DECREF(height);
  height = PyFloat_FromDouble(h);
  return h;
}


static PyMemberDef NodeData_members[] = {
  {(char*)"taxon", T_OBJECT_EX, offsetof(TreeNodeDataObject, taxon), 0,
    (char*)"taxon"},
  {(char*)"branchlength", T_OBJECT_EX, offsetof(TreeNodeDataObject, branchlength), 0,
    (char*)"branchlength"},
  {(char*)"height", T_OBJECT_EX, offsetof(TreeNodeDataObject, height), 0,
    (char*)"height"},
  {NULL}  /* Sentinel */
};

static TreeNodeDataObject *
TreeNodeData_new(PyTypeObject* type, PyObject *args, PyObject *kwds)
{
  TreeNodeDataObject* self = (TreeNodeDataObject *)type->tp_alloc(type, 0);

  if( self != NULL ) {
    self->taxon = NULL;
    self->branchlength = NULL;
    self->height = NULL;
    self->allData = PyDict_New();
  }

  return self;
}

static int
TreeNodeData_init(TreeNodeDataObject*  self,
		  const char*          taxon,
		  const double*        branchlength,
		  const double*        height,
		  const Attributes*    attributes)
{
  if( taxon ) {
    self->taxon = PyString_FromString(taxon);
  } else {
    Py_INCREF(Py_None);
    self->taxon = Py_None;
  }
  if( branchlength ) {
    self->branchlength = PyFloat_FromDouble(*branchlength);
  } else {
    Py_INCREF(Py_None);
    self->branchlength = Py_None;
  }
  if( height ) {
    self->height = PyFloat_FromDouble(*height);
  } else {
    Py_INCREF(Py_None);
    self->height = Py_None;
  }
  if( attributes ) {
    int o = PyDict_SetItemString(self->allData, "attributes", attributesAsPyObj(attributes));
    if( o == -1 ) {
      return -1;
    }
  }
  return 0;
}

static void
TreeNodeData_dealloc(TreeNodeDataObject* self)
{
  Py_XINCREF(self->taxon);
  Py_XINCREF(self->branchlength);
  Py_XINCREF(self->height);
  Py_XINCREF(self->allData);
  self->ob_type->tp_free((PyObject*)self);
}

PyObject*
TreeNodeData_getattr(TreeNodeDataObject* self, PyObject* aname)
{
  PyObject* item = PyDict_GetItem(self->allData, aname);
  if( item ) {
    Py_INCREF(item);
    return item;
  }
  return PyObject_GenericGetAttr(self, aname);
}

int
TreeNodeData_setattr(TreeObject* self, PyObject* name, PyObject* value)
{
  const char* const cname = PyString_AS_STRING(name);
  if( ! strcmp(cname, "branchlength") || ! strcmp(cname, "height") ) {
    PyErr_SetString(PyExc_RuntimeError, "Please set node branchlength/height via tree method.") ;
    return -1;
  }
  return PyObject_GenericSetAttr(self, name, value);
}

static PyTypeObject TreeNodeDataType = {
  PyObject_HEAD_INIT(NULL)
  0,				/* ob_size        */
  "treesset.NodeData",		/* tp_name        */
  sizeof(TreeNodeDataObject),	/* tp_basicsize   */
  0,				/* tp_itemsize    */
  (destructor)TreeNodeData_dealloc,		/* tp_dealloc     */
  0,				/* tp_print       */
  0,				/* tp_getattr     */
  0,				/* tp_setattr     */
  0,				/* tp_compare     */
  0,				/* tp_repr        */
  0,				/* tp_as_number   */
  0,				/* tp_as_sequence */
  0,				/* tp_as_mapping  */
  0,				/* tp_hash        */
  0,				/* tp_call        */
  0,				/* tp_str         */
  0/*TreeNodeData_getattr*/,		/* tp_getattro    */
  (setattrofunc)TreeNodeData_setattr,	/* tp_setattro    */
  0,				/* tp_as_buffer   */
  Py_TPFLAGS_DEFAULT,		/* tp_flags       */
  "tree Node data.",		/* tp_doc         */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  0,             		/* tp_methods */
  NodeData_members,             /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  offsetof(TreeNodeDataObject,allData),  /* tp_dictoffset */
  0,                         /* tp_init */
  0,                         /* tp_alloc */
  (newfunc)TreeNodeData_new,         /* tp_new */
};


struct TreeNodeObject : PyObject {
  PyObject* prev;
  PyObject* succ;
  PyObject* data;
};

static PyMemberDef Node_members[] = {
  {(char*)"succ", T_OBJECT_EX, offsetof(TreeNodeObject, succ), 0,
    (char*)"succ"},
  {(char*)"prev", T_OBJECT_EX, offsetof(TreeNodeObject, prev), 0,
    (char*)"prev"},
  {(char*)"data", T_OBJECT_EX, offsetof(TreeNodeObject, data), 0,
    (char*)"data"},
  {NULL}  /* Sentinel */
};

static TreeNodeObject *
TreeNode_new(PyTypeObject* type, PyObject *args, PyObject *kwds)
{
  TreeNodeObject* self = static_cast<TreeNodeObject*>(type->tp_alloc(type, 0));

  if( self != NULL ) {
    self->succ = NULL;
    self->prev = NULL;
  }

  return self;
}

static void
TreeNode_init(TreeNodeObject* self, int prev, uint nSons, uint* sons, PyObject* data)
{
  if( prev >= 0 ) {
    self->prev = PyInt_FromLong(prev);
  } else {
    Py_INCREF(Py_None);
    self->prev = Py_None;
  }
  if( nSons == 0 ) {
    Py_INCREF(Py_None);
    self->succ = Py_None;
  } else {
    self->succ = PyTuple_New(nSons);
    for(uint k = 0; k < nSons; ++k) {
      PyTuple_SET_ITEM(self->succ , k, PyInt_FromLong(sons[k]));
    }
  }

  self->data = data;
}

static void
TreeNode_dealloc(TreeNodeObject* self)
{
  Py_XINCREF(self->succ);
  Py_XINCREF(self->prev);
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject TreeNodeType = {
  PyObject_HEAD_INIT(NULL)
  0,				/* ob_size        */
  "treesset.Node",		/* tp_name        */
  sizeof(TreeNodeObject),	/* tp_basicsize   */
  0,				/* tp_itemsize    */
  (destructor)TreeNode_dealloc,	/* tp_dealloc     */
  0,				/* tp_print       */
  0,				/* tp_getattr     */
  0,				/* tp_setattr     */
  0,				/* tp_compare     */
  0,				/* tp_repr        */
  0,				/* tp_as_number   */
  0,				/* tp_as_sequence */
  0,				/* tp_as_mapping  */
  0,				/* tp_hash        */
  0,				/* tp_call        */
  0,				/* tp_str         */
  0,				/* tp_getattro    */
  0,				/* tp_setattro    */
  0,				/* tp_as_buffer   */
  Py_TPFLAGS_DEFAULT,		/* tp_flags       */
  "Tree Node.",	                /* tp_doc         */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  0,             		/* tp_methods */
  Node_members,                 /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  0/*(initproc)TreeNode_init*/,      /* tp_init */
  0,                         /* tp_alloc */
  (newfunc)TreeNode_new,         /* tp_new */
};


class TreesSet {
public:
  TreesSet(bool isCompressed, uint _precision, bool s);
  ~TreesSet();

  // Add a tree from text in NEWICK format.
  int add(const char* txt, PyObject* kwds);
  
  uint nTrees(void) const { return trees.size(); }
  
  TreeRep const& getTree(uint const i) const {
                                 assert( i < nTrees() );
    return *trees[i];
  }
  
  // taxon name (existing,string) from internal index.
  string const& taxonString(uint const k) const {
                                 assert(k < taxaList.size());
    return taxaList[k];
  }

  // Index of taxon if exists, -1 otherwise
  int hasTaxon(const char* taxon) const;
  
  // Populate hs/txhs with internal node/taxa heights for nt'th tree.
  void getHeights(uint nt, vector<double>& hs, vector<double>& txhs) const;

  void setTreeAttributes(uint nt, TreeObject* to) const;
  
  // Compressed/non-Compressed tree
  bool const compressed : 8;
  bool const store      : 8;
  // floating point precision (float or double)
  uint const precision  : 8;

  vector< vector<ParsedTreeNode> > asNodes;

  void add(TreesSet const& ts, uint const nt, vector<uint> const& filteredTaxa);

private:
  TreeRep*  repFromData(bool const                  cladogram,
			vector<uint> const&         taxa,
			uint const                  maxTaxaIndex,
			vector<double> const&       heights,
			vector<double>* const       taxaHeights,
			vector<uint>* const         labels,
			vector<const Attributes*>*  atrs);
  
  // Encodes a parsed tree 
  TreeRep*	nodes2rep(vector<ParsedTreeNode>& nodes);

  // taxon index (inserts new ones). 
  uint 		getTaxon(string const& taxon);

  // All trees
  vector<TreeRep*>		trees;

  vector<PyObject*>		treesAttributes;
  
  // taxa across all trees
  vector<string>		taxaList;
  // mapping from taxon to its position in taxaList 
  unordered_map<string,uint>	taxaDict;
};

TreesSet::TreesSet(bool isCompressed, uint _precision, bool s) :
  compressed(isCompressed),
  store(s),
  precision(_precision)
{}

TreesSet::~TreesSet()
{
  for(auto t = trees.begin(); t != trees.end(); ++t) {
    delete *t;
  }
  
  for(auto a = treesAttributes.begin(); a != treesAttributes.end(); ++a) {
    Py_XDECREF(*a);
  }
}

int
TreesSet::hasTaxon(const char* taxon) const {
  auto const i = taxaDict.find(taxon);
  if( i == taxaDict.end() ) {
    return -1;
  }
  return i->second;
}

uint
TreesSet::getTaxon(string const& taxon)
{
  auto const i = taxaDict.find(taxon);
  if( i == taxaDict.end() ) {
    int const k = taxaList.size();
    taxaList.push_back(taxon);
    taxaDict.insert( std::pair<string,uint>(taxon, k) );
    return k;
  }
  return i->second;
}

void
TreesSet::getHeights(uint nt, vector<double>& hs, vector<double>& txhs) const
{
  TreeRep const& r = getTree(nt);
  if( r.isCladogram() ) {
    CladogramRep const& c = static_cast<CladogramRep const&>(r);
    vector<uint> const& h = c.heights();
    hs.assign(h.begin(), h.end());
  } else {
    if( precision == 8 ) {
      hs = static_cast<PhylogramRep<double> const&>(r).heights();
      auto tx = static_cast<PhylogramRep<double> const&>(r).txheights();
      if( tx ) {
	txhs = *tx;
      }
    } else {
      auto const& h = static_cast<PhylogramRep<float> const&>(r).heights();
      hs.assign(h.begin(), h.end());
      auto tx = static_cast<PhylogramRep<float> const&>(r).txheights();
      if( tx ) {
	txhs.assign(tx->begin(), tx->end());
      }
    }
  }
}

void
TreesSet::setTreeAttributes(uint nt, TreeObject* to) const
{
  PyObject* a = treesAttributes[nt];
  if( a && PyDict_Size(a) > 0 ) {
    PyDict_Update(to->dict, a);
  }
}

TreeRep*
TreesSet::repFromData(bool const                  cladogram,
		      vector<uint> const&         taxa,
		      uint const                  maxTaxaIndex,
		      vector<double> const&       heights,
		      vector<double>* const       taxaHeights,
		      vector<uint>* const         labels,
		      vector<const Attributes*>*  atrs)
{
  Packer<uint>* top = 0;
  
  if( compressed ) {
    int nbitsStoreTaxa = lg2i(maxTaxaIndex) + 1;
    top = new FixedIntPacker(nbitsStoreTaxa, taxa.begin(), taxa.end());
  } else {
    top = new SimplePacker<uint>(taxa);
  }

  FixedIntPacker* lb = 0;
  if( labels ) {
    int const nbitsStoreH = lg2i(*std::max_element(labels->begin(),labels->end())) + 1;
    lb = new FixedIntPacker(nbitsStoreH, labels->begin(), labels->end());
  }

  TreeRep* r;
  
  if( cladogram ) {
    vector<uint> hs(heights.size());
    for(uint i = 0; i < hs.size(); ++i) {
      hs[i] = static_cast<uint>(heights[i] + 0.5);
    }
    Packer<uint>* hsb;
    // single leaf trees has no heights
    if( compressed && heights.size() > 0 ) {
      int const nbitsStoreH = lg2i(*std::max_element(hs.begin(),hs.end())) + 1;
      hsb = new FixedIntPacker(nbitsStoreH, hs.begin(), hs.end());
    } else {
      hsb = new SimplePacker<uint>(hs);
    }
    r = new CladogramRep(*top, lb, hsb, atrs);
  } else {
    if( precision == 8 ) {
      SimplePacker<double>* hsb = new SimplePacker<double>(heights);
      SimplePacker<double>* txhs = 0;
      if( taxaHeights ) {
	txhs = new SimplePacker<double>(*taxaHeights);
      }
      r = new PhylogramRep<double>(*top, lb, hsb, txhs, atrs);
    } else {
      SimplePacker<float>* hsb = new SimplePacker<float>(heights);
      SimplePacker<float>* txhs = 0;
      if( taxaHeights ) {
	txhs = new SimplePacker<float>(*taxaHeights);
      }
      r = new PhylogramRep<float>(*top, lb, hsb, txhs, atrs);
    }
  }
  return r;
}

TreeRep*
TreesSet::nodes2rep(vector<ParsedTreeNode>& nodes)
{
  // tree taxa (as indices into global table)
  vector<uint> taxa;
  bool cladogram = true;
  bool hasAttributes = false;
  bool hasInternalLabels = false;
  uint maxTaxaIndex = 0;
  
  for(auto n = nodes.begin(); n != nodes.end() ; ++n) {
    if( n->sons.size() == 0 ) {
      // can have un-named taxon
      uint const k = getTaxon(n->taxon);
      maxTaxaIndex = std::max(maxTaxaIndex, k);
      taxa.push_back(k);
    } else if( n->taxon.size() ) {
      hasInternalLabels = true;
    }
    if( n->branch ) {
      cladogram = false;
    }
    if( n->attributes ) {
      hasAttributes = true;
    }
  }

  uint const nTaxa = taxa.size();
  
  // heights of internal nodes (between taxa)
  vector<double> heights(nTaxa-1, -1);
  
  // if nodes[x] is a tip, locs[x]+1 is the ordinal number of the tip in the
  // array of tips. For internal nodes, locs[x]+1 is the ordinal number of the
  // node in the heights array, where x is the first son of that internal node.
  vector<int> locs(nodes.size());
  uint iloc = 0;
  
  // taxa heights (when tips are not contemporaneous) 
  vector<double>* taxaHeights = 0;
  
  for(auto n = nodes.begin(); n != nodes.end() ; ++n) {
    if( n->sons.size() == 0 ) {
      if( ! n->branch ) {
	n->branch = new double(1);
      }
                                               assert( 0 <= iloc && iloc < locs.size() );
      locs[iloc] = iloc == 0 ? -1 : locs[iloc-1]+1;
      ++iloc;
    } else {
      // node own height
      double h = -1;
      for(auto s = n->sons.begin(); s != n->sons.end(); ++s) {
	// After processing the node, we use the branch to store the height.
	h = std::max(*nodes[*s].branch, h);
      }
      if( ! cladogram ) {
	for(auto s = n->sons.begin(); s != n->sons.end(); ++s) {
	  double const h1 = *nodes[*s].branch;
	  double const dh = h - h1;
	  if( dh > 0 && !areSame(h,h1) ) {
	    if( ! taxaHeights ) {
	      taxaHeights = new vector<double>(nTaxa, 0.0);
	    }
	    // Add 'dh' to the clade: increase all taxa and internal node heights
	    std::stack<uint> hp;
	    hp.push(*s);
	    // a clade is consecutive
	    int lo = 4*nTaxa,hi = -1;
	    while( ! hp.empty() ) {
	      uint const x = hp.top(); hp.pop();
	      if( nodes[x].sons.size() == 0 ) {
		int const l = locs[x]+1;
		lo = std::min(l, lo);
		hi = std::max(l, hi);
		(*taxaHeights)[l] += dh;
	      } else {
		auto const& b = nodes[x].sons;
		for(auto a = b.begin(); a != b.end(); ++a) {
		  hp.push(*a);
		}
	      }
	    }
	    while( lo < hi ) {
	      heights[lo] += dh;
	      ++lo;
	    }
	  }
	}
      }
      // node height stored between all sons
      for(auto s = n->sons.begin(); s != n->sons.end()-1; ++s) {
	uint const i = locs[*s]+1;            assert( 0 <= i && i < heights.size() && heights[i] == -1);
	heights[i] = h;
      }
                                              assert( 0 < iloc && iloc < locs.size() );
      locs[iloc] = locs[iloc-1]; ++iloc;
      // Use the branch to store the node height
      if( n->branch ) {
	*n->branch += h;
      } else if( n+1 != nodes.end() ) {
	n->branch = new double(h+1);
      }
    }
  }

  vector<const Attributes*>* atrs = 0;
  vector<uint>* labels = 0;
  if( hasAttributes || hasInternalLabels ) {
    atrs = hasAttributes ? new vector<const Attributes*>(nodes.size(),0) : 0;
    labels = hasInternalLabels ? new vector<uint>(nTaxa-1, 0) : 0;

    for(auto n = nodes.begin(); n != nodes.end() ; ++n) {
      bool const isTip = n->sons.size() == 0;
    
      if( n->attributes ) {
	int const l = isTip ? locs[n - nodes.begin()]+1 : locs[n->sons[0]]+1+nTaxa;
	assert( atrs->at(l) == 0 );
	// steal
	(*atrs)[l] = n->attributes;
	n->attributes = 0;
      }
      
      if( !isTip && n->taxon.size() ) {
	uint const k = getTaxon(n->taxon);
	int const l = locs[n->sons[0]]+1;
	assert( labels->at(l) == 0 );
	(*labels)[l] = k+1;
      }
    }
  }

  TreeRep* const r = repFromData(cladogram, taxa, maxTaxaIndex, heights, taxaHeights, labels, atrs);
  
  if( taxaHeights ) {
    delete taxaHeights;
  }
  if( labels ) {
    delete labels;
  }

  return r;  
}

int
TreesSet::add(const char* treeTxt, PyObject* kwds)
{
  vector<ParsedTreeNode> nodes;

  int const txtLen = strlen(treeTxt);
  int nc = readSubTree(treeTxt, nodes);

  if( nc > 0 ) {
    nc += skipSpaces(treeTxt + nc);
  }
  
  if( ! (nc == txtLen || (nc+1 == txtLen && treeTxt[nc] == ';')) ) {
    if( nc < 0) {
      int const where = txtLen+(nc+1);
      PyErr_Format(PyExc_ValueError, "failed parsing around %d (%10.10s ...).", where, treeTxt+where);
    } else {
      PyErr_Format(PyExc_ValueError, "extraneous characters at tree end: '%s'",
		   string(treeTxt+nc,std::max(5,txtLen-nc)).c_str());
    }
    return -1;
  }

  if( kwds ) {
    Py_INCREF(kwds);
  }
  treesAttributes.push_back(kwds);
  
  if( store ) {
    asNodes.push_back(nodes);
    return asNodes.size()-1;
  } else {
    trees.push_back(nodes2rep(nodes));
    return trees.size()-1;
  }
}

void
TreesSet::add(TreesSet const& ts, uint const nt, vector<uint> const& filteredTaxa)
{
  vector<double> hs;
  vector<double> txhs;

  TreeRep const& rep = ts.getTree(nt);
  bool isp;
  vector<uint> const& t = rep.tips(isp);
  vector<uint> const& tax = isp ? t : vector<uint>(t);
  
  uint const nTaxa = tax.size();
  ts.getHeights(nt, hs, txhs);
  bool const hasTXheights = txhs.size() > 0;

  auto atrb = rep.getAttributes();
  vector<uint>* labels = rep.labels();

  vector<double> newhs;
  vector<double> newtxhs;
  vector<uint>   newTips;
  vector<uint>* newLabels = labels ? new vector<uint> : 0;

  // populate newTips with indices (into 'tax') of kept taxa 
  for(uint k = 0; k < nTaxa; ++k) {
    if( std::find(filteredTaxa.begin(), filteredTaxa.end(), tax[k]) == filteredTaxa.end() ) {
      newTips.push_back(k);
    }
  }

  vector<const Attributes*>* newAtrbs =
    atrb ? new vector<const Attributes*>(2*newTips.size()-1,0) : 0;

  uint const newNtaxa = newTips.size();
  uint s;

  for(uint k = 0; k < newNtaxa; ++k) {
    uint const e = newTips[k];

    newTips[k] = getTaxon(ts.taxonString(tax[e]));

    if( hasTXheights ) {
      newtxhs.push_back(txhs[e]);
    }

    if( newAtrbs ) {
      auto a = (*atrb)[e];
      (*newAtrbs)[k] = a ? new Attributes(*a) : 0;
    }
    
    if( k > 0 ) {
      uint const i = std::max_element(hs.begin()+s, hs.begin()+e) - hs.begin();
    
      newhs.push_back(hs[i]);

      if( newAtrbs ) {
	auto a = (*atrb)[i+nTaxa];
	(*newAtrbs)[k-1 + newNtaxa] = a ? new Attributes(*a) : 0;
      }
      
      if( newLabels ) {
	uint l = (*labels)[i];
	if( l > 0 ) {
	  l = getTaxon(ts.taxonString(l-1))+1;
	}
	newLabels->push_back(l);      
      }
    }
    s = e;
  }
  
  uint const maxTaxaIndex = *std::max_element(newTips.begin(), newTips.end());

  TreeRep* r = repFromData(rep.isCladogram(), newTips, maxTaxaIndex, newhs,
			   hasTXheights ? &newtxhs : 0, newLabels, newAtrbs);
  trees.push_back(r);
  PyObject* a = ts.treesAttributes[nt];
  Py_XINCREF(a);
  treesAttributes.push_back(a);
}

Tree::Expanded::Expanded(int                _itax,
			 uint               _nSons,
			 uint*              _sons,
			 double*            _branch,
			 double             _height,
			 int                _prev,
			 const Attributes*  _attributes) :
  itax(_itax),
  nSons(_nSons),
  sons(_sons),
  branch(_branch),
  height(_height),
  prev(_prev),
  attributes(_attributes)
{}
  

uint Tree::rep2treeInternal(vector<Expanded>&                       nodes,
			    uint                                    low,
			    uint const                              hi,
			    vector<uint> const&                     tax,
			    vector<double> const&                   htax,
			    vector<double> const&                   hs,
			    const vector<const Attributes*>* const  atrbs,
			    const vector<uint>* const               labels,
			    uint*&                                  sonsBlock,
			    uint* const                             curiScratch,
			    uint                                    bleft) const
{
  if( low == hi ) {
    nodes.push_back(Expanded(tax[low],0,0,0,htax[low],-1, atrbs ? (*atrbs)[low] : 0));
  } else {
    uint* curi = 0;
    double curh = -1;
    for(uint k = low; k < hi; ++k) {
      double const h = hs[k];
      if( h >= curh ) {
	if( h > curh ) {
	  curh = h;
	  *curiScratch = k;
	  curi = curiScratch+1;
	} else {
	  *curi = k; ++curi;
	}
      }
    }
    uint const nSplits = (uint)(curi-curiScratch);        assert( nSplits <= bleft );
    bleft -= nSplits;
    
    uint* sons = sonsBlock;
    uint nSons = 1 + nSplits;
    sonsBlock += nSons;
    ;                                                      assert( bleft > 0 );
    *curi = hi; ++curi; --bleft;
    
    for(uint* x = curiScratch; x < curi; ++x) {
      uint const k = rep2treeInternal(nodes, low, *x, tax, htax,
				      hs, atrbs, labels, sonsBlock, curi, bleft);
      double const hs = nodes[k].height;
      nodes[k].branch = new double(curh - hs);
      *sons = k; ++sons;
      low = *x+1;
    }

    sons -= nSons;
    for(uint* s = sons; s < sons+nSons; ++s) {
      nodes[*s].prev = nodes.size();
    }
    int iTax = -1;
    if( labels ) {
      uint nl = (*labels)[*curiScratch];
      if( nl > 0 ) {
	iTax = nl - 1;
      }
    }
      
    nodes.push_back(Expanded(iTax, nSons, sons, 0, curh, -1,
			     atrbs ? (*atrbs)[*curiScratch + tax.size()] : 0));
  }
  return nodes.size()-1;
}
    

void Tree::setup(void) const {
  if( internals ) {
    return;
  }
  
  vector<double> hs;
  vector<double> txhs;

  TreeRep const& rep = ts.getTree(nt);
  bool isp;
  vector<uint> const& t = rep.tips(isp);
  vector<uint> const& tax = isp ? t : vector<uint>(t);
  vector<uint> const* labels = rep.labels();
  
  nTaxa = tax.size();
  ts.getHeights(nt, hs, txhs);
  if( txhs.size() == 0 ) {
    txhs.resize(nTaxa, 0.0);
  }

  internals = new vector<Expanded>;
  sonsBlockSave = new uint[2*nTaxa];  // allocate in one contiguous block 
  uint block[2*nTaxa];     // scratch only
  uint* s = sonsBlockSave; // keep sonsBlockSave safe
  rep2treeInternal(*internals, 0, hs.size(), tax, txhs, hs,
		   rep.getAttributes(), labels, s, block, 2*nTaxa);
  if( labels ) {
    delete labels;
  }
  if( rep.isCladogram() ) {
    for(auto i = internals->begin(); i != internals->end(); ++i) {
      Expanded& x = *i;
      if( x.branch ) {
	delete x.branch;
	x.branch = 0;
      }
      x.height = -1;
    }
  }
}

	
vector<uint> const&
Tree::tips() const
{
  return ts.getTree(nt).tips();
}

void
Tree::getTerminals(vector<uint>& terms) const
{
  setup();
  terms.reserve(nTaxa);
  for(auto i = internals->begin(); i != internals->end(); ++i) {
    Expanded& x = *i;
    if( x.itax >= 0 ) {
      terms.push_back(x.itax);
    }
  }
}

uint
Tree::nNodes(void) const
{
  setup();
  return internals->size();
}

bool
Tree::isCladogram(void) const
{
  return ts.getTree(nt).isCladogram();
}

void
Tree::tostr(vector<string>&  s,
	    uint const       nodeId,
	    bool const       topoOnly,
	    bool const       includeStem,
	    bool const       withAttributes) const
{
  Expanded const& n = (*internals)[nodeId];
  if( n.nSons == 0 ) {
    if( n.itax >= 0 ) {
      s.push_back(ts.taxonString(n.itax));
    } else {
      s.push_back("");
    }
  } else {
    for(uint i = 0; i < n.nSons; ++i) {
      tostr(s, n.sons[i], topoOnly, true, withAttributes);
    }
    auto first = s.end() - n.nSons;
    std::sort(first, s.end());
    first->insert(first->begin(),'(');
    for(uint i = 1; i < n.nSons; ++i) {
      first->append(",").append(*(first + i));
    }
    first->append(")");
    s.erase(first+1, s.end());

    if( n.itax >= 0 ) {
      first->append(ts.taxonString(n.itax));
    }
  }

  if( withAttributes && n.attributes ) {
    Attributes const& atr = *n.attributes;
    string& a = *(s.end()-1);
    a.append("[&");
    for(uint k = 0; k < atr.size(); ++k) {
      auto const p = atr[k];
      if( k > 0 ) {
	a.append(",");
      }
      a.append(p.first).append("=").append(p.second);
    }
    a.append("]");
  }
  
  if( ! topoOnly && n.branch && includeStem ) {
    char* const b = PyOS_double_to_string(*n.branch, 'r', 0, Py_DTSF_ADD_DOT_0, 0);
    (s.end()-1)->append(":").append(b);
    PyMem_Free(b);
  }
}

void
Tree::toNewick(string& s, int nodeId, bool const topoOnly,
	       bool const includeStem, bool const withAttributes) const
{
  if( nodeId == -1 ) {
    nodeId = getRootID();
  } else {
    setup();
  }
  vector<string> c;
  tostr(c, static_cast<uint>(nodeId), topoOnly, includeStem, withAttributes);
  s.assign(*c.begin());
}

struct TreesSetObject : PyObject {
  TreesSet* ts;
  
  void del(void);
  void init(void);
  PyObject* taxon(uint k);
  PyObject* newNodeData(const char* tx, const double* branch, const double *height);

private:
  vector<PyObject*>* taxa;
};


void
TreesSetObject::del(void)
{
  delete ts;
  for(auto s = taxa->begin(); s != taxa->end(); ++s) {
    Py_XDECREF(*s);
  }
  delete taxa;
}

void
TreesSetObject::init(void)
{
  ts = NULL;
  taxa = new vector<PyObject*>();
}

PyObject*
TreesSetObject::taxon(uint k)
{
  if( taxa->size() <= k ) {
    taxa->resize(k+1, 0);
  }
  if( (*taxa)[k] == 0 ) {
    string const& s = ts->taxonString(k);
    (*taxa)[k] = PyString_FromString(s.c_str());
  }
    
  PyObject* s = (*taxa)[k];
  Py_INCREF(s);
  return s;
}

static void
TreesSet_dealloc(TreesSetObject* self)
{
  self->del();
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject*
TreesSet_new(PyTypeObject* type, PyObject* /*args*/, PyObject* /*kwds*/)
{
  TreesSetObject* self = (TreesSetObject *)type->tp_alloc(type, 0);

  if (self != NULL) {
    self->init();
  }

  return self;
}

static int
TreesSet_init(TreesSetObject* self, PyObject* args, PyObject* kwds)
{
  static const char *kwlist[] = {"compressed", "precision", "store",
				 static_cast<const char*>(0)};
  PyObject* comp = 0;
  PyObject* sto = 0;
  int precision = 4;
  
  if (! PyArg_ParseTupleAndKeywords(args, kwds, "|OiO", (char**)kwlist,
				    &comp,&precision,&sto)) {
    return -1;
  }

  if( ! (precision == 4 || precision == 8) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args (precision)");
    return -1;
  }
    
  bool const compressed = (! comp || PyObject_IsTrue(comp));
  bool const store = (sto && PyObject_IsTrue(sto));
  
  self->ts = new TreesSet(compressed, precision, store);
  return 0;
}

static Py_ssize_t
TreesSet_len(TreesSetObject* self) 
{
  if( self->ts->store ) {
    return self->ts->asNodes.size();
  }
  return self->ts->nTrees();
}


PyObject*
TreesSet_getItemInternals(TreesSetObject* self, uint i)
{
  auto const& ts = *self->ts;
  if( ts.store ) {
    auto const& ns = ts.asNodes[i];
    PyObject* n = PyTuple_New(ns.size());
    for(uint k = 0; k < ns.size(); ++k) {
      PyTuple_SET_ITEM(n, k, ns[k].asPyObject());
    }
    return n;
  }
    
  if( static_cast<uint>(i) >= ts.nTrees() ) {
    PyErr_SetNone(PyExc_IndexError);
    return 0;
  }
  TreeRep const& r = ts.getTree(i);
  PyObject* n = PyTuple_New(5);
  bool const isc = r.isCladogram();
  PyTuple_SET_ITEM(n, 0, PyBool_FromLong(isc));
  vector<uint> const& topo = r.tips();
  PyObject* t = PyTuple_New(topo.size());
  for(uint k = 0; k < topo.size(); ++k) {
    PyTuple_SET_ITEM(t, k, self->taxon(topo[k]));
  }
  PyTuple_SET_ITEM(n, 1, t);
  if( isc ) {
    CladogramRep const& c = static_cast<CladogramRep const&>(r);
    vector<uint> const& hs = c.heights();
    PyObject* h = PyTuple_New(hs.size());
    for(uint k = 0; k < hs.size(); ++k) {
      PyTuple_SET_ITEM(h, k, PyInt_FromLong(hs[k]));
    }
    PyTuple_SET_ITEM(n, 2, h);

    Py_INCREF(Py_None);
    PyTuple_SET_ITEM(n, 3, Py_None);

    Py_INCREF(Py_None);
    PyTuple_SET_ITEM(n, 3, Py_None);
  } else {
    if( self->ts->precision == 4 ) {
      typedef PhylogramRep<float> T;
      
      T const& p = static_cast<T const&>(r);
      PyTuple_SET_ITEM(n, 2, dvector2tuple(p.heights()));
      PyTuple_SET_ITEM(n, 3, dvector2tuple(*p.txheights()));
    } else {
      typedef PhylogramRep<double> T;
      
      T const& p = static_cast<T const&>(r);
      PyTuple_SET_ITEM(n, 2, dvector2tuple(p.heights()));
      PyTuple_SET_ITEM(n, 3, dvector2tuple(*p.txheights()));
    }
  }
  auto a = r.getAttributes();
  if( a ) {
    PyObject* ap = PyTuple_New(a->size());
    for(uint k = 0; k < a->size(); ++k) {
      PyTuple_SET_ITEM(ap, k, attributesAsPyObj((*a)[k]));
    }
    PyTuple_SET_ITEM(n, 4,ap);
  } else {
    Py_INCREF(Py_None);
    PyTuple_SET_ITEM(n, 4, Py_None);
  }
  return n;
}

PyObject*
TreesSet_getItem(TreesSetObject* self, Py_ssize_t i)
{
  auto const& ts = *self->ts;
  uint const nt = static_cast<uint>(i);
  if( nt >= ts.nTrees() ) {
    PyErr_SetNone(PyExc_IndexError);
    return 0;
  }
  TreeObject* to = Tree_new(&TreeType,0,0);
  to->init(self, new Tree(ts, nt));
  ts.setTreeAttributes(nt, to);
  return to;
}


PyObject*
TreeObject::getTaxa(void) const
{
  if( taxa == 0 ) {
    auto const& topo = tr->tips();
    taxa = PyTuple_New(topo.size());
    for(uint k = 0; k < topo.size(); ++k) {
      PyTuple_SET_ITEM(taxa, k, ts->taxon(topo[k]));
    }
    // for efficiency: keep one refcount on creation, so this tuple stays
    // around until tree is deleted.
  }
  Py_INCREF(taxa);
  return taxa;
}

PyObject*
TreeObject::getTerminals(void) const
{
  vector<uint> terms;

  tr->getTerminals(terms);

  PyObject* t = PyTuple_New(terms.size());
  for(uint k = 0; k < terms.size(); ++k) {
    PyTuple_SET_ITEM(t, k, PyInt_FromLong(terms[k]));
  }
  return t;
}

PyObject*
TreeObject::allIds(void) const 
{
  uint const n = tr->nNodes();
  PyObject* t = PyTuple_New(n);
  for(uint k = 0; k < n; ++k) {
    PyTuple_SET_ITEM(t, k, PyInt_FromLong(k));
  }
  return t;
}

PyObject*
TreeObject::getNode(uint nt) const
{
  if( nt >= tr->nNodes() ) {
    PyErr_SetNone(PyExc_IndexError);
    return 0;
  }
  
  if( treeNodes ) {
    TreeNodeObject* n = (*treeNodes)[nt];
    if( n ) {
      Py_INCREF(n);
      return n;
    }
  } else {
    treeNodes = new vector<TreeNodeObject*>(tr->nNodes(),0);
  } 
  
  Tree::Expanded const& e = tr->getNode(nt);

  const char* tx = 0;
  if( e.itax >= 0 ) {
    tx = tr->ts.taxonString(e.itax).c_str();
  }
  bool const isc = tr->isCladogram();
  TreeNodeObject* node = TreeNode_new(&TreeNodeType, 0, 0);

  TreeNodeDataObject* d = TreeNodeData_new(&TreeNodeDataType, 0, 0);
  int succ = TreeNodeData_init(d, tx, isc ? 0 : e.branch, isc ? 0 : &e.height, e.attributes);
  if( succ == -1 ) {
    // leaks, not expected to happen
    PyErr_SetString(PyExc_RuntimeError, "memory error.") ;
    return 0;
  }
  
  TreeNode_init(node, e.prev, e.nSons, e.sons, d);

  // keep it around
  Py_INCREF(node);
  (*treeNodes)[nt] = node;
  
  return node;
}


void
TreeObject::getInOrder(bool preOrder, vector<int>& ids,
		       int nodeId, bool includeTaxa)
{
  uint nSons;
  if( treeNodes && (*treeNodes)[nodeId] ) {
    TreeNodeObject const& no = *(*treeNodes)[nodeId];
    nSons = no.succ == Py_None ? 0 : PySequence_Size(no.succ);
    if( nSons > 0 && preOrder ) {
      ids.push_back(nodeId);
    }
    for(uint i = 0; i < nSons; ++i) {
      int d = PyLong_AsLong(PySequence_GetItem(no.succ, i));
      getInOrder(preOrder,ids, d,includeTaxa);
    }
  } else {
    Tree::Expanded const& no = tr->getNode(nodeId);
    nSons = no.nSons;
    if( nSons > 0 && preOrder ) {
      ids.push_back(nodeId);
    }
    for(uint i = 0; i < nSons; ++i) {
      getInOrder(preOrder, ids, no.sons[i], includeTaxa);
    }
  }
  if( (nSons > 0 && !preOrder) || (nSons == 0 && includeTaxa) ) {
    ids.push_back(nodeId);
  }
}

void
TreeObject::setBranch(uint const nodeId, double const newBranchLen)
{
  TreeNodeObject& no = *static_cast<TreeNodeObject*>(getNode(nodeId));

  TreeNodeDataObject& data = *static_cast<TreeNodeDataObject*>(no.data);

  double const brlen = data.setBranch(newBranchLen);
  
  double const dif = newBranchLen - brlen;
  if( dif == 0 ) {
    return; // None??
  }
  double minNewHeight = PyFloat_GetMax();
  vector<int> subt;
  getInOrder(false, subt, nodeId, true);
  for(auto n = subt.begin(); n != subt.end(); ++n) {
    TreeNodeObject& x = *static_cast<TreeNodeObject*>(getNode(*n));
    TreeNodeDataObject& d = *static_cast<TreeNodeDataObject*>(x.data);
    if( d.hasHeight() ) {
      double const h = d.adjustHeight(-dif);
      minNewHeight = std::min(h, minNewHeight);
    }
  }
  if( minNewHeight < 0 ) {
    uint const n = tr->nNodes();
    for(uint k = 0; k < n; ++k) {
      TreeNodeObject& x = *static_cast<TreeNodeObject*>(getNode(k));
      TreeNodeDataObject& d = *static_cast<TreeNodeDataObject*>(x.data);
      if( d.hasHeight() ) {
	d.adjustHeight(-minNewHeight);
      }
    }
  }
}

string
strOfObject(PyObject* o)
{
  if( PyString_Check(o) ) {
    return PyString_AsString(o);
  }
  PyObject* s = PyObject_Str(o);
  string tmp(PyString_AsString(s));
  Py_DECREF(s);
  return tmp;
}
  
void
TreeObject::tostr(vector<string>& s, int const nodeId,
		  bool const topoOnly, bool const includeStem,
		  bool const withAttributes) const
{
  TreeNodeObject* np = (*treeNodes)[nodeId];
  if( ! np ) {
    // shall we delete this before returning??
    np = static_cast<TreeNodeObject*>(getNode(nodeId));
  }
  TreeNodeObject const& no = *np;
  uint const nSons = no.succ == Py_None ? 0 : PySequence_Size(no.succ);
  TreeNodeDataObject const& nd = *static_cast<TreeNodeDataObject*>(no.data);
  if( nSons == 0 ) {
    char* t = PyString_AsString(nd.taxon);
    s.push_back(t);
  } else {
    for(uint i = 0; i < nSons; ++i) {
      int d = PyLong_AsLong(PySequence_GetItem(no.succ, i));
      tostr(s, d, topoOnly, true, withAttributes);
    }
    auto first = s.end() - nSons;
    std::sort(first, s.end());
    first->insert(first->begin(),'(');
    for(uint i = 1; i < nSons; ++i) {
      first->append(",").append(*(first + i));
    }
    first->append(")");
    s.erase(first+1, s.end());

    if( nd.taxon != Py_None ) {
      char* t = PyString_AsString(nd.taxon);
      first->append(t);
    }
  }

  if( withAttributes && PyObject_HasAttrString(no.data, "attributes") ) {
    PyObject* const atrs = PyObject_GetAttrString(no.data, "attributes");
    if( atrs && PyMapping_Check(atrs) ) {
      string& a = *(s.end()-1);
      a.append("[&");
      // Python macro - cant do anything here to remove gcc warning
      PyObject* keys = PyMapping_Keys(atrs);
      PyObject* vals = PyMapping_Values(atrs);
      uint nkeys = PySequence_Size(keys);
      for(uint k = 0; k < nkeys; ++k) {
	PyObject* const name = PySequence_GetItem(keys, k);
	PyObject* const value = PySequence_GetItem(vals, k);
	if( k > 0 ) {
	  a.append(",");
	}
	a.append(strOfObject(name)).append("=").append(strOfObject(value));
      }
      a.append("]");
      
      Py_DECREF(keys);
      Py_DECREF(vals);
      Py_XDECREF(atrs);
    }
  }

  PyObject* const br = nd.branchlength;
  
  if( ! topoOnly && br != Py_None && includeStem ) {
    char* const b = PyOS_double_to_string(PyFloat_AsDouble(br), 'r', 0, Py_DTSF_ADD_DOT_0, 0);
    (s.end()-1)->append(":").append(b);
    PyMem_Free(b);
  }
}

PyObject*
TreeObject::toNewick(int nodeId, bool const topoOnly,
		     bool const includeStem, bool const withAttributes) const
{
  string s;
  if( ! treeNodes ) {
    tr->toNewick(s, nodeId, topoOnly, includeStem, withAttributes);
  } else {
    // go through node, in case node was changed
    if( nodeId == -1 ) {
      nodeId = tr->getRootID();
    }
    vector<string> c;
    tostr(c, nodeId, topoOnly, includeStem, withAttributes);
    s.assign(*c.begin());
  }
  return PyString_FromString(s.c_str());
}

static PyMemberDef treesSet_members[] = {
    {NULL}  /* Sentinel */
};

static PyObject*
treesSet_filterTaxa(TreesSetObject* self, PyObject* args)
{
  PyObject* pTaxaSeq;
  if( !PyArg_ParseTuple(args, "O", &pTaxaSeq) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  if( ! PySequence_Check(pTaxaSeq) || PySequence_Size(pTaxaSeq) == 0 ) {
    PyErr_SetString(PyExc_ValueError, "wrong args (sequences of taxa expected).") ;
    return 0;
  }

  TreesSet& ts = *self->ts;

  if( ts.store ) {
    PyErr_SetString(PyExc_ValueError, "Sorry, not implemeted for 'store'.") ;
    return 0;
  }
  
  vector<uint> taxaIndices;

  int nt = PySequence_Size(pTaxaSeq);
  for(int k = 0; k < nt; ++k) {
    PyObject* const s1 = PySequence_GetItem(pTaxaSeq, k);
    if( ! PyString_Check(s1) ) {
      PyErr_SetString(PyExc_ValueError, "wrong args (sequences of taxa expected).") ;
      return 0;
    }
      
    const char* const s1c = PyString_AS_STRING(s1);
    int const taxIndex = ts.hasTaxon(s1c);
    if( taxIndex < 0 ) {
      PyErr_Format(PyExc_ValueError, "Unknown taxon (%s).", s1c) ;
      return 0;
    }
    taxaIndices.push_back(taxIndex);
  }

  TreesSet* const nts = new TreesSet(ts.compressed, ts.precision, ts.store);
  
  for(uint k = 0; k < ts.nTrees(); ++k) {
    nts->add(ts, k, taxaIndices);
  }

  PyTypeObject* const type = self->ob_type;
  TreesSetObject* const n = static_cast<TreesSetObject *>(TreesSet_new(type, 0, 0));
  n->ts = nts;
  return n;
}

static PyObject*
treesSet_treei(TreesSetObject* self, PyObject* args)
{
  int n;
  if( !PyArg_ParseTuple(args, "i", &n) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  return TreesSet_getItemInternals(self, n);
}

static PyObject*
treesSet_add(TreesSetObject* self, PyObject* args, PyObject* kwds)
{
  const char* treeTxt;

  if( !PyArg_ParseTuple(args, "s", &treeTxt) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  int const k = self->ts->add(treeTxt, kwds);
  if( k < 0 ) {
    return 0;
  }
  
  return PyInt_FromLong(k);
}

// static PyObject*
// treesSet_dat(TreesSetObject* self)
// {
//   PyObject* t = PyTuple_New(0);
//   PyObject* d = PyDict_New();
//   PyObject* p = ((PyTypeObject*)self->data)->tp_new(self->data, t, d);

//   int x = ((PyTypeObject*)self->data)->tp_init(p, t, d);
  
//   Py_DECREF(t);
//   Py_DECREF(d);
//   return p;
// }

// PyObject*
// TreesSetObject::newNodeData(const char* tx, const double* branch, const double *height)
// {
//   PyObject* t = PyTuple_New(0);
//   PyObject* d = PyDict_New();
//   PyObject* p = ((PyTypeObject*)data)->tp_new(data, t, d);

//   if( tx ) {
//     PyDict_SetItemString(d, "taxon", PyString_FromString(tx));
//   }
//   if( branch ) {
//     PyDict_SetItemString(d, "branchlength", PyFloat_FromDouble(*branch));
//   }
//   if( height ) {
//     PyDict_SetItemString(d, "height", PyFloat_FromDouble(*height));
//   }
    
//   int x = ((PyTypeObject*)data)->tp_init(p, t, d);
  
//   Py_DECREF(t);
//   Py_DECREF(d);
//   return p;
// }

static PyMethodDef treesSet_methods[] = {
  {"add", (PyCFunction)treesSet_add, METH_VARARGS|METH_KEYWORDS,
   "Add a tree to set."
  },

  {"filterTaxa", (PyCFunction)treesSet_filterTaxa, METH_VARARGS,
   "Clone set while removing the given taxa list from each tree."
  },

  {"treei", (PyCFunction)treesSet_treei, METH_VARARGS,
   "Internals of tree (debugging)."
  },

  {NULL}  /* Sentinel */
};


static PySequenceMethods TreesSet_sequence_methods = {
  (lenfunc)TreesSet_len,                  /* sq_length */
  NULL,
  NULL,
  (ssizeargfunc)TreesSet_getItem,        /* sq_item */    
};

static PyTypeObject TreesSetType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "treesset.TreesSet",        /*tp_name*/
    sizeof(TreesSetObject),    /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)TreesSet_dealloc,  /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    &TreesSet_sequence_methods, /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash*/
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "TreesSet objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    treesSet_methods,             /* tp_methods */
    treesSet_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)TreesSet_init,      /* tp_init */
    0,                         /* tp_alloc */
    TreesSet_new,                 /* tp_new */
};

static PyObject*
parseTree(PyObject*, PyObject* args)
{
  const char* treeTxt;

  if( !PyArg_ParseTuple(args, "s", &treeTxt) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  vector<ParsedTreeNode> nodes;

  int const txtLen = strlen(treeTxt);
  int nc = readSubTree(treeTxt, nodes);

  if( nc > 0 ) {
    nc += skipSpaces(treeTxt + nc);
  }
  
  if( ! (nc == txtLen || (nc+1 == txtLen && treeTxt[nc] == ';')) ) {
    if( nc < 0) {
      int const where = txtLen+(nc+1);
      PyErr_Format(PyExc_ValueError, "failed parsing around %d (%10.10s ...).", where, treeTxt+where);
    } else {
      PyErr_Format(PyExc_ValueError, "extraneous characters at tree end: '%s'",
		   string(treeTxt+nc,std::max(5,txtLen-nc)).c_str());
    }
    return 0;
  }

  PyObject* n = PyTuple_New(nodes.size());
  for(uint k = 0; k < nodes.size(); ++k) {
    PyTuple_SET_ITEM(n, k, nodes[k].asPyObject());
  }
  return n;
}

static PyMethodDef module_methods[] = {
  {"parsetree", parseTree, METH_VARARGS, ""},
  
  {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
inittreesset(void) 
{
  PyObject* m;

  PyTypeObject* t[] = {&TreesSetType, &TreeType, &TreeNodeType, &TreeNodeDataType};
  for(uint i = 0; i < sizeof(t)/sizeof(t[0]); ++i) {
    if (PyType_Ready(t[i]) < 0) {
      return;
    }
  }

  m = Py_InitModule3("treesset", module_methods,
		     "Trees Set storage.");

  if (m == NULL)
    return;

  Py_INCREF(&TreesSetType);
  PyModule_AddObject(m, "TreesSet", (PyObject*)&TreesSetType);

  Py_INCREF(&TreeType);
  PyModule_AddObject(m, "Tree", (PyObject*)&TreeType);

  Py_INCREF(&TreeNodeType);
  PyModule_AddObject(m, "Node", (PyObject*)&TreeNodeType);

  Py_INCREF(&TreeNodeDataType);
  PyModule_AddObject(m, "NodeData", (PyObject*)&TreeNodeDataType);

  // const char* cd = "class NodeData: pass";
  //   // "class NodeData(object):\n  def __init__(self,taxon=None,branchlength=0.0,support=None):\n" 
  //   // "    self.taxon=taxon\n    self.branchlength=branchlength\n    self.support=support\n";
  // PyObject* globals = PyEval_GetGlobals(); //PyDict_New();
  // PyObject* l = PyDict_New();
  // PyObject* x = PyRun_String(cd, Py_file_input, globals, l);
  // int h = 0;
}
