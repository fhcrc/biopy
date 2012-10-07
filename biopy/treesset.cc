#include <Python.h>
#include "structmember.h"

#include <cassert>

#include <string>
using std::string;
#include <algorithm>
#include <vector>
#include <stack>
#include <list>
#include <map>
#include <unordered_map>
using std::vector;
using std::stack;
using std::list;
using std::unordered_map;

#include <cmath>
#include <limits>

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
typedef list< std::pair<string,string> > Attributes;

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
  PyObject* psubs = PyList_New(ns);
  for(int k = 0; k < ns; ++k) { 
    PyList_SET_ITEM(psubs, k, PyLong_FromLong(sons[k]));
  }
  PyList_SET_ITEM(nodeData, 2, psubs);

  Py_INCREF(Py_None);
  PyList_SET_ITEM(nodeData, 3, Py_None);

  // PyObject* o = PyTuple_New(2);
  // PyTuple_SET_ITEM(o, 0, taxon);
  // PyTuple_SET_ITEM(o, 1, trimString(v));
  // vals.push_back(o);
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
    // a terminal
    const char* s = txt;
    while( ! isspace(*s) && *s != ':' && *s != '[' && *s != ','
	   && *s != '(' && *s != ')' && *s != ']' ) {
      ++s;
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
      if( isspace(*txt) || isdigit(*txt) || has(*txt, ":.+-Ee") ) {
	nodeTxt.append(txt, 1);
	txt += 1;
	eat += 1;
      } else {
	break;
      }
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
  if( *nTxt && *nTxt == ':' ) {
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
  } else {
#if PYINTERFACE 
    Py_INCREF(Py_None);
    branch = Py_None;
#endif
  }

#if PYINTERFACE 
  PyList_SET_ITEM(nodeData, 1, branch);
#endif
  nodes.push_back(nodeData);

  return eat;
}

template<typename T> class Packer {
public:
  // unpacked. can be
  virtual uint size(void) const = 0;
  virtual vector<T> const&  unpacked(void) const = 0;
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
  vector<T> const&  unpacked(bool& isPermanent) const { isPermanent = true; return vals; }
  
private:
  vector<T> vals;
};


class FixedIntPacker : public Packer<uint> {
public:
  FixedIntPacker(uint nBitsPerValue, vector<uint>::iterator from, vector<uint>::iterator to);
  ~FixedIntPacker() {
    delete [] bits;
  }

  vector<uint> const&  unpacked(void) const;
  vector<uint> const&  unpacked(bool& isPermanent) const;
  //get(vector<uint>& to) const;

  virtual uint size(void) const { return len; }

  uint const nBitsPerValue : 8;
  uint const len : 24;
private:
  static uint const usize = 8*sizeof(char);
  static vector<uint> temp;
  
  unsigned char* bits;
};

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
			       vector<uint>::iterator from,
			       vector<uint>::iterator to) :
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
  virtual bool isCladogram(void) const = 0;
  TreeRep(Packer<uint>* t) : ptopo(t) {}
  virtual ~TreeRep() { delete ptopo; }

  uint nTaxa(void) const { return ptopo->size(); }
  
  vector<uint> const& topology() const;
  vector<uint> const& topology(bool& permanent) const;
protected:
  Packer<uint>*   ptopo;
};

vector<uint> const&
TreeRep::topology(void) const
{
  return ptopo->unpacked();
}

vector<uint> const&
TreeRep::topology(bool& permanent) const
{ 
  return ptopo->unpacked(permanent);
}

class CladogramRep : public TreeRep {
public:
  CladogramRep(Packer<uint>* t, Packer<uint>* h);
  ~CladogramRep();
 
  bool isCladogram(void) const { return true; }
 
  vector<uint> const& heights() const { return pheights->unpacked(); }

private:
  Packer<uint>*   pheights;
};


CladogramRep::CladogramRep(Packer<uint>* t, Packer<uint>* h) :
  TreeRep(t),
  pheights(h)
{}

CladogramRep::~CladogramRep()
{
  delete pheights;
}

template<typename T>
class PhylogramRep : public TreeRep {
public:
  PhylogramRep(Packer<uint>* t, Packer<T>* h, Packer<T>* txh);
  ~PhylogramRep() {
    delete pheights;
    delete ptxheights;
  }

  bool isCladogram(void) const { return false; }
  
  vector<T> const& heights() const { return pheights->unpacked(); } 
  vector<T> const* txheights() const {
    return ptxheights ? &ptxheights->unpacked() : static_cast< vector<T>* >(0);
  } 

private:
  Packer<T>*	pheights;
  Packer<T>* 	ptxheights;
};

template<typename T>
PhylogramRep<T>::PhylogramRep(Packer<uint>* t, Packer<T>* h, Packer<T>* txh) :
  TreeRep(t),
  pheights(h),
  ptxheights(txh)
{}

class TreesSet;

class Tree {
public:
  Tree(TreesSet const& _ts, uint _nt) :
    ts(_ts),
    nt(_nt),
    internals(0),
    sonsBlockSave(0)
    {}

  ~Tree() {
    delete internals;
    delete [] sonsBlockSave;
  }
  
  vector<uint> const& topology() const;
  
  void getTerminals(vector<uint>& terms) const;
  
  uint getRootID(void) const {
    if( ! internals ) setup();
    return internals->size() - 1;
  }

  struct Expanded {
    Expanded(int itax,
	     uint nSons,
	     uint* sons,
	     double* branch,
	     double height,
	     int    prev);
    
    int itax;
    uint nSons;
    uint* sons;
    double* branch;
    double height;
    int    prev;
  };

  bool isCladogram(void) const;
  
  uint nNodes(void) const;
  
  Expanded const& getNode(uint n) const {
    if( ! internals ) setup();
    return (*internals)[n];
  }

  void toNewick(string& s, int nodeId, bool topoOnly, bool includeStem) const;

  TreesSet const& ts;
  uint const nt;
private:
  void setup() const;
  void tostr(vector<string>& s, uint nodeId, bool topoOnly, bool includeStem) const;
  
  uint rep2treeInternal(vector<Expanded>& nodes,
		      uint low, uint hi,
		      vector<uint> const& tax,
		      vector<double> const& htax,
		      vector<double> const& hs,
		      uint*& sonsBlock,
		      uint* curiScratch,
		      uint bleft) const;

  // only after setup  
  mutable vector<Expanded>* internals;
  mutable uint              nTaxa;
  mutable uint*		sonsBlockSave;
};


struct TreesSetObject;

struct TreeObject : PyObject {
  void del(void) {
    delete tr;
    Py_XDECREF(taxa);
  }

  void init(TreesSetObject* _ts = NULL, Tree* t = NULL) {
    tr = t;
    ts = _ts;
    taxa = 0;
  }

  PyObject* getTaxa(void);

  PyObject* getTerminals(void) const ;
  PyObject* allIds(void) const;
  uint      getRootID(void) { return tr->getRootID(); }
  PyObject* getNode(uint nt) const;

  PyObject* toNewick(int nodeId, bool topoOnly, bool includeStem) const;
  
private:
  Tree* 	  tr;
  TreesSetObject* ts;
  PyObject*       taxa;
};

static void
Tree_dealloc(TreeObject* self)
{
  self->del();
  self->ob_type->tp_free((PyObject*)self);
}

static TreeObject *
Tree_new(PyTypeObject* type, PyObject*, PyObject *)
{
  TreeObject* self = (TreeObject *)type->tp_alloc(type, 0);

  if (self != NULL) {
    self->init();
  }

  return self;
}

PyObject*
Tree_getattr(TreeObject *self, PyObject* aname)
{
  const char* const name = PyString_AS_STRING(aname);
  if( ! strcmp(name, "root") ) {
    return PyInt_FromLong(self->getRootID());
  }
  return PyObject_GenericGetAttr(self, aname);
}

PyObject*
Tree_str(TreeObject *self)
{
  return self->toNewick(-1, false, false);
}

PyObject *
tree_getTaxa(TreeObject* self)
{
  return self->getTaxa();
}

PyObject *
tree_getTerminals(TreeObject* self)
{
  return self->getTerminals();
}

PyObject *
tree_allIds(TreeObject* self)
{
  return self->allIds();
}

PyObject*
tree_node(TreeObject* self, PyObject* args)
{
  int n;
  if( !PyArg_ParseTuple(args, "i", &n) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  return self->getNode(n);
}

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
  
  return self->toNewick(nodeId, topoOnly, incS);
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
  {"node", (PyCFunction)tree_node, METH_VARARGS,
   "get tree node."
  },

  {"toNewick", (PyCFunction)tree_2newick, METH_VARARGS,
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
    0,                         /* tp_dictoffset */
    0/*(initproc)Tree_init*/,      /* tp_init */
    0,                         /* tp_alloc */
    (newfunc)Tree_new,         /* tp_new */
};



struct TreeNodeDataObject : PyObject {
  PyObject* taxon;
  PyObject* branchlength;
  PyObject* height;
  
};

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

  if (self != NULL) {
    self->taxon = NULL;
    self->branchlength = NULL;
  }

  return self;
}

static void
TreeNodeData_init(TreeNodeDataObject* self, const char* taxon, const double* branchlength, const double* height)
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
}

static void
TreeNodeData_dealloc(TreeNodeDataObject* self)
{
  Py_XINCREF(self->taxon);
  Py_XINCREF(self->branchlength);
  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject TreeNodeDataType = {
  PyObject_HEAD_INIT(NULL)
  0,				/* ob_size        */
  "treesset.NodeData",		/* tp_name        */
  sizeof(TreeNodeDataObject),		/* tp_basicsize   */
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
  0,				/* tp_getattro    */
  PyObject_GenericSetAttr,	/* tp_setattro    */
  0,				/* tp_as_buffer   */
  Py_TPFLAGS_DEFAULT,		/* tp_flags       */
  "tree Node data.",		/* tp_doc         */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  0,             			/* tp_methods */
  NodeData_members,                  /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  0/*(initproc)TreeNode_init*/,      /* tp_init */
  0,                         /* tp_alloc */
  (newfunc)TreeNodeData_new,         /* tp_new */
};


struct TreeNodeObject : PyObject {
  PyObject* prev;
  PyObject* succ;
  PyObject* data;
  /*TreeNodeDataObject*/
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
  TreeNodeObject* self = (TreeNodeObject *)type->tp_alloc(type, 0);

  if (self != NULL) {
    self->succ = NULL;
    self->prev = NULL;
  }

  return self;
}

static void
TreeNode_init(TreeNodeObject* self, uint prev, uint nSons, uint* sons, PyObject* data)
// , const char* taxon const double* branch, const double* height)
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
  
  // TreeNodeDataObject* d = TreeNodeData_new(&TreeNodeDataType, 0, 0);
  // TreeNodeData_init(d, taxon, branch, height);
  // self->data = d;
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
  sizeof(TreeNodeObject),		/* tp_basicsize   */
  0,				/* tp_itemsize    */
  (destructor)TreeNode_dealloc,		/* tp_dealloc     */
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
  "Simple objects are simple.",	/* tp_doc         */
  0,		               /* tp_traverse */
  0,		               /* tp_clear */
  0,		               /* tp_richcompare */
  0,		               /* tp_weaklistoffset */
  0,		               /* tp_iter */
  0,		               /* tp_iternext */
  0,             			/* tp_methods */
  Node_members,                  /* tp_members */
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
  TreesSet(bool c, uint p, bool s) :
    compressed(c),
    store(s),
    precision(p)
    {}

  ~TreesSet() {
    for(auto t = trees.begin(); t != trees.end(); ++t) {
      delete *t;
    }
  }

  int add(const char* txt);

  string const& taxonString(uint k) const {
    assert(k < taxaList.size());
    return taxaList[k];
  }
  
  uint nTrees(void) const { return trees.size(); }
  
  TreeRep const& getTree(uint i) const {
    assert( i < nTrees() );
    return *trees[i];
  }
  
  void getHeights(uint nt, vector<double>& hs, vector<double>& txhs) const {
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
    
  bool const compressed : 8;
  bool const store : 8;
  uint const precision : 8;

  vector< vector<ParsedTreeNode> > asNodes;
  
private:
  TreeRep*		nodes2rep(vector<ParsedTreeNode>& nodes);
  
  vector<TreeRep*>	trees;
  
  unordered_map<string,uint>	taxaDict;
  vector<string>		taxaList;
  
  uint getTaxon(string const& taxon);
};


uint
TreesSet::getTaxon(string const& taxon)
{
  auto i = taxaDict.find(taxon);
  if( i == taxaDict.end() ) {
    int const k = taxaList.size();
    taxaList.push_back(taxon);
    taxaDict.insert( std::pair<string,uint>(taxon, k) );
    return k;
  }
  return i->second;
}

TreeRep*
TreesSet::nodes2rep(vector<ParsedTreeNode>& nodes)
{
  vector<uint> taxa;
  bool cladogram = true;
  uint maxTaxaIndex = 0;
  
  for(auto n = nodes.begin(); n != nodes.end() ; ++n) {
    if( n->taxon.length() ) {
      uint k = getTaxon(n->taxon);
      maxTaxaIndex = std::max(maxTaxaIndex, k);
      taxa.push_back(k);
    }
    if( n->branch ) {
      cladogram = false;
    }
  }
  
  vector<double> heights(taxa.size()-1);
  vector<int> locs(nodes.size());
  uint iloc = 0;
  vector<double>* taxaHeights = 0;
  
  for(auto n = nodes.begin(); n != nodes.end() ; ++n) {
    if( n->taxon.length() ) {
      if( ! n->branch ) {
	n->branch = new double(1);
      }
      assert( 0 <= iloc && iloc < locs.size() );
      
      locs[iloc] = iloc == 0 ? -1 : locs[iloc-1]+1;
      ++iloc;
    } else {
      double h = -1;
      for(auto s = n->sons.begin(); s != n->sons.end(); ++s) {
	h = std::max(*nodes[*s].branch, h);
      }
      if( ! cladogram ) {
	for(auto s = n->sons.begin(); s != n->sons.end(); ++s) {
	  double const h1 = *nodes[*s].branch;
	  double const dh = h - h1;
	  if( dh > 0 && !areSame(h,h1) ) {
	    if( ! taxaHeights ) {
	      taxaHeights = new vector<double>(taxa.size());
	      std::fill(taxaHeights->begin(), taxaHeights->end(), 0.0);
	    }
	    std::stack<uint> hp;
	    hp.push(*s); 
	    while( ! hp.empty() ) {
	      uint const x = hp.top(); hp.pop();
	      if( nodes[x].taxon.length() ) {
		(*taxaHeights)[locs[x]+1] += dh;
	      } else {
		heights[locs[x]] += dh;
		//assert( nodes[x].branch ) ;
		//*(nodes[x].branch) += dh;
		auto const& b = nodes[x].sons;
		for(auto a = b.begin(); a != b.end(); ++a) {
		  hp.push(*a);
		}
	      }
	    }
	  }
	}
      }
      for(auto s = n->sons.begin(); s != n->sons.end()-1; ++s) {
	uint const i = locs[*s]+1;            assert( 0 <= i && i < heights.size() );
	heights[i] = h;
      }
      assert( 0 < iloc && iloc < locs.size() );
      
      locs[iloc] = locs[iloc-1]; ++iloc;
      h += (n->branch ? *n->branch : 1);
      if( n->branch ) {
	*n->branch = h;
      } else {
	n->branch = new double(h);
      }
    }
  }

  Packer<uint>* top = 0;
  
  if( compressed ) {
    int nbitsStoreTaxa = lg2i(maxTaxaIndex) + 1;
    top = new FixedIntPacker(nbitsStoreTaxa, taxa.begin(), taxa.end());
  } else {
    top = new SimplePacker<uint>(taxa);
  }
    
  if( cladogram ) {
    vector<uint> hs(heights.size());
    for(uint i = 0; i < hs.size(); ++i) {
      hs[i] = static_cast<uint>(heights[i] + 0.5);
    }
    Packer<uint>* hsb;
    // single leaf trees has no heights
    if( compressed && heights.size() > 0 ) {
      int nbitsStoreH = lg2i(*std::max_element(heights.begin(),heights.end())) + 1;
      hsb = new FixedIntPacker(nbitsStoreH, hs.begin(), hs.end());
    } else {
      hsb = new SimplePacker<uint>(hs);
    }
    return new CladogramRep(top, hsb);
  } else {
    if( precision == 8 ) {
      SimplePacker<double>* hsb = new SimplePacker<double>(heights);
      SimplePacker<double>* txhs = 0;
      if( taxaHeights ) {
	txhs = new SimplePacker<double>(*taxaHeights);
      }
      return new PhylogramRep<double>(top, hsb, txhs);
    } else {
      SimplePacker<float>* hsb = new SimplePacker<float>(heights);
      SimplePacker<float>* txhs = 0;
      if( taxaHeights ) {
	txhs = new SimplePacker<float>(*taxaHeights);
      }
      return new PhylogramRep<float>(top, hsb, txhs);
    }
  }
}

int
TreesSet::add(const char* treeTxt)
{
  vector<ParsedTreeNode> nodes;
  if( readSubTree(treeTxt, nodes) < 0 ) {
    PyErr_SetString(PyExc_ValueError, "failed parsing.") ;
    return -1;
  }

  if( store ) {
    asNodes.push_back(nodes);
    return asNodes.size()-1;
  } else {
    trees.push_back(nodes2rep(nodes));
    return trees.size()-1;
  }
}

Tree::Expanded::Expanded(int _itax,
			 uint _nSons,
			 uint* _sons,
			 double* _branch,
			 double _height,
			 int    _prev) :
  itax(_itax),
  nSons(_nSons),
  sons(_sons),
  branch(_branch),
  height(_height),
  prev(_prev)
{}
  

uint Tree::rep2treeInternal(vector<Expanded>& nodes,
			    uint low, uint hi,
			    vector<uint> const& tax,
			    vector<double> const& htax,
			    vector<double> const& hs,
			    uint*& sonsBlock,
			    uint* curiScratch,
			    uint bleft) const
{
  if( low == hi ) {
    nodes.push_back(Expanded(tax[low],0,0,0,htax[low],-1));
  } else {
    uint* curi;
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
    uint const nSplits = (uint)(curi-curiScratch);
    assert( nSplits <= bleft );
    bleft -= nSplits;
    
    uint* sons = sonsBlock;
    uint nSons = 1 + nSplits;
    sonsBlock += nSons;

    assert( bleft > 0 );
    *curi = hi; ++curi; --bleft;
    
    //uint xlow = low;
    for(uint* x = curiScratch; x < curi; ++x) {
      uint const k = rep2treeInternal(nodes, low, *x, tax, htax, hs, sonsBlock, curi, bleft);
      double const hs = nodes[k].height;
      nodes[k].branch = new double(curh - hs);
      *sons = k; ++sons;
      low = *x+1;
    }
    // uint k = rep2treeInternal(nodes, xlow, hi, tax, htax, hs, sonsBlock, curi, bleft);
    // double hs = nodes[k].height;
    // nodes[k].branch = new double(curh - hs);
    // *sons = k; ++sons;
    sons -= nSons;
    for(uint* s = sons; s < sons+nSons; ++s) {
      nodes[*s].prev = nodes.size();
    }
    nodes.push_back(Expanded(-1, nSons,sons,0,curh,-1));
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
  vector<uint> const& t = rep.topology(isp);
  vector<uint> const& tax = isp ? t : vector<uint>(t);
  
  nTaxa = tax.size();
  ts.getHeights(nt, hs, txhs);
  if( txhs.size() == 0 ) {
    txhs.resize(nTaxa, 0.0);
  }

  internals = new vector<Expanded>;
  sonsBlockSave = new uint[2*nTaxa];  // allocate in one contiguous block 
  uint block[2*nTaxa];     // scratch only
  uint* s = sonsBlockSave; // keep sonsBlockSave safe
  rep2treeInternal(*internals, 0, hs.size(), tax, txhs, hs, s, block, 2*nTaxa);

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
Tree::topology() const
{
  return ts.getTree(nt).topology();
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
Tree::tostr(vector<string>& s, uint nodeId, bool topoOnly, bool includeStem) const
{
  Expanded const& n = (*internals)[nodeId];
  if( n.itax >= 0 ) {
    s.push_back(ts.taxonString(n.itax));
  } else {
    for(uint i = 0; i < n.nSons; ++i) {
      tostr(s, n.sons[i], topoOnly, true);
    }
    auto first = s.end() - n.nSons;
    std::sort(first, s.end());
    first->insert(first->begin(),'(');
    for(uint i = 1; i < n.nSons; ++i) {
      first->append(",").append(*(first + i));
    }
    first->append(")");
    s.erase(first+1, s.end());
  }

  if( ! topoOnly && n.branch && includeStem ) {
    char* const b = PyOS_double_to_string(*n.branch, 'r', 0, Py_DTSF_ADD_DOT_0, 0);
    (s.end()-1)->append(":").append(b);
    PyMem_Free(b);
  }
}

void
Tree::toNewick(string& s, int nodeId, bool topoOnly, bool includeStem) const
{
  if( nodeId == -1 ) {
    nodeId = getRootID();
  } else {
    setup();
  }
  vector<string> c;
  Tree::tostr(c, static_cast<uint>(nodeId), topoOnly, includeStem);
  s.assign(*c.begin());
}

struct TreesSetObject : PyObject {
  TreesSet* ts;
  //PyTypeObject* data;
  
  void del(void) {
    delete ts;
    for(auto s = taxa->begin(); s != taxa->end(); ++s) {
      Py_XDECREF(*s);
    }
    delete taxa;
  }

  void init(void) {
    ts = NULL;
    taxa = new vector<PyObject*>();
  }

  PyObject* taxon(uint k) {
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

  PyObject* newNodeData(const char* tx, const double* branch, const double *height);

private:
  vector<PyObject*>* taxa;
};


static void
TreesSet_dealloc(TreesSetObject* self)
{
  self->del();
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
TreesSet_new(PyTypeObject* type, PyObject *args, PyObject *kwds)
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
  static const char *kwlist[] = {"compressed", "precision", "store"/*, "data"*/,
				 static_cast<const char*>(0)};
  PyObject* comp = 0;
  PyObject* sto = 0;
  //PyObject* data = 0;
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
  //self->data = (PyTypeObject*)data;
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

template<typename T>
PyObject*
vector2tuple(vector<T> const& v)
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
  PyObject* n = PyTuple_New(4);
  bool const isc = r.isCladogram();
  PyTuple_SET_ITEM(n, 0, PyBool_FromLong(isc));
  vector<uint> const& topo = r.topology();
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
      PyTuple_SET_ITEM(n, 2, vector2tuple(p.heights()));
      PyTuple_SET_ITEM(n, 3, vector2tuple(*p.txheights()));
    } else {
      typedef PhylogramRep<double> T;
      
      T const& p = static_cast<T const&>(r);
      PyTuple_SET_ITEM(n, 2, vector2tuple(p.heights()));
      PyTuple_SET_ITEM(n, 3, vector2tuple(*p.txheights()));
    }
  }  
  return n;
}

// PyObject*
// TreesSet_getItemx(TreesSetObject* self, Py_ssize_t i)
// {
//   return TreesSet_getItemInternals(self, i);
// }

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
  to->init(self, new Tree(*self->ts, nt));
  return to;
}



PyObject*
TreeObject::getTaxa(void)
{
  if( taxa == 0 ) {
    auto const& topo = tr->topology();
    taxa = PyTuple_New(topo.size());
    for(uint k = 0; k < topo.size(); ++k) {
      PyTuple_SET_ITEM(taxa, k, ts->taxon(topo[k]));
    }
    // we keep one refcount on creation, so this tuple stays
    // around until tree goes for efficiency.
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
  
  Tree::Expanded const& e = tr->getNode(nt);

  const char* tx = 0;
  if( e.itax >= 0 ) {
    tx = tr->ts.taxonString(e.itax).c_str();
  }
  bool const isc = tr->isCladogram();
  TreeNodeObject* node = TreeNode_new(&TreeNodeType, 0, 0);

  //TreeNode_init(node, tx, e.prev, e.nSons, e.sons, isc ? 0 : e.branch, isc ? 0 : &e.height);
  // PyObject* d = ts->newNodeData(tx, isc ? 0 : e.branch, isc ? 0 : &e.height); 
  // self->data = d;

  TreeNodeDataObject* d = TreeNodeData_new(&TreeNodeDataType, 0, 0);
  TreeNodeData_init(d, tx, isc ? 0 : e.branch, isc ? 0 : &e.height);
 
  TreeNode_init(node, e.prev, e.nSons, e.sons, d);
  
  return node;
}

PyObject*
TreeObject::toNewick(int nodeId, bool topoOnly, bool includeStem) const
{
  string s;
  tr->toNewick(s, nodeId, topoOnly, includeStem);
  return PyString_FromString(s.c_str());
  //  return 0;
}

static PyMemberDef treesSet_members[] = {
  // {"name", T_OBJECT_EX, offsetof(TreeObject, name), 0,
  //  "tree name"},
  // {"weight", T_OBJECT_EX, offsetof(TreeObject, weight), 0,
  //  "tree weight"},
  // {"rooted", T_OBJECT_EX, offsetof(TreeObject, rooted), 0,
  //  "rooted/unrooted"},
    {NULL}  /* Sentinel */
};

static PyObject *
treesSet_treei(TreesSetObject* self, PyObject* args)
{
  int n;
  if( !PyArg_ParseTuple(args, "i", &n) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }
  return TreesSet_getItemInternals(self, n);
}

static PyObject *
treesSet_add(TreesSetObject* self, PyObject* args)
{
  const char* treeTxt;

  if( !PyArg_ParseTuple(args, "s", &treeTxt) ) {
    PyErr_SetString(PyExc_ValueError, "wrong args.") ;
    return 0;
  }

  int const k = self->ts->add(treeTxt);
  if( k < 0 ) {
    return 0;
  }
  
  return PyLong_FromLong(k);
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
     {"add", (PyCFunction)treesSet_add, METH_VARARGS,
      "Add a tree to set."
     },

     {"treei", (PyCFunction)treesSet_treei, METH_VARARGS,
      "Internals of tree."
     },

     // {"dat", (PyCFunction)treesSet_dat, METH_NOARGS,
     //  "Add a tree to set."
     // },
     
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

static PyMethodDef module_methods[] = {
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
  PyModule_AddObject(m, "TreesSet", (PyObject *)&TreesSetType);

  Py_INCREF(&TreeType);
  PyModule_AddObject(m, "Tree", (PyObject *)&TreeType);

  Py_INCREF(&TreeNodeType);
  PyModule_AddObject(m, "Node", (PyObject *)&TreeNodeType);

  Py_INCREF(&TreeNodeDataType);
  PyModule_AddObject(m, "NodeData", (PyObject *)&TreeNodeDataType);

  // const char* cd = "class NodeData: pass";
  //   // "class NodeData(object):\n  def __init__(self,taxon=None,branchlength=0.0,support=None):\n" 
  //   // "    self.taxon=taxon\n    self.branchlength=branchlength\n    self.support=support\n";
  // PyObject* globals = PyEval_GetGlobals(); //PyDict_New();
  // PyObject* l = PyDict_New();
  // PyObject* x = PyRun_String(cd, Py_file_input, globals, l);
  // int h = 0;
}
