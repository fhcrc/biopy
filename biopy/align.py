## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

import sys
import calign
from collections import namedtuple

# messy code, still in various states of development

__all__ = ["parseScores", "MatchScores", "seqMultiAlign", "mpc", "mpa", "orderSeqsbyTree", "cons"]

MatchScores = namedtuple('MatchScores', "match mismatch gap gape freeEnds")
defaultMatchScores = MatchScores(match = 10, mismatch = -5, gap = -6, gape = -6,
                                 freeEnds = True)

_tab = [None]*5
for x in 'AGCTN' :
  _tab[getattr(calign, x)] = x 
_tab = "".join(_tab)
_gtab = _tab + '-';                  assert _gtab[calign.GAP] == '-'

def parseScores(strScores) :
  if strScores is None :
    return defaultMatchScores
  
  m = strScores.split(',')
  if not 4 <= len(m) <= 5 :
    raise RuntimeError("Expecting 4 or 5 comma separated values")
  if len(m) == 4 :
    m.insert(3,m[2])
  ms = [float(x) if x != "" else dflt for x,dflt in zip(m,defaultMatchScores)[:4]]

  f = bool(int(m[4])) if m[4] else defaultMatchScores.freeEnds
  ms.append(f)
  
  if not (ms[0] > 0 and ms[1] <= 0 and ms[2] <= 0 and ms[3] <= 0) :
    raise RuntimeError("Expecting a non positive mismatch score and gap penalty.")
  return MatchScores(*ms)

def ntoi(n) :
  i = _tab.find(n)
  if i == -1 :
    i = calign.N
  return i

def ngtoi(n) :
  i = _gtab.find(n.upper())
  if i == -1 :
    return calign.N
  return i

def iton(s) :
  ## ints to string
  if isinstance(s, str) :
    return s
  if isinstance(s, (int,long)) :
    return _gtab[s]
  if isinstance(s, (tuple,list)) :
    return "".join([_gtab[x] for x in s])

def sasn(seq) :
  if isinstance(seq[0], str) :
    seq = tuple(ngtoi(x) for x in seq)
  return seq

def stripseq(s) :
  if isinstance(s[0], (int,long)) :
    return [x for x in s if x != calign.GAP]
  return "".join([x for x in s if x != '-'])

#matchScore = 10
#misMatchScore = -5
#_gapPenalty = -6

def score(nSite, n, cnts) :
  _gapPenalty = defaultMatchScores.gap
  c = cnts[nSite]
  if n == calign.GAP :
    return _gapPenalty * (sum(c) - c[calign.GAP])
  return matchScore * c[n] + misMatchScore * (sum(c[:4]) - c[n]) + _gapPenalty * c[calign.GAP]

def ascore(s, cnts) :
  s = sasn(s)
  return sum([score(k, n, cnts) for k,n in enumerate(s)])

def cons(profile) :
  return iton([p.index(max(p)) for p in profile])
  
def createProfile(seqs) :
  return calign.createProfile(seqs)

def alignToProfile(seq, profile) :
  wasNucs = isinstance(seq[0], str)
  s = calign.profileAlign(seq, profile)
  s = iton(s) if wasNucs else s
  
  return s

def _sbreak(s, n = 100) :
  while len(s) :
    print s[:n]
    s = s[n:]

def _amc(x,y) :
  if x==y and x!='-' and x!= ' ' :
    return ':'
  if x!=y and x!='-' and y!='-' :
    return '*'
  return ' '

def showAlignment(s1, s2, n = 100, stripGaps = True) :
  s1 = iton(s1)
  s2 = iton(s2)

  if stripGaps:
    s1,s2 = zip(*[(x,y) for x,y in zip(s1,s2) if not (x=='-' and y=='-')])
    s1,s2 = ["".join(x) for x in s1,s2]
 
  stats = sum([x==y and x !='-' and y !='-' for x,y in zip(s1,s2)]),\
          sum([x!=y and x !='-' and y !='-' and x != 'N' and y != 'N' for x,y in zip(s1,s2)]),\
          sum([ x =='-' or y =='-' for x,y in zip(s1,s2)])
   
  while s1 or s2 :
    s1n,s2n = s1[:n],s2[:n]
    if len(s1n) < len(s2n) :
      s1n = s1n + ' '*(len(s2n)-len(s1n))
    if len(s2n) < len(s1n) :
      s2n = s2n + ' '*(len(s1n)-len(s2n))
    print s1n
    print "".join([_amc(x,y) for x,y in zip(s1n,s2n)])
    print s2n,'\n'
    s1 = s1[n:]
    s2 = s2[n:]
  
  print stats

def mshowAlignment(sqs, n = 100, stripGaps = True) :
  sqs = [iton(x) for x in sqs]
  ns = len(sqs[0])
  if stripGaps:
    allg = []
    for i in range(ns) :
      if all([s[i]=='-' for s in sqs]) :
        allg.append(i)
    for i in reversed(allg) :
      sqs = [s[:i] + s[i+i:] for s in sqs]
   
  while any(sqs) :
    sqsn = [s[:n] if len(s) >= n else s + ' '*(n-len(s))  for s in sqs]
    for i in range(len(sqs)-1) :
      print sqsn[i]
      print "".join([_amc(x,y) for x,y in zip(sqsn[i],sqsn[i+1])])
    print sqsn[-1]
    print
    sqs = [s[n:] for s in sqs]
  
## def mp(seqs) :
##   a = calign.globalAlign(seqs[0], seqs[1], freeEndGaps=1)
##   for ss in seqs[2:] :
##     p = createProfile(a)
##     s2 = [ngtoi(x) for x in ss]
##     x = alignToProfile(s2, p + [[0,0,0,0,0,sum(p[0])],]*20)
##     i = -1
##     while x[i] == calign.GAP :
##       i -= 1
##     n = max(len(x)+(i+1),len(p))
##     x = x[:n]
##     if n > len(p) :
##       a = [x + (calign.GAP,)*(n-len(p)) for x in a]
##     a = a + (x,)
##   return a

def seqMultiAlign(seqs, scores = defaultMatchScores, report=False) :
  if len(seqs) < 2:
    return seqs
  
  a = calign.globalAlign(seqs[0], seqs[1], scores=scores)
  ns = 2
  p = calign.createProfile(a)

  a = tuple(iton(x) for x in a)
  
  for kk,s2 in enumerate(seqs[2:]) :
    #print ns
    # assert p == calign.createProfile(a)
    assert len(a[0]) == len(p) and \
           p[0][calign.GAP] != ns and p[-1][calign.GAP] != ns
    
    pad = 20
    if len(p)+2*pad < len(s2) :
      # enough for sequences start to align
      pad = (len(s2) - len(p))  ;                assert len(p)+2*pad >= len(s2)

    pa,extendLeft,extendRight = calign.profileAlign(s2, p, pad=pad, chop=True,
                                                    gapPenalty=defaultMatchScores.gap)

    if extendLeft > 0 or extendRight > 0 :
      gapProfile = [0,0,0,0,0,ns]
      p1 = tuple((list(gapProfile) for k in range(extendLeft))) + p \
           + tuple((list(gapProfile) for k in range(extendRight)))
    else :
      p1 = p
    
    for k,n in enumerate(pa) :
      p1[k][n] += 1
    p = p1
      
    if extendLeft > 0 or extendRight > 0 :
      fr = '-'*extendLeft
      bk = '-'*extendRight
      a = tuple((fr + x + bk for x in a))
      
    a = a + (iton(pa),)
    ns += 1

    if report and (kk+1) % 1000 == 0 :
      import sys
      print kk+1, len(a[0]),
      sys.stdout.flush()
  if report: print
  
  return a


from math import log
from itertools import count

def sortByEntropy(seqs) :
  alen = len(seqs[0])
  trans = not (seqs[0][0] in range(6))
  
  #print "sorting,"
  counts = calign.createProfile(seqs)
  
  #print "stage 1,"
  cc = [sum([log(x) * x for x in c if x > 0])/len(seqs) - log(len(seqs))
        for c in counts]
  scc = sorted(zip(cc, count()), reverse=0)
  occ = [x[1] for x in scc]
  #print "stage 2,"

  if trans :
    sx = sorted(zip([[ngtoi(x[k]) for k in occ] for x in seqs], count()))
  else :
    sx = sorted(zip([[x[k] for k in occ] for x in seqs], count()))
    
  #print "stage 3,"
  seqsSorted = [k for s,k in sx]
  #print "sorted"
  
  return seqsSorted


import os.path
def saveAlignment(a, names, fname, sortit = False, overw = False) :
  assert overw or not os.path.exists(fname)

  fs = file(fname, 'w')
  if sortit:
    a1s = sortByEntropy(a)
    
    for i in a1s :
      print >> fs, names[i]
      print >> fs, iton(a[i])
  else :
    for i in range(len(a)):
      print >> fs, names[i]
      print >> fs, iton(a[i])
  fs.close()

def toRanges(r) :
  if len(r) :
    r = sorted(r)
    r1 = [ [r[0],r[0]+1] ]
    for i in r[1:]:
      if r1[-1][1] == i :
        r1[-1][1] += 1
      else :
        r1.append( [i,i+1] )
    return r1
  return []

def removeColumns(l, r) :
  r1 = toRanges(r)
  for a,b in r1[::-1] :
    l = l[:a] + l[b:]
  return l,r1

def restoreColumns(l, r, fill) :
  for a,b in r :
    l = l[:a] + [fill,]*(b-a) + l[a:]
  return l
  

def refineAlignment(al, ci = [0], drop = False, mx = -1, rev = False, verbose = False) :
  if not isinstance(al, list) :
    al = list(al)
    
  p = calign.createProfile(al)
  if rev:
    can = [(k,s) for k,s in enumerate(al) if all([s[x] == '-' for x in ci])]
  else :
    can = [(k,s) for k,s in enumerate(al) if any([s[x] != '-' for x in ci])]

  if verbose:
    print len(al[0]), len(can)
    import sys
    sys.stdout.flush()
    
  if mx > 0 and len(can) > mx:
    return
  
  changed = 0
  if len(can) :
    
    pr = calign.createProfile([s for n,s in can])
    p2 = tuple([[a-b for a,b in zip(x,y)] for x,y in zip(p,pr)])
    if drop:
      p2,rx = removeColumns(p2, ci)
    
    for n,s in can:
      if len(s.replace('-','')) <= len(p2) :
        ra = calign.profileAlign(s, p2)
        if drop :
          ra = restoreColumns(list(ra),rx, calign.GAP)
        al[n] = iton(ra)
        if tuple(ra) != sasn(s) :
          changed += 1

  if changed :
    p = calign.createProfile(al)
    r = [k for k,i in enumerate(p) if i[calign.GAP] == len(al)]

    r1 = toRanges(r)
    if r1 :
      for a,b in r1[::-1] :
        for k in range(len(al)) :
          x = al[k]
          al[k] = x[:a] + x[b:]

  return al,len(al[0]), changed


def refineSingle(al, gapPenalty = defaultMatchScores.gap) :
  if not isinstance(al, list) :
    al = list(al)
    
  p = calign.createProfile(al)

  q = []
  for k,s in enumerate(al) :
    ps = calign.createProfile([s])
    p2 = tuple([[a-b for a,b in zip(x,y)] for x,y in zip(p,ps)])

    q.append( calign.profileAlign(s, p2, gapPenalty = gapPenalty) )

  return q

def mpc(seqs, nRefines = 4, gapPenalty = defaultMatchScores.gap) :
  al = seqMultiAlign(seqs, scores = defaultMatchScores._replace(gap = gapPenalty))

  c0 = stripseq(cons(calign.createProfile(al)))
  r = refineSingle(al, gapPenalty = gapPenalty)
  c1 = stripseq(cons(calign.createProfile(r)))
  cnt = 0
  while c0 != c1 and cnt < nRefines:
    c0 = c1
    r = refineSingle(r, gapPenalty = gapPenalty)
    c1 = stripseq(cons(calign.createProfile(r)))
    cnt += 1
  return c1, r

# muscleExec = "/home/joseph/bin/muscle"
import StringIO,  tempfile, os
from bioutils import readFasta

def writeAlignmentAsNexus(al, nm, fs, tr = None) :
  print >> fs, '#NEXUS'
  print >> fs, """Begin data;
  Dimensions ntax=%d nchar=%d;
  Format datatype=dna symbols="ACTG" missing=? gap=-;
  Matrix
  """ % (len(al),len(al[0][1]))

  for h,b in al :
    print >> fs,h[1:],b
  print >> fs, ';'
  print >> fs, 'End;'

  if tr :
    print >> fs, """Begin Trees;
  Tree %s =""" % nm
    print >> fs, str(tr),';'
    print >> fs, 'End;'

  
def musAl(tr, getSeqs, muscleExec, outTo = None,
          tmpDir = None, reportProgress = False) :
  if tmpDir is None :
    tmpDir = tempfile.gettempdir()
    keep = False
  else :
    keep = True

  nm = tr.name if tr.name else "name"
  ftin = tmpDir + '/mu%s.tre' % nm
  fs = file(ftin, 'w')
  print >> fs, str(tr),';'
  fs.close()

  ffin = tmpDir + '/mu%s.fas' % nm
  fs = file(ffin , 'w')
  for x in tr.get_terminals() :
    i = int(tr.node(x).data.taxon)
    print >> fs, '>' + str(i)
    print >> fs, getSeqs([i])[0]
  fs.close()

  fout = tmpDir + '/mu%s.afa' % nm
  pargs = [muscleExec,'-in', ffin,'-out', fout,
           '-usetree_nowarn',ftin, '-maxiters','3']
  #return pargs
  
  import subprocess
  p = subprocess.Popen(pargs, stdout = subprocess.PIPE,
                       stderr = subprocess.PIPE)
  for x in p.stderr :
    if reportProgress :
      print >> sys.stderr,"*",
    sys.stderr.flush()
  
  #  while next(p.stderr) :
  #  pass
  returncode = p.wait()
  if returncode != 0 :
    raise returncode

  al = [(h,b) for h,b in readFasta(file(fout))]
  if outTo is not None :
    writeAlignmentAsNexus(al, nm, outTo, tr)

  if not keep :
    for f in (fout, ftin, ffin) :
      os.remove(f)

  if outTo is None :
    return al

def geto(tr, i) :
  n = tr.node(i)
  if not n.succ:
    return [n]
  else :
    lf, rt = geto(tr, n.succ[0]) , geto(tr, n.succ[1])
    if len(lf) < len(rt) :
      return rt + lf
    return lf + rt
    #return geto(tr, n.succ[0]) + geto(tr, n.succ[1])
  
def orderSeqsbyTree(tr, seqs) :
  dseqs = dict(seqs)
  o = geto(tr, tr.root)
  s = [dseqs[n.data.taxon.strip("'")] for n in o]
  return s

def trmsa(tr, seqs) :
  s = orderSeqsbyTree(tr, seqs)
  return seqMultiAlign(s)


from biopy.treeutils import getPostOrder

def trimendsp(p, th = .6) :
  n = int(sum(p[0])*th)
  
  for i in range(len(p)-1, 1, -1) :
    if p[i][calign.GAP] < n :
      break
  for j in range(len(p)-1) :
    if p[j][calign.GAP] < n :
      break
  #print len(p), len(p[j:i+1])
  return p[j:i+1]

def mpa(tr, seqs, scores = defaultMatchScores, trimEnd = None) :
  dseqs = dict(seqs)
  #scores = (None,None,gapPenalty,feg)
  for n in getPostOrder(tr) :
    data = n.data
    if not n.succ :
      data.seq = (None,dseqs[n.data.taxon.strip("'")])
    else :
      s1,s2 = [tr.node(x).data.seq for x in n.succ]
      if s1[1] :
        if s2[1] :
          a = calign.globalAlign(s1[1],s2[1], scores = scores)
          data.seq = (calign.createProfile(a),None)
        else :
          p1,p2 = calign.createProfile(s1[1:]), s2[0]
          pa = calign.prof2profAlign(p1,p2, scores = scores)
          data.seq = (trimendsp(pa, trimEnd) if trimEnd is not None else pa,None)
          #print len(pa)
      else :
        p1 = s1[0]
        if s2[1] :
          p2 = calign.createProfile(s2[1:])
        else :
          p2 = s2[0]
        pa = calign.prof2profAlign(p1,p2, scores = scores)
        data.seq = (trimendsp(pa, trimEnd) if trimEnd is not None else pa,None)
        #print len(pa)
        #import pdb; pdb.set_trace()
  assert n.id == tr.root
  return n.data.seq[0]


## def trimp(p) :
##   n = sum(p[0])*8//10
  
##   for i in range(len(p)-1, 1, -1) :
##     if p[i][calign.GAP] < n :
##       break
##   for j in range(len(p)-1) :
##     if p[j][calign.GAP] < n :
##       break
##   return p[j:i+1]

def trimp(p, th = .6) :
  n = int(sum(p[0])*th)

  np = []
  for x in p :
    if x[calign.GAP] < n :
      np.append(x)
  return np

def mpat(tr, seqs) :
  dseqs = dict(seqs)
  for n in getPostOrder(tr) :
    data = n.data
    if not n.succ :
      data.seq = (None,dseqs[n.data.taxon])
    else :
      s1,s2 = [tr.node(x).data.seq for x in n.succ]
      if s1[1] :
        if s2[1] :
          a = calign.globalAlign(s1[1],s2[1])
          data.seq = (calign.createProfile(a),None)
        else :
          p1,p2 = calign.createProfile(s1[1:]), s2[0]
          #assert all([sum(x)==sum(p1[0]) for x in p1])
          #assert all([sum(x)==sum(p2[0]) for x in p2])
          pa = calign.prof2profAlign(p1,p2)
          data.seq = (trimp(pa),None)
          #print len(pa)
      else :
        p1 = s1[0]
        if s2[1] :
          p2 = calign.createProfile(s2[1:])
        else :
          p2 = s2[0]
        #assert all([sum(x)==sum(p1[0]) for x in p1])
        #assert all([sum(x)==sum(p2[0]) for x in p2])
        pa = calign.prof2profAlign(p1,p2)
        data.seq = (trimp(pa),None)
        #print len(pa)
        #import pdb; pdb.set_trace()
  assert n.id == tr.root
  return n.data.seq[0]
