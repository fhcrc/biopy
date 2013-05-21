## This file is part of biopy.
## Copyright (C) 2013 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

from __future__ import division

import time, gzip, random, math

from collections import defaultdict, namedtuple
from itertools import ifilter, count

from genericutils import tohms, fileFromName

import calign, cclust
import align

from parseNewick import parseNewick
from treeutils import CAhelper, getPostOrder

# for hierarchy.to_tree, can we get rid of that??
import scipy.cluster

__all__ = ["findDuplicates", "deClutter", "deClutterDown", "treeFromSeqs",
           "MatchScores", "defaultMatchScores", "declutterToTrees"
           "saveDistancesMatrix", "getDistanceMatrix",
           "clusterFromTree", "assembleTree"]

MatchScores = namedtuple('MatchScores', "match mismatch gap freeEnds")
defaultMatchScores = MatchScores(match = 10, mismatch = -5, gap = -6, freeEnds = True)
                                      
def nPairs(x) :
  if isinstance(x, (list,tuple)) :
    x = len(x)
  return (x*(x-1))//2

def saveDistancesMatrix(fname, dsm, labels, compress=-1) :
  if fname :
    suf = '.dists.gz' if compress else '.dists'
    if not fname.endswith(suf):
      fname = fname + suf

    assert nPairs(labels) == len(dsm)
    if compress >= 0 :
      fs = gzip.open(fname, 'wb', max(compress, 9))
    else :
      fs = file(fname, 'w')
    print >> fs, '#',' '.join(labels)
    for x in dsm :
      print >> fs, x
    fs.close()

def getDistanceMatrix(saveName) :
  fl = fileFromName(saveName)
  labs = [x.strip() for x in next(fl)[2:].split(' ')]
  dists = [float(x) for x in fl]
  fl.close()
  return dists,labs

def findDuplicates(allSeqs, verbose = None) :
  """ Locate sequences s1,s2 such that s1 is a prefix of s2.
  return the merge of those pairs."""
  
  if verbose:
    tstart = time.clock()

  byLen = defaultdict(lambda : [])
  for ns,s in enumerate(allSeqs):
    byLen[len(s)].append(ns)
  cc = []
  for b in byLen:
    d = defaultdict(lambda : [])
    for si in byLen[b] :
      d[allSeqs[si]].append(si)
    for s,si in ifilter(lambda x : len(x[1])>1 , d.iteritems()):
      cc.append(si)
  
  if verbose:
    print >> verbose, "find duplicates in ", tohms(time.clock() - tstart)
    
  return cc

def mergeGroupings(grps) :
  """ Merge any partitions with a non empty intersetions. """
  # must be a more efficient way, but what?

  grps = sorted([sorted(x) for x in grps])
  
  newGrps = [set(x) for x in grps]
  for i in range(len(newGrps)-1, 0,-1) :
    merged = False
    for j in range(i-1,-1,-1) :
      if not newGrps[i].isdisjoint(newGrps[j]) :
        newGrps[j].update(newGrps[i])
        newGrps[i] = None
        break

  sGrps = sorted([sorted(x) for x in newGrps if x])
  return sGrps

# C code will be faster?
def _buildLookup(seqs, asSet=True) :
  bl = dict()
  
  for j,s in enumerate(seqs) :
    for i in range(len(s)-11) :
      p = s[i:i+11]
      m = bl.get(p)
      if m is not None :
        m.append(j)
      else :
        bl[p] = [j]
  matches = dict([(x,set(y) if asSet else y) for x,y in bl.iteritems() if len(y) > 1])
        
  return matches

import heapq

def hiter(heap) :
  while heap:
    yield heapq.heappop(heap)

_failsHardTH = 20

def getPotentials(iSeq, seqs, lseqs, th, fdis, matches, elem2grps, rs = None) :
  seq = seqs[iSeq]

  stats = True
  ## cclust.counts does this counting faster
  ##
  ##   cans = [0]*len(seqs)
  ##   for i in range(len(seq)-11) :
  ##     p = seq[i:i+11]
  ##     ms = matches.get(p)
  ##     if ms :
  ##       for m in ms:
  ##         cans[m] += 1
  ##   cans[iSeq] = 0

  cans = [0]*len(seqs)
  cclust.counts(seq, matches, cans)
  cans[iSeq] = 0

  lseq = len(seq)
  # min is expensive (seriously! the volume of calls is staggering)
  # len too, that is why we pre-compute them.
  scans = [(-x/(lseq if lseq <= lseqs[k] else lseqs[k]), k)
           for k,x in enumerate(cans) if x > 0]
  
  # heap instead of full sorting, we usually need only a very small fraction of
  # the elements. Should move to faster C code. heapq is pure python
  heapq.heapify(scans)
      
  if rs:
    pfails,pok = 0,0
  tries = 0

  # Sequences we do not need to consider
  alreadyMatched = set()
  
  fails = 0
  matched = []
  distances = []

  # should we make the q sub-order on length (same nmatches, bring shorter one first)
  # Or another way to handle multiples from longer sequences?
  
  for cn,k in hiter(scans):
    if k in alreadyMatched :
      continue
    
    tries += 1

    djk = fdis(iSeq,k)
    if djk <= th :
      matched.append(k)
      distances.append(djk)

      fails = 0

      g = elem2grps.get(k)
      if g :
        for x in g :
          alreadyMatched.add(x)
          
      if rs:
        pok = cans[k]/min(lseqs[k],lseq)
        if pok < rs[0] :
          rs[0] = pok
        rs[1] += 1
        
    else :
      if rs and fails == 0 :
        pfails = cans[k]/min(lseqs[k],lseq)
      # rs[0] minimal accepted rs[1] n-accepted rs[2] total fails, rs[3] n-fails
      fails += 1
      
      # One mismatch can cause such a drop. dont give up before that.
      # heigher scores are possible due to multiple matchings.
      if cans[k] <= lseq-21 :
        if fails > _failsHardTH :
          break
        if fails > 1 and rs and rs[1] > 500 and rs[3] > 500 :
          p = cans[k]/min(lseqs[k],lseq)
          if p < .99*rs[0] and p < (rs[2]/rs[3] + rs[0])/2 :
            break
  if rs :
    rs[2] += pfails
    rs[3] += 1

  return matched,distances, (fails,tries,len(scans))


def deClutter(seqs, th, correction, failsTH = 20,
              matchScores = defaultMatchScores, verbose = None) :
  if verbose:
    tmain = time.clock()
    print >> verbose, "building lookup table ..., ",

  _failsHardTH = failsTH
    
  matches = _buildLookup(seqs, asSet = False)
  N = len(seqs)

  if verbose:
    print >> verbose, "done."
    print >> verbose, "n n-matched #singles #pairs #groups #matched-now #new-matched"

  # Duplicate code for speed 
  if correction :
    fdis = lambda i,j : calign.globalAlign(seqs[i], seqs[j], scores = matchScores,
                                           report = calign.JCcorrection)
  else :
    fdis = lambda i,j : calign.globalAlign(seqs[i], seqs[j],scores = matchScores,
                                           report = calign.DIVERGENCE)
    
  doStats = True
  if doStats :
    totTries = totFails = 0

  emptys = set()
  
  elem2grps = dict()
  lseqs = [len(s) for s in seqs]
  getp = lambda j : getPotentials(j, seqs, lseqs, th, fdis, matches, elem2grps, rs = rs)
  # 0,1 : min( seq-ident producing a match below th ), count
  # 2,3 - sum of p(first fail before terminating), count
  rs = [2,0,1,0]
  
  grps = []
  
  singles = set()
  pairs = list()
  pairsassingles = set()
  paired = set()
  for jSeq in range(N) :
    if jSeq in paired :
      continue

    pp = getp(jSeq)
    if doStats :
      iCans,sd,stat = pp
      totTries += stat[1]
      totFails += stat[0]
    else:
      iCans,sd = pp

    done = False
    if len(iCans) == 0 :
      singles.add(jSeq)
      done = True
      
    # Matching a previous single happens
    elif len(iCans) == 1 and iCans[0] not in paired \
             and iCans[0] not in singles :
      partner = iCans[0]
      pp = getp(partner)
      if doStats :
        jCans,sd,stat1 = pp
        totTries += stat1[1]
        totFails += stat1[0]
        stat = (stat[0],stat[1], stat[2] + stat1[2])
      else :
        jCans,sd = pp
      
      if len(jCans) == 1 and jCans[0] == jSeq :
        pairs.append([jSeq,partner,sd[0]])

        # probably not safe: avoid duplicates but those will not be examined
        # again for possible removal
        paired.add(jSeq)
        paired.add(partner)

        assert all([x not in singles and x not in pairsassingles for x in
                    [jSeq,partner]])
        
        pairsassingles.add(jSeq)
        pairsassingles.add(partner)
        
        done = True
      else :
        iCans = set(jCans + iCans)
        iCans.discard(jSeq)
        iCans = list(iCans)
        assert iCans is not None

    if verbose and not done:
      print >> verbose, jSeq, len(paired), len(singles), len(pairs), len(grps), len(iCans)+1, \
            sum([x not in paired for x in iCans])+1,\
            ("(%d %d  %d)" % (totTries,totFails,stat[2])) if doStats else ""
        
    if not done:
      g = [jSeq] + iCans
        
      paired.add(jSeq)
      sRemoved = []
      for x in iCans:
        if x in singles:
          sRemoved.append(x)
        singles.discard(x)
        paired.add(x)

      if verbose and len(sRemoved) > 0:
        print >> verbose, len(sRemoved),"removed from singles",\
              '('+",".join([str(x) for x in sRemoved])+')'

      for x in g :
        if x in pairsassingles :
          for k,(a,b,d) in enumerate(pairs):
            if a == x or b == x :
              del pairs[k]
              if a not in g :
                g.append(a)
              if b not in g :
                g.append(b)
              break
      if True:
        for x in g :
          if x in elem2grps:
            if not isinstance(g, set):
              g = set(g)
            g.update(elem2grps[x])
        if isinstance(g, set) :
          g = list(g)
        for x in g :
          elem2grps[x] = g
      grps.append(g)
        
  if verbose:
    print >> verbose, "merging", len(grps), "groups"
    
    if doStats:
      print >> verbose, totTries, "matching alignments,", totFails,"fails."
    
  mgrps = sorted(mergeGroupings(grps), key = lambda x : len(x), reverse=0)
  return list(singles),pairs,mgrps

# 10% [s,p g_1 g_2 ...  ] g_s small
# 3% [k g_k_s g_k_p  g_k_1 g_k_2 ... ] [l g_l_s g_l_p g_l_1 g_l_2 ...]
# 1% [k x_1_1 x_1_2

# 10% main 3%_1 %3_2 ....
# 3%_1 = 3% main 1%_1 1%_2 ...

def trans(ss, px, gx, grp) :
  "map partition elements back from indices using 'grp'"
  ss = [grp[x] for x in ss]
  px = [[grp[x],grp[y],h] for x,y,h in px]
  gx = [[grp[x] for x in g] for g in gx]
  return ss,px,gx
  
  
def breakGroup(grp, ths, dc, seqs, cnum, maxClade, result) :
  assert len(ths)
  
  assert len(grp) > maxClade > 0

  sq = [seqs[x] for x in grp]
  th = ths[0]
  #import pdb; pdb.set_trace()
  ss,px,gx = dc(sq, th)
  ss = sorted(ss)
  gx = sorted(gx, key = list.__len__)

  if len(ths) == 1 :
    ss,px,gx = trans(ss,px,gx, grp)
    c = [cnum + [0],ss,px,gx]
    result.append(c)
  else :
    i = len(gx)
    while i > 0 and len(gx[i-1]) > maxClade:
      i -= 1
    ss,px,gx = trans(ss,px,gx, grp)
    c = [cnum + [0],ss,px,gx[:i]]
    result.append(c)
    k1 = sum([len(x) for x in c[1:]])
    for k,g in enumerate(gx[i:]) :
      breakGroup(g, ths[1:], dc, seqs, cnum + [k1+k], maxClade, result)
    
def deClutterDown(seqs, ths, maxClade, correction, failsTH = 20,
                  matchScores = defaultMatchScores, verbose = None) :

  if len(seqs) <= maxClade :
    return [ [[],[],[],[range(len(seqs))] ] ]

  dc = lambda seqs, th : \
       deClutter(seqs, th, correction = correction, failsTH = failsTH,
                 matchScores = matchScores,  verbose = verbose)

  result = []
  breakGroup(range(len(seqs)), ths, dc, seqs, [], maxClade, result)
  return result

def thstr(th) :
  assert 0 < th < 1
  p = "%g" % th
  return 'p' + p[p.index('.')+1:]

def thfromstr(p) :
  assert p[0] == 'p' and p[1:].isdigit()
  return float('.' + p[1:])

def declutterToTrees(breakdown, ths) :
  trs = []
  cnt = 0

  addit = lambda t : trs.append((t, "cluster_%s_%s" %
                                 ('_'.join([namesuf,str(cnt)]) if namesuf else str(cnt) ,p)))

  for inds,singles,pairs,grps in breakdown:
    th = ths[len(inds)-1]
    p = thstr(th)
    
    namesuf = '_'.join([str(x) for x in inds[:-1]])
    cnt = inds[-1]
    for g in grps :
      d = th / (len(g)-1)
      h = d
      t = '(' * (len(g)-2)
      t = t + "(%d:%f,%d:%f)" % (g[0],d,g[1],d)
      for i in range(2, len(g)) :
        h += d
        t = t + (":%f,%d:%f)" % (d,g[i],h))
      addit(t)
      cnt += 1

    for i,j,d in pairs:
      d /= 2
      t = "(%s:%f,%s:%f)" % (i,d,j,d)
      addit(t)
      cnt += 1
      
    for s in singles:
      t = '(%d)' % s
      addit(t)
      cnt += 1
  return trs

def _cln2newick(n, tax) :
  if n.is_leaf() :
    return tax[n.id] if tax else str(n.id)
  
  ch = (n.get_left(), n.get_right())
  c = [_cln2newick(x, tax) for x in ch]
  return '(' + ",".join([(c[i] + ':' + ("%f" % (n.dist - ch[i].dist)))
                         for i in (0,1)]) + ')'

def upgma2tree(upgma, tax = None) :
  tr = scipy.cluster.hierarchy.to_tree(upgma)
  return _cln2newick(tr, tax)

def treeFromDists(ds, tax = None, weights = None, asString = False) :
  #up = scipy.cluster.hierarchy.average([x/2. for x in ds])
  #up = calign.upgma([x/2. for x in ds], weights = weights)
  up = calign.upgma(ds, weights = weights)
  if weights and any([x>1 for x in weights]) :
    lw = len(weights)
    #wt = lambda i : weights[i] if i < lw else nup[i-lw][3]
    wt = lambda i : 1 if i < lw else nup[i-lw][3]
    nup = []
    for i,j,d,w in up:
      w = wt(i) + wt(j)
      nup.append([i,j,d/2,w])
    up = nup
  else :
    up = tuple((a,b,c/2,d) for a,b,c,d in up)  # check that equiv
  
  tr = upgma2tree(up, tax)
  return tr if asString else parseNewick(tr)


def treeFromSeqs(seqs, tax = None, matchScores = None, correction = True,
                 weights = None, asString = False, verbose = None) :
  assert tax is None or len(seqs) == len(tax)
  
  if verbose :
    m = sum([len(x) for x in seqs])/len(seqs) # and??
    secs = (nPairs(seqs) * m**2)/(5/(20*1e-9))
    print >> verbose, len(seqs),tohms(secs),time.strftime("%T"),
    verbose.flush()
    tnow = time.clock()

  sseqs,order = zip(*sorted(zip(seqs,count())))

  ds = calign.distances(sseqs, align = True, scores = matchScores, reorder = order,
                        report = calign.JCcorrection if correction else calign.DIVERGENCE )

  tr = treeFromDists(ds, tax = tax, weights = weights, asString = asString)
  
  if verbose :
    print >> verbose, tohms(time.clock() - tnow), time.strftime("%T")
    
  return tr,ds

# th - half distance (i.e. height) 
def clusterFromTree(tr, th, caHelper = None) :
  if not caHelper :
    caHelper = CAhelper(tr)
  topNodes = []
  po = getPostOrder(tr)
  for nd in po:
    if nd.data.rh <= th :
      if nd.id == tr.root or nd.data.rh + nd.data.branchlength > th :
        topNodes.append(nd)
  return topNodes

# th distance, not node height  (i.e. **not** /2)
def cutForextAt(forest, th, helpersPool) :
  tops = []
  for it,tr in enumerate(forest) :
    if len(tr.get_terminals()) == 1 :
      tops.append([tr, tr.node(1)])
    else :
      ca = helpersPool(tr)
      nn = clusterFromTree(tr, th/2, ca)
      tops.extend([[tr,x] for x in nn])
  return tops


def assembleTree(trees, thFrom, thTo, getSeqForTaxon,
                 nMaxReps = 20, maxPerCons = 100,
                 lowDiversity = 0.02, thxx = 1.1, fantasyLand = .15,
                 verbose = None) :
  # cut trees at thFrom

  cahelpers = dict()
  cahelper = lambda t : cahelpers.get(t.name) or \
             (cahelpers.update([(t.name,CAhelper(t))])
              or cahelpers.get(t.name))

  if verbose:
    print >> verbose, "cutting",len(trees),"trees at %g" % thFrom
  
  pseudoTaxa  = cutForextAt(trees, thFrom, cahelper)
  nReps = len(pseudoTaxa)

  reps = [None]*nReps
  def getReps(k) :
    if not reps[k] :
      t,n = pseudoTaxa[k]
      nc = len(n.data.terms)
      if nc > 2:
        nc = min(max(int(math.log(nc,3)), 2), nMaxReps)
        r = random.sample(n.data.terms, nc)
      else :
        r = n.data.terms
      reps[k] = [getSeqForTaxon(x.data.taxon) for x in r]
    return reps[k]

  cons = [None]*nReps
  def getCons(k) :
    if not cons[k] :
      t,n = pseudoTaxa[k]
      nc = len(n.data.terms)
      if nc > maxPerCons :
        i = random.sample(n.data.terms, maxPerCons)
      else :
        i = n.data.terms
      sq = [getSeqForTaxon(x.data.taxon) for x in i]
      s, r = align.mpc(sq, nRefines=0)
      del r
      cons[k] = s
    return cons[k]

  mhs = []
  for t,n in pseudoTaxa:
    cahelper(t) # populate rh
    mhs.append(n.data.rh)

  # if both low diversity - use representative. If not valid or close to cluster height, do the
  # means thing. If not low diversity, use log representatives
  # low less then 4%??
  ## lowDiversity = 0.02
  ## thxx = 1.1
  ## fantasyLand = .15
  
  # counts how many alignments done (for display)
  global acnt
  acnt = 0
  
  def getDist(i,j) :
    mi,mj = mhs[i],mhs[j]
    anyCons = False
    if mi < lowDiversity :
      ri = [getCons(i)]
      anyCons = True
    else :
      ri = getReps(i)
    if mj < lowDiversity :
      rj = [getCons(j)]
      anyCons = True
    else :
      rj = getReps(j)

    nhs = len(ri)*len(rj)
    if nhs == 1 :
      h = calign.globalAlign(ri[0], rj[0], scores = defaultMatchScores,
                             report = calign.JCcorrection)
    else :
      ap = calign.allpairs(ri, rj, align=True, scores = defaultMatchScores,
                           report = calign.JCcorrection)

      h = sum([sum(x) for x in ap])/nhs
      
    global acnt
    acnt += nhs
    
    lowLim = 2*max(mi,mj)
    
    if anyCons and (h < lowLim or (h < fantasyLand and h < lowLim*thxx)) :
      xri = getReps(i) if len(ri) == 1 else ri
      xrj = getReps(j) if len(rj) == 1 else rj

      if ri != xri or rj != xrj :
        ap1 = calign.allpairs(xri, xrj, align=True, scores = defaultMatchScores,
                              report = calign.JCcorrection)
        h1 = sum([sum(x) for x in ap1])

        xnhs = (len(xri)*len(xrj))
        acnt += xnhs
      
        h = (h * nhs + h1)/(nhs + xnhs)
    return max(h, lowLim)

  if verbose :
    print >> verbose, "assembling",nReps,"sub-trees into one tree",time.strftime("%T")
    print "n-sub-tree #pair-only-alignments #alignments time"
    verbose.flush()
    tnow = time.clock()

  # use array??
  ds = [-1]*nPairs(nReps)

  pos = 0
  for i in range(nReps-1) :
    for j in range(i+1, nReps) :
      ds[pos] = getDist(i,j)
      pos += 1

    if verbose :
      dn = sum(range(nReps-1, nReps-i-2,-1))
      print >> verbose, i, dn, "%3.2g%%" % ((100.*dn)/len(ds)), acnt, time.strftime("%T")

  if verbose :
    print >> verbose, tohms(time.clock() - tnow), time.strftime("%T")

  # Using correct weights can throw off the height guarantee, or not?
  wt = [len(n.data.terms) for t,n in pseudoTaxa]

  tnew = treeFromDists(ds, tax = [str(x) for x in range(nReps)], weights = wt)

  saveDistancesMatrix('/home/joseph/research/Projects/AWCCoRE/data/iii/CO1/work/dists-inter1', ds, [str(x) for x in range(nReps)], compress=9)
  del ds

  for n in getPostOrder(tnew) :
    if not n.succ :
      t,nd = pseudoTaxa[int(n.data.taxon)]
      if len(nd.data.terms) == 1 :
        n.data.taxon = nd.data.taxon
        n.data.rtree = "%s:%f" % (n.data.taxon, n.data.branchlength)
      else :
        # Insure heights are there
        cahelper(t)
        s = t.toNewick(nd.id)
        d = n.data.branchlength - nd.data.rh
        if (d < -1e-10) :
          print "***** ERROR", d
        n.data.rtree = "%s:%f" % (s, max(d,0.0))
    else :
      ch = [tnew.node(x).data.rtree for x in n.succ]
      n.data.rtree = "(%s,%s)" % (ch[0],ch[1])
      if n.id != tnew.root :
        n.data.rtree = n.data.rtree + (":%f" % n.data.branchlength)

  trec = tnew.node(tnew.root).data.rtree
  trec = parseNewick(trec)
  return trec
