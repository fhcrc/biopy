## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.

#
#
# Commenting your code is like cleaning your bathroom -- you never want to do
# it, but it really does create a more pleasant experience for you and your
# guests.  Ryan Campbell
#            (http://particletree.com/features/successful-strategies-for
#                                                         -commenting-your-code
#

import math, random, copy, sys
from math import log, ceil

from numpy import array, exp, dot, multiply, cumsum, prod, alltrue, matrix, identity
from numpy.linalg import eig, inv
from numpy import float as npfloat

from cchelp import seqevolve

def _siteProbs(n) :
  p = [0,0,0,0]
  p[n] = 1
  return p

def pick(p, k) :
  """ Randomly pick a integer in [0-3] according to cumulative probabilities
  in p[k,0..3]."""

  # Profiling suggests pick is used a lot during simulation. reason to
  # optimize.
  # slower
  # return [x < r for x in p].index(False)
  
  r = random.random()

  if r < p[k,1] :
    if r < p[k,0] : return 0
    return 1
  if r < p[k,2] : return 2
  return 3

class SubstitutionModel(object) :
  """ Base class for any General Time Reversible (GTR) substitution model """

  NUCS = "AGCT?"
  
  def __init__(self, m, pi, mu, name = "Substitution Model") :
    """ m    - six transition rate parameters.
        pi   - four stationary probabilities at time 0.
        mu   - mutation rate (scales time to substitutions).
        name - string for debug/development message."""
    
    self.name = name
    self.m = array(m)

    self.pi0 = array(pi)
    assert all(self.pi0 > 0) and abs(sum(self.pi0) - 1.0) < 1e-15
    
    self.setMu(mu)

  def __str__(self) :
    return self.name

  @staticmethod
  def toNucCode(seq) :
    """ Internal codes for char string."""
    return "".join([SubstitutionModel.NUCS[n] for n in seq])

  @staticmethod
  def toSeq(s) :
    """ Char string to internal codes."""
    return [SubstitutionModel.NUCS.index(x) for x in s]
  
  def setMu(self, mu) :
    """ Change mutation rate. """
    
    self.mu = float(mu)                                  ;assert self.mu > 0.0
    self.qInitial = self.qFromGTR()
    self.initialTotalSubRate = self.totalSubRate()
    self.qNorm = self.qInitial / self.initialTotalSubRate
    
  def totalSubRate(self) :
    """ Total substitution rate.
    Typically mu is normalized so that total rate is 1. We do the same for
    t == 0"""
    pi = self.pi0
    
    q = self.qInitial

    r = -sum([pi[i] * q[i,i] for i in range(4)])
    return r
  
  def qFromGTR(self) :
    """ Build explict Q matrix from GTR parameters at time 0"""
    
    p = self.m
    pi = self.pi0
    
    m = array( [ [0.0,        p[0]*pi[1], p[1]*pi[2], p[2]*pi[3]] ,
                 [p[0]*pi[0], 0.0,        p[3]*pi[2], p[4]*pi[3]] ,
                 [p[1]*pi[0], p[3]*pi[1], 0.0,        p[5]*pi[3]] ,
                 [p[2]*pi[0], p[4]*pi[1], p[5]*pi[2], 0.0] ] )
    for i in range(4) :
      m[i,i] = -sum(m[i,])

    q = m * self.mu
    return q

  def q(self) :
    """ Q (rate) matrix.
    Constant Q normalized to substitutions and total rate is 1. """
    
    return self.qNorm * self.mu
  
  def setTransitionSpeedup(self, tree) :
    for n in tree.all_ids() :
      node = tree.node(n)
      if n != tree.root :
        node.data.pMatrix = self.pExact(node.data.branchlength)

  def clearTransitionSpeedup(self, tree) :
    for n in tree.all_ids() :
      node = tree.node(n)
      if n != tree.root :
        del node.data.pMatrix
    
  def simulateOverInterval(self, branch, nucs, deltaT) :
    """ Evolve sites over [t,t+branch], using steps of size deltaT."""
    
    nTimes =  ceil(branch / deltaT)
    deltaT = branch / nTimes

    qmdt = self.q() * deltaT
    p = cumsum(identity(4) + qmdt, 1)
    for i in range(int(nTimes)) :
      nucs = [pick(p,n) for n in nucs]

    return nucs

  def pSimulated(self, branchLen, nSites = 1000, deltaT = 0.01, t = None) :
    """ P matrix over branch by counting substitution occuring as a result
    of simulation in small steps over the interval."""
    
    if t is None: t = -branchLen
    
    p = identity(4, npfloat)

    for i in range(4) :
      s = self.simulateOverInterval(branchLen, (i,)*nSites, deltaT)
      p[i,] = [sum([x == k for x in s]) / float(len(s)) for k in range(4)]

    return matrix(p)

  def p(self, branch) :
    """ P (transition probability) matrix over branch. """
    return self.pExact(branch)

  def evolve(self, seq, pMatrix) :
    """ Evolve sequence C{seq} according to transition probabilities in
    C{pMatrix}. """
    return seqevolve(pMatrix, list(seq))

    #p = cumsum(pMatrix, 1)
    #return [pick(p,n) for n in seq]
    
  def populateSubtreeSeq(self, tree, nid, seq) :
    """ Populate sub-tree with sequences evolved from ancestor C{nid} with data
    C{seq}. """
    
    node = tree.node(nid)
    # If leaf, simply take the sequences
    if node.data.taxon :
      node.data.seq = seq
    else :
      # Do each side recursively
      for n in node.succ:
        subSeq = self.evolve(seq, tree.node(n).data.pMatrix)
        self.populateSubtreeSeq(tree, n, subSeq)
      
  def populateTreeSeqBySimulation(self, tree, seqLen) :
    """ Populate tree tips, each with a sequences of length 'seqLen', evolved
    from a common ancestor drawn from the stationary probabilities."""
    
    self.setTransitionSpeedup(tree)

    # Draw root from stationary distribution
    p = array([cumsum(self.pi())])
    seq = [pick(p,0) for i in range(seqLen)]
    
    self.populateSubtreeSeq(tree, tree.root, seq)
    
    self.clearTransitionSpeedup(tree)

  def condProbsSubTree(self, tree, nid) :
    node = tree.node(nid)
    data = node.data
    # If leaf, simply take the sequences
    if data.taxon :
      subsCond = [_siteProbs(site) for site in data.seq]
    else :
      s = [self.condProbsSubTree(tree, n) for n in node.succ]
      subsCond = [[x*y for x,y in zip(l,r)] for l,r in zip(*s)]

    if nid != tree.root:
      pMat = data.pMatrix
      subsCond = [dot(pMat, cp) for cp in subsCond]
    return subsCond

  def logLike(self, tree) :
    self.setTransitionSpeedup(tree)
    subsCond = self.condProbsSubTree(tree, tree.root)
    self.clearTransitionSpeedup(tree)
    pi = self.pi0
    l = sum([log(dot(pi, x)) for x in subsCond])
    return l

#
# jc = msubmodels.JCSubstitutionModel(mu = 0.005)
# t = md2.INexus.Tree('((a:1,b:1):1,c:2)')
# jc.populateTreeSeqBySimulation(t, 1)
# r = []
# for x in itertools.product(*itertools.tee(range(4),3)):
#   for n,z in zip(t.get_terminals(), x) :
#     t.node(n).data.seq[0] = z
#   r.append(jc.logLike(t))
# print sum([exp(x) for x in r])
    
class StationaryGTR(SubstitutionModel) :
  def __init__(self, m, pi, mu, name = None) :
    if name is None:
      name = "GTR mu %g rates %s pi %s" % (mu, ",".join(["%g" % x for x in m]),
                                           ",".join(["%g" % x for x in pi]))
    SubstitutionModel.__init__(self, m = m, pi = pi, mu = mu, name = name)
    
    if self.pExact.im_func == StationaryGTR.pExact.im_func :
      self.w,self.v = eig(self.q())
      self.iv = inv(self.v)

  def pi(self) :
    return self.pi0
  
  def QonInterval(self, branch) :
    """ Element-wise integral of Q over [t,t+branch] """
    assert branch >= 0
    
    return self.q() * branch

  def pExact(self, branch) :
    """Exact solution for P."""
    assert branch >= 0

    return dot(multiply(self.v, exp(self.w*branch)), self.iv)

class JCSubstitutionModel(StationaryGTR) :
  def __init__(self, mu = 1, name = None) :
    if name is None :
      name = "JC69 mu %g" % mu
    SubstitutionModel.__init__(self, m = (1,)*6, pi = (1/4.0,)*4, mu = mu, \
                               name = name)

  def pExact(self, branch) :
    """Exact analytic solution for P."""

    nd = 0.25 - 0.25 * math.exp((-4.0/3.0) * self.mu * branch)
    d = 1 - 3*nd
    return array([[d,nd,nd,nd],[nd,d,nd,nd],[nd,nd,d,nd],[nd,nd,nd,d]])
  

class Kimura2P(object) :
  def ab(self, r) :
    # Get a,b from solving
    # r = a/2b , a + 2b = 1
    #
    # written in such a way which works for array r as well. (no big deal,
    # very natural in python).
    
    b = 0.5 / (r + 1.0)
    a = 1 - 2*b
    return (a, b)
  
class Kimura2PSubstitutionModel(StationaryGTR, Kimura2P) :
  def __init__(self, mu = 1, kappa = 1, name = None) :
    self.kappa = kappa
    r = float(kappa)/2
    a,b = self.ab(r)
    if not name :
      name = "Kimura 2P(R=%g)" % r
      
    StationaryGTR.__init__(self, m = (a,b,b,b,b,a), pi = (1/4.0,)*4, \
                            mu = mu, name = name)
  
  def pExact(self, branch) :
    """Exact analytic solution for P."""

    b = branch * self.mu
    alpha,beta = self.m[0:2]
    x1 = math.exp(-4 * beta * b)
    k = 0.25 * (1 - x1)
    m = 0.5 * ( x1 - math.exp(-2 * (beta+alpha) * b) )

    K = matrix(((-3,1,1,1),(1,-3,1,1),(1,1,-3,1),(1,1,1,-3)), npfloat)
    M = matrix(((-1,1,0,0),(1,-1,0,0),(0,0,-1,1),(0,0,1,-1)), npfloat)

    return array(identity(4) + k * K + m * M)

class HKYSubstitutionModel(StationaryGTR, Kimura2P) :
  def __init__(self, mu = 1, kappa = 1, pi = [.25,]*4, name = None) :

    if name is None :
      name = "HKY pi (%g %g %g %g) kappa %g mu %g" \
             % (pi[0], pi[1], pi[2], pi[3], kappa, mu)

    self.kappa = float(kappa)
    
    a,b = self.ab(kappa/2.0)
    
    StationaryGTR.__init__(self, m = (a,b,b,b,b,a), pi = pi,
                           mu = mu, name = name)
    
