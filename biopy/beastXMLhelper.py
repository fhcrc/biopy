## This file is part of biopy.
## Copyright (C) 2010 Joseph Heled
## Author: Joseph Heled <jheled@gmail.com>
## See the files gpl.txt and lgpl.txt for copying conditions.
#

""" Miscellaneous helpers for reading BEAST xml files.
"""

from __future__ import division
from lxml import etree

__all__ = ["readBeastFile"]

def getAlName(treeElementTag, gtreeName, root) :
  pats = None
  for tl in root.findall('treeLikelihood') :
    t = tl.find(treeElementTag)
    if t is not None and t.attrib['idref'] == gtreeName :
      pats = tl.find('patterns').attrib['idref']
      break
  alName = None
  for p in root.findall('patterns') :
    if p.attrib['id'] == pats :
      alName = p.find('alignment').attrib['idref']
      break
  return alName

def readBeastFile(path, what) :
  """ Read portions of BEAST XML file for use in scripts.

  The keys in the dictionary C{what} indicate which portions to read.  Currently
  supported are 'DNA' for the alignments, 'mcmc' for some basic chain
  information, 'logging' for output file names and 'species' for
  individuals/species/genes mappings in a *BEAST file.

  @param path: BEAST file name
  @type path: str
  @param what: Fill information 
  @type path: dict
"""
  doc = etree.parse(path)
  root = doc.getroot()

  data = dict()
  if "DNA" in what :
    for a in root.findall("alignment") :
      aName = a.get('id')
      seqs = dict()
      for s in a.findall("sequence") :
        taxon = s.find('taxon')

        n = s[0].tail.strip()
        assert len(n)
        assert all([c in "AGCT-MRWSYKVHDBXN?" for c in n.upper()])
        seqs[taxon.get('idref')] = n

      l = len(seqs[iter(seqs).next()])
      for s in seqs :
        assert len(seqs[s]) == l

      data[aName] = seqs

    what['DNA'] = data

  if "mcmc" in what :
    m = root.find("mcmc")
    cl = int(m.get('chainLength'))
    for l in m.findall('log') :
      if l.get('id') == 'fileLog' :
        le = int(l.get('logEvery'))
        break
  
    what['mcmc'] = {'chainLength' : cl, 'logEvery' : le}

  if "species" in what :
    xspecies = root.find('species')
    if xspecies is not None :
      species = dict()

      for sp in xspecies.findall('sp') :
        spName = sp.attrib['id']
        species[spName] = [tx.attrib['idref'] for tx in sp.findall('taxon')]

      genes = dict()
      gts = xspecies.find('geneTrees')

      for gt in gts.findall('gtree') :
        ploidy = gt.attrib.get('ploidy')
        ploidy = float(ploidy) if ploidy else None
        gtreeName = gt[0].attrib['idref']
        treeElementName = gt[0].tag
        alName = getAlName(treeElementName, gtreeName, root) 
        if alName:    
          genes[alName] = {'ploidy' : ploidy, 'tree' : gtreeName}
      for gt in gts.findall('treeModel') :
        gtreeName = gt.attrib['idref']
        alName = getAlName(gt.tag, gtreeName, root) 
        if alName:    
          genes[alName] = {'ploidy' : 1, 'tree' : gtreeName}
      
      what['species'] = { 'species' : species, 'genes' : genes }

  if "logging" in what :
    mcmc = root.find('mcmc')

    treeLogs = dict()
    for lt in mcmc.findall('logTree') :
      t = lt[0]
      treeLogs[lt.attrib['fileName']] = [t.attrib['idref'], t.tag]

    logs = dict()
    for l in mcmc.findall('log'):
      if 'fileName' in l.attrib:
        d = [e.attrib.get('idref') or e.attrib.get('id') for e in l]
        logs[l.attrib['fileName']] = d
      
    what['logging']  = { 'logs' : logs, 'treelogs' : treeLogs }