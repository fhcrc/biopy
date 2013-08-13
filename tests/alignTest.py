from __future__ import division
from math import log,exp

from calign import globalAlign, createProfile, profileAlign, distances, DIVERGENCE, IDENTITY, JCcorrection

#scores = (10,-5,-6,None,False)
scores = (10,-5,-6,-6,False)
fescores = (10,-5,-6,-6,True)

s1 = """AAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTAGTGATTGCCATCTTGGCTTAAACTATATCCATCTACACCTGTGAACTGTTTGATTGAATCTCACGATTCAATTCTTTACAAACATTGTGTAATGAACGTCATTAGATCATAACAAAAAAACTTTAACTAACGGATCTCTTGGCTCTCGCATCGATGAAGAACGCAGCGGTCATAGCTGTTTCC"""

s2 = """AAGTCGTAACGAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTAAAGGTTTCGGGTACCCAGTGCCCAAACTTCAACCCTATGTTTTAACTATACATGTTTCTTTGCCGGTTTCGGCCGGCAGAAGTTTTCTCAAACTCATTATAAATGTGTCTTCTGAATCAAAACTAAATAAATTATAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCGGTCATAGCTGTTTCCTG"""

s3 = """CGCAGTGGGGTGGTGTGCTCGATGTGGCCGGTGCCGATCGCGTGACCGTTGGCCGTCAACGTCACGGTGCCGCCTTTGCCGAAATCCTTAGGGCCGCCGTCGTATTTGAAGTTGAATCCCAGTGTGACCCTGCCGGTCGGAACTGACTGGCTCGACGTCACCCGGTAGCGCGCGGCATCGATGAAGAACGCAGCGGTCATAGCTGTTTCCTG"""

def test00() :
  """
>>> globalAlign("AAA", "AAA", report = IDENTITY)
1.0
>>> globalAlign("AAA", "AAA", report = DIVERGENCE)
0.0
>>> 3*globalAlign("AAA", "AAC", report = IDENTITY, scores=scores)
2.0
>>> (globalAlign("AAA", "AAC", report = JCcorrection, scores=scores) - (-3/4. * log(5/9.))) < 1e-12
True
>>> globalAlign("AAGG", "AAT", report = IDENTITY, scores=scores)
0.5
>>> 3*globalAlign("AAGG", "AAT", report = IDENTITY, scores=fescores)
2.0
>>> globalAlign("AAGG", "TAA", scores=fescores)
((5, 0, 0, 1, 1), (3, 0, 0, 5, 5))
>>> 265.0*globalAlign(s1, s2, scores=scores, report=IDENTITY)
188.0
>>> 263*globalAlign(s1, s2, scores=fescores, report=IDENTITY)
188.0
>>> a = globalAlign(s1, s2, scores=scores) ; p = createProfile(a) ; b = profileAlign(s3,p) ; distances((b,) + a, align=0)
(0.508, 0.4921875, 0.28517110266159695)
"""
  pass

if __name__ == '__main__':
  import doctest
  doctest.testmod()

