.. _randDists:

===============================
Specifiying random distrbutions
===============================

A number by itself specifies a constant (the delta
distrbution). Otherwise, the argument should be a comma separated list
(no spaces), where the first character specifies the distribution.

Supportd distributions:
 * u,L,H : Uniform between L and H (Example ``u,2,4``).
 * e,R   : Exponential with rate R (mean 1/R) (Example ``e,2``).
 * l,M,S : Log-normal with mean M (in real space) and std S (in log space).
 * g,S,L : Gamma with shape S and scale L.
 * i,A,B : Inverse-gamma with shape(alpha) A and scale(beta) B.
 * n,M,S : Normal with mean M and std S.

