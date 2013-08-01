===================
Sequence generation
===================

------------------
simulate_sequences
------------------

    Simulation sequences for tips of a gene tree using a standard GTR
    substitution model. Jukes-Cantor and HKY are provided for convenience.

    The substitution model is specified by the ``-m`` option. The parameter is
    a comma separated list (no spaces), where the first element is the model
    name and the rest are the model parameters. The Jukes-Cantor takes only a
    mutation rate. The HKY takes a mutation rate, kappa and up to 4 stationary
    probabilities. The GTR takes a mutation rate, up to 7 rates and 4 stationary
    probabilities. 

      * JC,mu
      * HKY,mu,kappa,A,pi_A,...
      * GTR,mu,AG,r_AG,CG,r_CG,...,A,pi_A
      
    The rates and probabilities can come at any order. An unspecified or missing
    parameter is assigned a default value. The default for the mutation rate,
    kappa or rates is 1. For the stationary probabilities, each of the
    unassigned probabilities is set to the same value so that all 4 add to one.
    For example, setting only pi_A to 0.7 will set C,G, and T to .1 = (1-0.7)/3.

Usage:
^^^^^^

::

  Use default Jukes-Cantor model with mutation rate 1. 5bp sequence
  length. Output an annotated tree.
  
| **$ simulate_sequences -n 5 -a "(a:1,b:1)"**
| (a[&seq=TGCCT]:1.0,b[&seq=AGCGG]:1.0)

::

  Create a NEXUS file with the alignment.

| **$ simulate_sequences.py -n 5 "(a:1,b:1)"**
| **$ cat alignment0.nex**

::

  #NEXUS
  Begin data;
   	Dimensions ntax=2 nchar=5;
   	Format datatype=dna gap=-;
   	Matrix
  a 	TGCCT
  b 	AGCGG
   	;
  End;

::

  Jukes-Cantor with mutation rate 0.05
  
| **$ simulate_sequences.py -m JC,0.05 -n 5 -a "(a:1,b:1)"**
| (a[&seq=ATACT]:1.0,b[&seq=ATACT]:1.0)

::

  HKY with the default mutation rate (1), Kappa 2 and equal stationary probabilities.

| **$ simulate_sequences.py -m HKY,,2 -n 10 -a "(a:1,b:1)"**
| (a[&seq=CCCTTACCGC]:1.0,b[&seq=GTGCTTTCAC]:1.0)

::

  HKY with mutation rate 0.1, Kappa 2 and stationary probabilities of p(A) = .4
  and P(C,G,T) = 0.2

| **$ simulate_sequences.py -m HKY,.1,2,A,.4 -n 10 -a "(a:1,b:1)"**
| (a[&seq=CCCTTACCGC]:1.0,b[&seq=GTGCTTTCAC]:1.0)


::

  GTR with mutation rate 0.2. Stationary probabilities are of p(A) = .4, P(T) =
  0.4 and P(C,G) = 0.1. The A<-->G rate in the Q matrix is 2, the C<-->T rate is
  3, all other rates are 1 (the default).

| **$ simulate_sequences.py -m GTR,.2,AG,2,CT,3,A,.4,T,.4 -n 10 -a "(a:1,b:1)"**
| (a[&seq=TCTCTTTACT]:1.0,b[&seq=TCTCTTTACT]:1.0)
