====================================
Proximate UPGMA Phylogeny (PUP) tree
====================================

Next Generation Sequencing produces a prodigious amount of relatively short
sequences. Those sequences are usually too numerous and too diverse to align,
the prerequisite for using your favourite tree-building method. Even vanilla
pairwise methods such as NJ and UPGMA struggle because computing :math:`O(n^2)`
alignments and storing the results is prohibitively expensive in both time and
space (and NJ is :math:`O(n^3)`). The *otool* utility can build an ultrametric,
UPGMA-like tree for a large (200,000 or even more) number of sequences. While
feasible, building those trees can take a long time (days) depending on size and
sequence diversity.

You can build a tree from ``myfile.fasta`` with the following command:

| $ ``otool --progress --max-single-tree 4000 --derep --declutter 0.2,0.1,0.03,0.01 --forest --stitch 1 myfile.fasta``

This will output a single tree to the NEXUS file ``myfile.derep.trees``, and
create a bunch of other files, so you want to run the above in an empty
directory. So far I was not able to figure out how to find good parameters
values from the data, so you may need to experiment with the decluttering
parameters to lower the total running time. See the comments for the
``--declutter`` option.

General information:
--------------------

The tree is UPGMA-like, based on genetic pairwise distances. The distance is
computed by aligning two sequences, computing the sequence identity from the
alignment and applying the Jukes-Cantor correction.

Building the tree, like any proper military activity, is split into 3 parts:
decluttering, building a forest and stitching. In addition there is one
pre-option (de-replicating) and two post options (shave and cluster).

You can build the tree in a single call or by stages. There is little difference
since *otool* runs the stages in a sequence anyway, but this setup provides more
flexibility. For example, you can omit stitching and generate a forest instead,
or select between several decluttering runs with different parameters.

There are no options for specifying output file names. Several files are
generated in the various stages, and their names are based on the original input
file. That way you can specify a single input file name and the program can
locate other files, if needed, without the need to specify them explicitly on
the command line. I find this convenient and I hope it will reduce your
confusion too.


Basic usage:
------------

::

 In most cases you want to start with de-replicating the data, that is, creating one "master" sequence for groups of
 identical sequences. Sequence names will be augmented by a suffix of the form ';size=NNN;' when 'NNN' is the duplicity
 number (same format as USEARCH).  You can skip this step but the program will take longer to complete and the end
 result will be the same. Skip de-replicating is you know the data contains no duplicates or if you wish to keep
 duplicate distinct (say sequences are tagged with a species name and identical sequences from different species are
 allowed). The size is used throughout the building process and you can think of the sequence as representing a
 "collapsed" sub-tree with zero branches. otool will create the output file 'myfile.derep.fasta'.
 
| $ ``otool --derep myfile.fasta``

::

 The first stage *declutters* the data, partitioning the sequences into e-islands, where sequence from two different
 islands are at least e apart. This is a hierarchical process, i.e. '--declutter 0.2,0.1,0.03' partitions the data into
 0.2-islands, where each 0.2-island may be partitioned further into 0.1-islands, which in turn may be partitioned into
 0.03-islands. The islands are split based on size: islands larger than 4000 sequences are split (in the example below),
 but smaller islands are not.
 The total running time, and to some extent the "quality" of the tree depends on the maximum tree size and the declutter
 levels. Generally you want larger islands: The maximum size on my 4GB laptop is around 11,000, but you probably want to
 keep it between 3000 and 7000, since larger islands tend to have a root age greater than half the level and would be
 split into a forest before being stitched into a larger island. As to the levels, you probably want the least number of
 levels up to 0.1, but above that additional levels are less of a worry since the phylogenetic signal is much weaker at
 those time scales and sequence lengths. With '-p', *otool* will print a very rough estimate of the time required to
 complete building the tree. While the time may be inaccurate, it is very useful for comparing different settings.
 Decluttering outputs the island in a NEXUS file ``myfile.derep.decl_p2_p1_p03.trees`` (where the _pXX stands for the
 levels). The branch lengths of the trees in this file arbitrary, and the island affiliation is coded in the tree name.

| $ ``otool [--max-single-tree=NN] --declutter TH1,TH2,... [--max-single-tree=NN] myfile.derep.fasta``

::

 Build a UPGMA tree for each of islands in the file. Creates a NEXUS file ``myfile.derep.pXX.trees``, where XX is the
 forest *level*, that is, the smallest declutter size (``myfile.derep.p03.trees`` for the example above).
 
| $ ``otool --forest myfile.derep.decl_pXX_pXX_...trees``

::

 Merge trees built from islands with level below TH to one tree, which represents an TH-island. For example,
 ``--stitch 0.1`` will merge groups of 0.03-islands back into the 0.1-island that was split during declutting.
 ``--stitch 1`` would build a single tree in the NEXUS file  ``myfile.derep.trees``. 
 
| $ ``otool --stitch TH myfile.derep.pXX_pXX_...trees``


::

  Cluster the sequences into putative OTUs. Clusters are formed from maximal monophyletic clades whose with root age
  less than TH/2. That is, the distance between any pair of sequences in a cluster is less than TH and the distance
  between sequence from different clusters is greater than TH -- according to the distances defined by the tree. Of
  course, any specific pairwise distance might violate this condition, since the distances induced by the tree are an
  ultrametric reconciliation of the non-ultrametric pairwise relationships. The text output file 'myfile.derep.clusters'
  contains one line per cluster, and each line contains the sequence names in the cluster separated by tabs.
  
| $ ``otool --cluster TH myfile.derep.trees``

::

  Partition the sequences into putative OTUs at the ``TH`` level (same logic as for cluster above), and generate a FASTA
  file with one consensus sequence for each OTU. The names in the 'myfile.derep.shave_pTH.fasta' file include the names
  of two tips, where the common ancestor of those two is the OTU root.
  
| $ ``otool --shave TH myfile.derep.trees``


Generic options:
----------------

::

 Help on usage and commands.
 
| $ ``otool --help``

::

 Output progress messages for your entertainment during the build process. Repeating the option will generate lots of
 output useful for developers.
 
| $ ``otool --progress ...``
| $ ``otool --progress --progress ...``


Advanced options:
-----------------

::

  The search for matching sequences during declutter are terminated after N failures. Increasing the limit increases the
  search sensitivity and increasing running time cost. You probably want to keep this under 40.
		    
| ``--terminate-search-count N[=20]``

::

  Control saving of pairwise distances during forest building. Those saved distances can save time if the program is
  unexpectedly terminated during the second stage. The default is 'compressed-fast', equivalent to compressing with
  'gzip -1 ...'. Set to 'no' to disable saving. The distance files, '*.dists.gz', can be safely deleted at the end of
  the run.
  
| ``--save-distances no|plain|compressed-fast|compressed-normal|compressed-best``

::

  By default otool will not overwrite output files. Use this option to allow overwriting existing output files.

| ``--overwrite``

::

 Set the aligning parameters. Sequences are globally aligned to maximize the score. Each match contributes an M, each
 mismatch an X, the first gap a G and each subsequent gap E. This is all fairly standard, and F=1 sets the free-end-gaps
 option, that is, a series of gaps at the alignment end (or the beginning) have no cost (zero score). This option is
 especially relevant for NGS data, since reads of the same sequence can have various lengths, and we don't want to pay
 the cost of mismatches or gaps when aligning those reads. The default is 10,-5,-6,-6,1.
 
| ``--match-scores M,X,G,E,F``

::

 Use sequence identity, that is, do *not* apply the Jukes-Cantor correction.
  
| ``--use-sequence-identity``

::

 For the experts. Leave those at the default settings.

| ``--stitch-parameters``

..  LocalWords:  ultrametric decluttering declutter declutting fasta clades

.. --max-single-tree 3000 --declutter 0.25,0.2,0.1,0.03,0.01 
.. Rough build time 265:54:42 - 266:18:21 (333090381-333584308),

.. --max-single-tree 3000 --declutter 0.2,0.15,0.1,0.06,0.03,0.01 
.. Rough build time 110:36:51 - 118:15:14 (138559395-148129146)

.. --max-single-tree 3000 --declutter 0.2,0.15,0.125,0.1,0.08,0.06,0.03,0.01
.. Rough build time 67:52:09 - 75:30:32 (85015537-94585288),

.. --max-single-tree 4000 --declutter 0.2,0.15,0.125,0.1,0.08,0.06,0.03,0.01
.. Rough build time 67:52:09 - 75:30:32 (85015537-94585288),

.. --max-single-tree 3000 --declutter 0.18,0.15,0.125,0.1,0.08,0.06,0.04,0.02,0.01
.. Rough build time 66:01:04 - 89:29:42 (82696292-112104854),
..  LocalWords:  monophyletic
