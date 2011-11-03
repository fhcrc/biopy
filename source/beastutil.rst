===============
BEAST utilities
===============

-------------------
starbeast_demog_log
-------------------

    A specialized utility for \*BEAST users. The population size numbers as
    reported in the BEAST log can be misleading, since they may refer to
    diffrent populations at diffrent stages of the MCMC chain. This utility
    generates a new trace file, viewable in tracer, tracing population size(s)
    for each population encountered during the run. In addition, a trace is
    created for the species tree topology using a simply heuristic.

------------
summary_tree
------------

    Create a summary tree from posterior trees. BEAST tree annotator utility
    first picks up a topology, then independently assigns heights to each
    internal node. Sometimes this can be suboptimal - in the worst case the tree
    contains negative branch lengths. summary_tree takes a more global approach
    and looks for a tree which minimizes the overall distance to the whole set
    of posterior trees, using a tree distance measure.

    For more details see *Posterior summary of trees* (page 68 of my
    :ref:`thesis<myThesis>`).

----------------------------
starbeast_posterior_popsizes
----------------------------

    Annotate a \*BEAST summary tree (such as the one from summary_tree
    above) with posterior estimates of population sizes.

    On UNIX it is easy to get the tree directly on the command line, e.g

| **$ starbeast_posterior_popsizes `summary_tree trees.nexus` trees.nexus**

--------------
sptree_plot
--------------

    Generate a figure showing all posterior species trees from a
    \*BEAST analysis. The figure can be used to visually examine the
    uncertainty of both divergence times and population sizes.

    Inspired by `DensiTree
    <www.cs.auckland.ac.nz/~remco/DensiTree/DensiTree.html>`_. Some of the ideas
    developed here became part of version 2.

.. image:: run02tmc20.png

  
-----------------------------------
multispecies_coalescent_estimate
-----------------------------------

    Heuristically estimate the birth rate and effective population
    size from a multispecies data in a \*BEAST XML file.

    See `Average Sequence dissimilarity under simple multi-species coalescent <http://arxiv.org/abs/1104.0727>`_.
