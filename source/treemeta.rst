
--------------
Tree meta data
--------------

Meta data is associated with tree nodes, and is placed inside matching square
brackets between the node and the branch length. The first character must be a
``&``, and each node may contain several comma separate items, where each item
has the form ``NAME=VALUE``. The value may be enclosed with a matching ``{}`` or
``""`` pair. This is especially handy when the value contains commas. For example,

::

  (a[&colors={red,blue}]:1,b[&color=green,location="7,8.1"]:1)

.. _MetapopSize:
  
Setting population sizes
^^^^^^^^^^^^^^^^^^^^^^^^

The ``dmv`` attribute specifies population sizes. For example,


::

  (a[&dmv=2]:1,b[&dmv=1]:1)[&dmv=3]


assigns a constant population of 2 to ``a``, 1 to ``b``, and 3 to the ancestor
of ``a`` and ``b``. Now, 


::

  (a[&dmv={3,1}]:1,b[&dmv=1]:1)[&dmv={2,1},dmt=3]


assigns linear population size to ``a``, which declines from 3 to 1 over the
branch of length 1. The root declines from 2 to 1 over 3 time units, as given by
the ``dmt`` attribute. In general, any continuous piecewise linear function can
be given using ``k`` values and ``k-1`` times (or ``k-2`` times, branch length
is used for the last time point). For example,

::

  (a[&dmv={4,3,2}, dmt={.2,1}]:1,b[&dmv=1]:1)[&dmv=2]

assigns to ``a`` a population size decreasing from 4 to 3 over [0,0.2], then
decreasing to 2 from 0.2 to 1 (end of branch). Note that population size
specifications are relative to the node, i.e. the node has the implicit time of
zero.
