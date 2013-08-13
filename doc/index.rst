.. biopy documentation master file, created by
   sphinx-quickstart on Wed Sep 14 18:54:02 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Biopy documentation
===================

biopy is a small Python library for exploring and developing phylogenetic
ideas. biopy is useful for generating simulated data and performing
off-line analysis.

Most of the code is of general interest, but some is BEAST/\*BEAST specific.

biopy is **not** a high-duty, high-performance library, but some code
which proved painfully slow in pur Python has been converted to "C++".

The library requires Python 2.7 and the following packages:
`numpy <http://numpy.scipy.org//>`_, `scipy <www.scipy.org//>`_, lxml,
matplotlib and biopython.

If you use biopy in a publication, please cite it using
`<http://dx.doi.org/10.6084/m9.figshare.761224>`_.

.. toctree::
   :maxdepth: 5

   gen
   biopyscripts
   biopylib
      
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

