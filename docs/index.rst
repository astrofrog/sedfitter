.. SED Fitter documentation master file, created by
   sphinx-quickstart on Tue Sep  3 17:16:10 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python SED Fitter documentation
===============================

This package is an experimental Python port of the SED Fitting tool described in
`Robitaille et al. (2007) <http://adsabs.harvard.edu/abs/2007ApJS..169..328R>`_
and is still under active development.

.. important::

  **No guarantees are made** about the accuracy of the results, and users are
  responsible for ensuring that the results are reasonable and that they
  understand all the limitations inherent to SED fitting.

Please report any issues with the package `here
<https://github.com/astrofrog/sedfitter/issues>`_ (not by email).

User documentation
------------------

.. toctree::
   :maxdepth: 1

   introduction.rst
   installation.rst
   fitting.rst
   convolution.rst
   oo.rst

Appendices
----------

.. toctree::
   :maxdepth: 1

   data.rst
   select_syntax.rst
   common_model_packages.rst
   creating_model_packages.rst
   filter_names.rst
   api/api.rst


