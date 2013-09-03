============
Installation
============

Requirements
============


The SED fitter requires Python 2.6, 2.7, 3.2, or 3.3, and the following Python
packages to be installed:

* `Numpy <http://www.numpy.org>`_

* `Scipy <http://www.scipy.org>`_

* `Matplotlib <http://www.matplotlib.org>`_

* `ATpy <http://atpy.github.io>`_ 0.9.7 or later.

* `Astropy <http://www.astropy.org>`_ (this should be the latest developer
  version, as it will not work with Astropy 0.2 or earlier)

First-time Python users may want to consider an all-in-one Python installation
package, such as the `Anaconda Python Distribution
<http://continuum.io/downloads>`_ which provides Numpy, Scipy, and Matplotlib.

.. _installation:

Installation
============

To obtain the SED fitting package, you will need to clone the git repository::

    git clone http://github.com/astrofrog/sedfitter

You can then install the package with::

    cd sedfitter
    python setup.py install

Testing
=======

If you want to check that all the tests are running correctly with your Python
configuration, you can also run::

    python setup.py test

in the source directory. If there are no errors, you are good to go!    