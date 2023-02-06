============
Installation
============

Requirements
============

The SED fitter requires Python 3.8 or later, and the following Python packages
to be installed:

* `Numpy <http://www.numpy.org>`_

* `Scipy <http://www.scipy.org>`_

* `Matplotlib <http://www.matplotlib.org>`_

* `Astropy <http://www.astropy.org>`_

First-time Python users may want to consider an all-in-one Python installation
package, such as the `Anaconda Python Distribution
<http://continuum.io/downloads>`_ which provides all of the above dependencies.

.. _installation:

Installation
============

You can install the SED fitting package using::

    pip install sedfitter

or you can get the latest tar file from `PyPI
<https://pypi.python.org/pypi/sedfitter>`_.

If you want to install the latest developer version, you can do:

    pip install git+https://github.com/astrofrog/sedfitter

Testing
=======

If you want to check that all the tests are running correctly with your Python
configuration, you can also run::

    pytest sedfitter

in the source directory. If there are no errors, you are good to go!
