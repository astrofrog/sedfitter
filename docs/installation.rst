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

* `Astropy <http://www.astropy.org>`_ 0.3 or later

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

If you want to install the latest developer version, you will need to clone
the git repository::

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