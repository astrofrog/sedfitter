============
Installation
============

Requirements
============


The SED fitter requires Python 2.5 or 2.6, and the following Python packages to be installed:

* `Numpy <http://numpy.scipy.org/>`_

* `Scipy <http://www.scipy.org/>`_

* `Matplotlib <http://matplotlib.sourceforge.net/>`_

First-time Python users may want to consider an all-in-one Python installation package, such as the `Enthought Python Distribution <http://www.enthought.com/products/getepd.php>`_ which contains all three of the above packages (and many more).

.. _installation:

Installation
============


To obtain the SED fitting package, there are two options - either downloading a `stable' tar file, or using subversion to access the most up-to-date code.

Tar file installation
---------------------

Download the latest tar file from [not yet available], then install the package using::

    tar xvzf sedfitter-x.x.x.tar.gz
    cd sedfitter-x.x.x
    python setup.py install

Once the package is installed, you can safely remove the ``sedfitter-x.x.x.tar.gz`` file and the ``sedfitter-x.x.x`` directory.

SVN installation
----------------

Use `subversion <http://subversion.tigris.org/>`_ to check out the latest development version::

    svn co svn://caravan.astro.wisc.edu/pysedfitter/sedfitter  --username fitteruser

then install the package using::

   cd sedfitter
   python setup.py install
   
Contact me if you do not know the svn password. The advantage of using subversion is that it is extremely easy to update the source code when bugs are fixed and features added.

Once the package is installed, you can safely remove the ``sedfitter`` directory (although you may want to keep it for easy :ref:`updating`).

.. _updating:

Updating
========

To update an existing installation, just follow the :ref:`installation` instructions. If you are using subversion, and you have kept the original ``sedfitter`` directory, you can update by going to that directory and typing::

    svn update
    python setup.py install

Installing a set of models
==========================

To install a set of models (e.g YSO SEDs, stellar photospheres, etc.) simply download the appropriate tar file from the fitter webpage, and expand it in a central repository on your computer::

   cd /my/models/repository
   tar xvzf models_r06.tgz
