============
Installation
============

Requirements
============


The SED fitter requires Python 2.5 or 2.6, and the following Python packages to be installed:

* `Numpy <http://numpy.scipy.org/>`_

* `Scipy <http://www.scipy.org/>`_

* `Matplotlib <http://matplotlib.sourceforge.net/>`_

.. _installation:

Installation
============


To obtain the SED fitting package, there are two options - either downloading a `stable' tar file, or using subversion to access the most up-to-date code.

Tar file installation
---------------------

Download the latest tar file from `here <http://target>`_, then install the package using::

    tar xvzf pysedfitter-x.x.x.tar.gz
    cd pysedfitter-x.x.x
    python setup.py install

Once the package is installed, you can safely remove the ``pysedfitter-x.x.x.tar.gz`` file and the ``pysedfitter-x.x.x`` directory.

SVN installation
----------------

Use `subversion <http://subversion.tigris.org/>`_ to check out the latest development version::

    svn co svn://caravan.astro.wisc.edu/pysedfitter  --username fitteruser

then install the package using::

   cd pysedfitter
   python setup.py install
   
Contact me if you do not know the svn password. The advantage of using subversion is that it is extremely easy to update the source code when bugs are fixed and features added.

Once the package is installed, you can safely remove the ``pysedfitter`` directory (although you may want to keep it for easy updating, see :ref:`updating`).

.. _updating:

Updating
========

To update an existing installation, just follow the instructions for :ref:`installation`. If you are using subversion, and you have kept the original ``pysedfitter`` directory, you can update by going to that directory and typing::

    svn update
    python setup.py install

Installing a set of models
==========================

To install a set of models (e.g YSO SEDs, stellar photospheres, etc.) simply download the appropriate tar file from the fitter webpage, and expand it in a central repository on your computer::

   cd /my/models/repository
   tar xvzf models_r06.tgz
