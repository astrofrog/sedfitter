========================
Available Model Packages
========================

Below you will find links to existing model packages that can be downloaded.

Robitaille et al. (2006) YSO SED Models
=======================================

These are models that were published in `Robitaille et al. (2006)
<http://adsabs.harvard.edu/abs/2006ApJS..167..256R>`_.

.. important::

  **Please read this  disclaimer before downloading the models!**

  You are responsible for reading about the models and ensuring that you
  understand all the limitations inherent to YSO SED models in general, and
  this specific set of models. A number of caveats related to the models are
  listed in Section 2.2.4 of `Robitaille et al. (2006)
  <http://adsabs.harvard.edu/abs/2006ApJS..167..256R>`_ and a number of more
  general limitations are also discussed in `Robitaille (2008)
  <http://adsabs.harvard.edu/abs/2008ASPC..387..290R>`_.

  In particular, it is *very* important to not over-interpret results carried
  out by fitting these models. The models are simple axisymmetric 2-D models,
  which are likely only a very rough approximation of reality, and the models
  should only be used accordingly. The :math:`\chi^2` value is thus not
  meaningful in the traditional statistical sense of being linked to a
  well-defined probability. The models also include strong correlations
  between parameters due to the way the parameters were sampled (this is
  described in more detail in the two publications linked above) and can
  therefore not be reliably used to search for real trends in datasets.

  An example of a study using the models in the way they were intended (as
  rough probes of evolutionary stage), see `Forbrich et al. (2010)
  <http://adsabs.harvard.edu/abs/2010ApJ...716.1453F>`_. In addition, an
  example of a study of the reliability of SED modeling for protostars with
  these models is presented in `Offner et al. (2012)
  <http://adsabs.harvard.edu/abs/2012ApJ...753...98O>`_.

*By downloading the models, you acknowledge that you have read and agree with the above disclaimer*

**Download:** `models_r06_17jun08.tgz <ftp://ftp.astro.wisc.edu/outgoing/tom/model_packages/models_r06_17jun08.tgz>`_ [5.3Gb]

.. note:: the default ``parameters.fits`` file contained in the above model
          package does not include the envelope mass. If you are interested
          in this parameter, you will need to download
          `this <ftp://ftp.astro.wisc.edu/outgoing/tom/model_packages/parameters.fits.gz>`_
          file, decompress it, then put it inside the ``models_r06``
          directory and make sure it is named ``parameters.fits``.

Castelli & Kurucz (2004) stellar photosphere models
===================================================

These are the `Castelli and Kurucz (2004)
<http://arxiv.org/abs/astro-ph/0405087>`_ models downloaded from `here
<http://kurucz.harvard.edu/grids.html>`_ and bundled into the model package
format as a convenience.

.. important::

  **Please read this  disclaimer before downloading the models!**

  These models are provided as a convenience, and no warranty is made on the
  accuracy of the models or of the conversion into model package format. It is
  your responsability to ensure that you understand the models and their
  limitations.

*By downloading the models, you acknowledge that you have read and agree with the above disclaimer*

**Download:** `models_kurucz_05sep11.tgz <ftp://ftp.astro.wisc.edu/outgoing/tom/model_packages/models_kurucz_05sep11.tgz>`_ [87Mb]