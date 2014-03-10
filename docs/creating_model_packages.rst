.. _modelpackages:

=========================
Creating model packages
=========================

This page describes the format required of a *model package* in order to be used by the fitter. The format adopted for all files is the FITS standard. All FITS files presented here should have ``BITPIX=-32`` and ``EXTEND=T`` in the primary header. The model package should have a name starting with ``model_`` and containing the following:

* ``models.conf`` - a configuration file for the models.

* ``seds/`` - a directory containing all the model SEDs as ``fits`` or ``fits.gz`` files, in the format described below. In order to avoid huge numbers of files in a single directory, sub-directories can be created with a section of the model names. The size of this section can later be specified in the models :ref:`config`.

* ``convolved/`` - a directory containing all the convolved models, one file per filter. The name of the files should be ``filter.fits`` where ``filter`` should be replaced with the filter code (e.g. ``2J``, ``2H``, etc.).

* ``parameters.fits.gz`` - the parameters for the models, in the same order as the convolved fluxes

.. _config:

Configuration file
==================

The ``models.conf`` configuation file should be a plain-ASCII file containing the following keyword/value pairs:

* The name of the model package::

    name = YSO models

* The length of the sub-directories for SED files::

    length_subdir = 5

  In cases where hundreds of thousands of SEDs are present, it is recommended
  to split them into sub-directories inside the ``seds`` directory. This is
  done by taking the first few letters of all the models and putting all the
  models that have the same first few letters in the same sub-directory, with
  the name of these first letters. The above option is the length of the
  sub-string to use for this. For example, if ``length_subdir`` is set to 5,
  all models starting with ``30001`` should be in a directory called ``30001``
  inside ``seds/``.

* Whether the models are aperture/scale dependent::

    aperture_dependent = yes

  This should only be set to ``no`` if the models are not absolutely scaled to
  1kpc, and if they are only specified in one aperture. In this case, the
  ``distance_range`` option to :func:`~sedfitter.fit` is ignored.

* The step in the distance, in log space. Model SEDs are computed for a range
  of distances between ``mind`` and ``maxd``, uniformly spaced in log space,
  and separated by this value::

    logd_step = 0.025

SED files
=========

HDU 0
-----

The primary header must contain following keywords::

    VERSION  = 1   (integer) (the current version, in case we make big changes in future)
    MODEL    = value (character) (the model name)
    IMAGE    = T/F (logical) (whether images are available in HDU 0)
    WAVLGHTS = T/F (logical) (whether wavelengths are available in HDU 1)
    APERTURS = T/F (logical) (whether apertures are available in HDU 2)
    SEDS     = T/F (logical) (whether SEDs are available in HDU 3)

If ``WAVLGHTS = T``, then the header must contain::

    NWAV = value (integer) which is the number of wavelengths

If ``APERTURS = T``, then the header must contain::

    NAP = value (integer) which is the number of apertures

If ``IMAGE = T``, then this **HDU** should contain a data cube which is the model image as a function of wavelength. This will be specified in more detail in future.

HDU 1
-----

if ``WAVLGTHS = T``, this **HDU** should contain a 2-column binary table with
``NWAV`` values in each column. The first column should have the title
``WAVELENGTH``, format ``1E``, and unit ``MICRONS``. The second column should
have the title ``FREQUENCY``, format ``1E``, and unit ``HZ``. Optionally, a
third column can be included with the title ``STELLAR_FLUX``, format ``1E``,
and units specified appropriately. This column can be used to specify the
model stellar photosphere used to compute YSO models for example. Only ``mJy``
and ``ergs/cm\^2/s`` are supported as units at this time. If specified, the
stellar photosphere should have the correct scaling relative to the SEDs in
``HDU 3``.

HDU 2
-----

if ``APERTURS = T``, this **HDU** should contain a 1-column binary table with
``NAP`` values. The column should have the title ``APERTURE``, format ``1E``,
and units ``AU``. These are the apertures for which the SEDs in **HDU 3**
are tabulated.

HDU 3
-----

if ``SED = T``, this **HDU** should contain a binary table with at least one column, and **HDU 1** and  **HDU 2**  should also contain data. Each column should consist of ``NAP`` rows of real vectors with dimension ``NWAV``. Thus, each cell contains an SED. The format of each column should be ``nE``, where ``n=NWAV``. The title and units of each column should be specified. The columns can contains SEDs such as the total flux, the stellar flux, the disk flux, etc. or related uncertainties. The following column titles are examples of ones that can be used::

    TOTAL_FLUX
    TOTAL_FLUX_ERR
    STELLAR_FLUX
    STELLAR_FLUX_ERR
    DISK_FLUX
    DISK_FLUX_ERR
    ENVELOPE_FLUX
    ENVELOPE_FLUX_ERR
    DIRECT_FLUX
    DIRECT_FLUX_ERR
    SCATTERED_FLUX
    SCATTERED_FLUX_ERR
    THERMAL_FLUX
    THERMAL_FLUX_ERR
    etc.

The order of the columns is not important as there are ``FITS`` routines to search for a specific column.

.. note::
    The SED fitter requires a column ``TOTAL_FLUX`` to be present, and will
    return an error otherwise. Only ``mJy`` and ``ergs/cm^2/s`` are supported
    as units at this time.

Convolved fluxes file
=====================

HDU 0
-----

The primary header must contain following keywords::

    FILTWAV  = value (real) (the characteristic wavelength of the filter)
    NMODELS  = value (integer) (the number of models)
    NAP      = value (integer) (the number of apertures)

HDU 1
-----

This **HDU** should contain a 5-column binary table. The column titles should be::

    MODEL_NAME
    TOTAL_FLUX
    TOTALF_FLUX_ERR
    RADIUS_SIGMA_50
    RADIUS_CUMUL_99

The first column should have format ``30A`` and should contain the name of each model. No units are required. The second and third columns should have format ``nE`` where ``n=NAP``, with each cell containing a vector with the fluxes in the different apertures. The fourth and fifth column should have format ``1E`` and contain the outermost radius at which the surface brightness falls to 50% of the maximum surface brightness, and the radius inside which 99% of the flux is contained respectively. These two columns should have units ``AU``.

HDU 2
-----

This **HDU** should contain a 1-column binary table with ``NAP``
values. The column should have the title ``APERTURE``, format ``1E``, and units ``AU``. These are the apertures for which the fluxes in **HDU 1** are tabulated.

Model parameters
================

HDU 0
-----

The primary header must contain following keywords::

    NMODELS  = value (integer) (the number of models)

HDU 1
-----

This **HDU** should contain a binary table with the model parameters. Any number of columns can be included, in any order. Only parameters with format ``1E`` will be usable by the programs to plot parameters, but text parameters with format ``nA`` can also be included (e.g. dust model filenames, etc.). One column is compulsory, with title ``MODEL_NAME`` and format ``30A``. It should contain the same names as the convolved fluxes file, and in the same order.

