==============================
Traditional command-line tools
==============================

During the installation of the SED fitter, a number of executables are added to the ``bin/`` directory of the Python installation, and are thus available anywhere on the command-line.

.. note::
   Since the commands are available anywhere on the system regardless of
   directory, general names such as ``fit`` and ``plot`` (which were used in
   the Fortran version of the fitter) have been renamed to have the ``sed_``
   prefix. Thus, ``fit`` is not ``sed_fit``, ``plot`` is now ``sed_plot``, and
   so on.
   
The command-line interface to the SED fitting tool uses parameter files (identical to the Fortran version) to describe data, fitting, and plotting parameters. These are described in the following sub-sections.

Fitting
=======

.. _datadescr:

The data description file
-------------------------

The data description file is one containing information about the data you wish to fit, i.e. the number of wavelengths, and the filters and apertures the fluxes were measured in.

An example of such a file is the following::

   ############### DATA-DEPENDANT PARAMETERS ###############
   7                       = Number of wavelengths
   2J 2H 2K I1 I2 I3 I4    = Filters (format: xx xx xx xx xx)
   3. 3. 3. 3. 3. 3. 3.    = Aperture radii in arcseconds

This file contains the parameters one would need for a datafile containing 2MASS and *Spitzer*/IRAC fluxes. The first line is just used for comments and is ignored. For the subsequent lines, the information to the right of the ``=`` sign is ignored. The second line should contain the number of wavelengths/filters specified. The third line should be used to specify the filters in which the data are given, in the same order as in the data file. Filters are specified by two or more characters separated by a space. For a list of all the filters available, see :ref:`filters`. The fourth line gives the aperture radii that were used to compute the fluxes (in arcseconds).

The aperture information can be used for example to eliminate SED models of YSOs which would have been well resolved in a given aperture when using the aperture-dependent SED fitter. If this option is set in the fitter parameter file (see :ref:`fitterpar`), then models where :math:`r_{1/2}\,>\,{\rm aperture}` in at least one aperture are discarded, where :math:`r_{1/2}` is the outermost radius at which the surface brightness falls to 50% of the peak flux. This is because if the surface brightness of a model is half of the peak surface brightness at the edge of the aperture, then this means that it is clearly resolved with respect to the aperture size. For these reasons, the aperture size should be specified to the best of one's knowledge. If the flux was measured with an irregular aperture, then the aperture should correspond roughly to half of the maximum length of the aperture, and if the flux was measured using PSF photometry, it should be set to a few times the half-width at half maximum.

.. _parnotes:

Note on the fitting and plotting parameter files
------------------------------------------------

The fitter and plotting parameter files is made up of lines with syntax::

    value = name = comment

When reading in parameters, the various command looks for the required
parameters based on the ``name``, so it is possible to re-order the
parameters, to insert blank lines or comment lines.

One ubiquitous pair of parameters is the output format (``oform``) and number
(``onumb``), which consist of a single character and a real or integer number
respectively. Setting ``oform`` to ``N`` indicates that the number of fits to
select is given by ``onumb``. Setting ``oform`` to ``C``, ``D``, ``E``, or
``F`` indicates that the number of fits should be selected according to the
:math:`\chi^2` value of the fits. The ``C`` option selects all fits with a
:math:`\chi^2` value below ``onumb``. The ``D`` option selects all fits with a
:math:`\chi^2-\chi^2_{\rm best}` value below ``onumb``. The ``E`` and ``F``
options are equivalent to the ``C`` and ``D`` options respectively, but with
the \chi^2 value given per datapoint.

The fitter parameter file
-------------------------

The fitter parameter file is one containing parameters for the fitting process. Two examples of such a file - ``fitter_stellar_ex.par`` and ``fitter_yso_ex.par` - are given in the ``examples`` directory in the SED fitting tool tar file. The first is an example of parameter file for aperture and scale-independent models (e.g. Kurucz stellar photospheres) and the latter is an example of parameter file for aperture and scale-dependent models  (e.g. R06 YSO models). The parameters required are:

* The filenames of the data parameter file and the input data file in fitter
  format::

    data_glimpse_ex.par     = dform = data format
    data_comparison         = dfile = data file

* The minimum number of datapoints that you require to be able to fit a
  source::

    4                       = drequ = minimum number of datapoints

  Any source with a number of valid datapoints less than this will be
  dropped, and will not be considered for fitting. I would suggest never to
  set this value any lower than ``3``.
  
* The directory containing the models for the fitting::

    models_r06              = modir = models directory

* The minimum and maximum :math:`A_{\rm V}` to use in the fitting::

    3.                      = minav = Minimum interstellar extinction
    22.                     = maxav = Maximum interstellar extinction

  These extinction values refer to interstellar extinction, and in the case of
  the YSO models do not take into account the extinction caused by the
  circumstellar environment (disks, envelopes, etc...) which is intrinsic to
  the models.
  
* If the aperture-dependent SED fitter is used, a distance range can be
  specified. This distance range (specified in kpc) is used to rule out any
  models which do not have the right luminosity to fit the data::

    2.                      = mind  = minimum distance in kpc
    2.5                     = maxd  = maximum distance in kpc

  These options are not required for the aperture-independent fitter as shown
  in the ``fitter_stellar_ex.par`` file provided.

* The extinction law to use.::

    ex_law_gl.par    = exlaw = extinction law

* Output options - the ``oform`` and ``onumb`` options specify the maximum
  number of fits to be stored in the fitter output file::
  
    output.fits             = ofile = output filename
    F                       = oform = what to output (N/C/P/D/E/F)   
    6.                      = onumb = number relating to the above option

  See :ref:`parnotes` for more details.


* Whether to output convolved fluxes to the output file (used for plotting the
  convolved fluxes with ``src/plot``::

    Y                       = oconv = Output convolved fluxes

* Whether the source is larger than the aperture at any wavelength::

    N                       = kpext = Any sources larger than aperture? (Y/N)

  If set to ``N``, models with :math:`r_{1/2}\,>\,{\rm aperture}` in at least
  one aperture are discarded (see :ref:`datadescr`). This option is not
  required for the aperture-independent fitter as shown in the
  ``fitter_stellar_ex.par`` file provided.

Running the fitter
------------------

To run the fitter, use the ``sed_fit`` or ``sed_fit_stellar`` command::

  Usage: fit [arguments]
  
  Required arguments :
    par_file=filename  - the parameter file (e.g. config/sn_requirements.par)
  
  Optional arguments :
    output=filename    - the output directory
    input=filename     - the input fitter data file
    models=directory   - the models to use
    best=yes/no        - whether to output only the best fit

For example if your parameter file is named ``fitter.par``, and you are fitting stellar photosphere models, type::

    sed_fit_stellar par_file=fitter.par

The fitter should now run! An output file will be created, containing information about the sources and the model fits. This output can be further processed with the programs described in the following sections.

To make it easier to write pipeline scripts to fit a region several times with different models, several options can be specified after the fit command. These options **override** the values set in the ``fitter.par`` file:

* ``output=filename`` tells the fitter what to call the output FITS file (the
  filename should end in ``.fits``)

* ``input=filename`` tells the fitter which input file to use

* ``models=directory`` tells the fitter which models to use

* ``best=yes`` tells the fitter to output only the best SED, effectively
  forcing ``oform=N`` and ``onumb=1``. Setting best=no defauts to the
  ``oform`` and ``onumb`` values in the fitter parameter file

You can use any combination of these. If you do not specify one of these options, the value in the parameter file will be used instead.

.. note::
    The output file is no longer in FITS format, and thus can no longer be
    viewed by external programs.

Processing the output file
==========================

Splitting the output file
-------------------------

Once the above output file has been produced, it is possible to split it to separate well and badly fit sources. To do this, the ``filter_output`` command can be used::

  Usage: filter_output [arguments]
  
  Required arguments :
    input=filename     - the input file
  
  then one and only one of the following:
    chi=value          - constrain using total chisquared
    cpd=value          - constrain using chisquared per datapoint
 
The various ways to split the output are:

* ``chi=value``: this means that any source whose best fit (on the basis of
  the :math:`\chi^2` value) has a total (not reduced) :math:`\chi^2` less than
  ``value`` will be considered as well fit, and any other source will be
  considered as badly fit.

* ``cpd=value``: this means that any source whose best fit (on the basis of
  the :math:`\chi^2` value) has a :math:`\chi^2` per data point value less
  than ``value`` will be considered as well fit, and any other source will be
  considered as badly fit.

Only one option should be specified.

This command will create two new files in the original output directory, one for the well-fitted sources, and one for the badly-fitted sources (indicated by the ``_good`` and ``_bad`` suffix).

Making a data file from an output file
--------------------------------------

A utility named ``fits2data`` is available to produce a data file suitable for input to the fitter from an output file. It is used as follows:

  Usage: fits2data [arguments]
  
  Required arguments :
    input=filename     - the input fitter FITS file
    output=filename    - the output data file
    
Plotting the results
====================

.. _plotcommon:

Overview
--------

.. note::
    ``plot_params_1d`` and ``plot_params_2d`` are not implemented to date.

There are three tools available to plot the results. ``plot`` produces plots of the SED fits, ``plot_params_1d`` produces histograms of parameter values, and ``plot_params_2d`` produces plots with one parameter on each axis (useful to distinguish real from artificial correlations between parameters). Two examples of plotting parameter files are provided in the ``example/`` directory: ``plot_stellar_ex.par`` is set to show only the best fit (useful for stellar photosphere models for example) and ``plot_yso_ex.par`` is set to show a large number of good fits, making use of the greyscale feature and interpolating the SED to different apertures at different wavelengths (useful for the YSO models for example).

Both example parameter files contains all the parameters required by the three plotting programs. Either of the files can technically be split into three different files, but there are advantages of not doing so, for example to make it easy to consistently change the number of model fits shown in the different plots (via the ``oform`` and ``onumb`` parameters). The parameters are split into \textbf{General} and \textbf{Advanced} parameters, and in each case, whether they apply to all programs, or one in particular.

The **general** options required by all programs are::

    N       = oform = what to plot (N/C/P/D/E/F)
    1000    = onumb = number relating to the oform option

These are the output options - the ``oform`` and ``onumb`` options specify the maximum number of fits to be plotted (see :ref:`parnotes`).

The **advanced** options required by all programs are::

    Y        = pname = show source name on plot (Y/N)
    .eps/VPS = devic = PGPLOT device (change to .eps/VCPS for color)
    0.75     = chlab = title and labels character height
    0.60     = chaxi = numerical axis character height  
    2        = lwbox = line width of box, numbers, and labels
    1        = lwsed = line width of SEDs
    0.75     = shade = shade of grey
    2.5      = labdx = x-label displacement
    3.0      = labdy = y-label displacement

These are very technical options used to change the appearance of the plots. The default values are sensible, but you can always tweak them if needed.

SED plotting with ``sed_plot``
------------------------------

This is the main plotting program which produced SED plots. The ``plot`` specific parameters are:

* Options to control the general type of plot::

    A       = pmode = mode - I(ndividual) or A(ll SEDs on one plot)
    Y       = pgrey = if in A mode, whether to plot a greyscale

  ``I`` mode means make one plot per source per fit. ``A`` mode means make one
  plot per source with all the fits on the same plot. For ``A`` mode, the best
  fit is shown in black, and the subsequent fits are shown in grey. However,
  this may produce ``eps`` files which take a while to draw - the ``pgrey``
  option fixes this by rasterizing all the grey SEDs into a single greyscale
  map, which makes loading almost instantaneous.

* Options to control the minimum and maximum values shown::

    A       = xmode = X values - A(uto)/M(anual)
    0.1     = xminm = Xmin (if manual)
    100.    = xmaxm = Xmax (if manual)
    1.      = xmina = Xmin margin in orders of magnitude (if auto)
    1.      = xmaxa = Xmax margin in orders of magnitude (if auto)

    A       = ymode = Y values (A(uto)/M(anual))
    1.e-14  = yminm = Ymin (if manual)
    1.e-9   = ymaxm = Ymax (if manual)
    1.      = ymina = Ymin margin in orders of magnitude (if auto)
    2.      = ymaxa = Ymax margin in orders of magnitude (if auto)

  These options control the minimum and maximum x and y values for the plots.
  If the ``xmode`` or ``ymode`` options are set to ``M`` for example, then the
  ``xmin`` and ``xmax`` or ``ymin`` and ``ymax`` values are used. If the
  options are set to ``A``, then the plotting tool automatically finds the
  minimum and maximum values from the data, and adds a margin specified in the
  above parameters. For example, if ``xmode=A``, and if the minimum datapoint
  is at 1&mu;m, and the ``xmina`` margin is set to ``1``, then the minimum
  wavelength on the plot will be 0.1&mu;m.

* Options to control what to plot::
 
    Y       = pseds = plot SEDs (Y/N)
    N       = pconv = plot convolved fluxes (Y/N)
    Y       = patmo = overplot stellar atmosphere (Y/N)
    interp  = stype = what type of SED to show (see manual)

    Y       = pinfo = show fit info on plot (Y/N)

    N       = pcapt = create color caption (Y/N)

  These are extra options which allow you to further customize the plots. You
  can either plot SEDs, or the convolved fluxes (requires ``oconv=Y`` in
  ``fitter.par``), or both, and you can show the stellar photosphere (if
  ``pmode=A`` only the stellar photosphere for the best fit is shown).The
  ``stype`` parameter sets what kind of SED to plot. The options are:

    * ``largest`` - the SED for the largest aperture specified for the data

    * ``largest+smallest`` - the SEDs for the largest and smallest apertures
      specified for the data

    * ``interp`` - an SED interpolated to the various apertures as a function
      of wavelength. This works best if there are no sudden jumps in the
      aperture as a function of wavelength.

    * ``all`` - the SEDs for all the different apertures with color coded
      lines. Note that this is no longer possible when ``pmode=A`` as this
      produced confusing plots.

  In the case where ``stype=all``, a color caption can be produced showing the
  color of each SED as a function of aperture. Finally, ``pinfo`` allows the
  :math:`\chi^2`, A$_{\rm V}$ and scalefactor to be overplotted on the SED
  plot.

* Options for customizing the look of the plots::

    400     = greyx = greyscale resolution (x direction)
    300     = greyy = greyscale resolution (y direction)
    2.5     = greyc = greyscale clipping value
    8.0     = greym = greyscale stretch

The ``plot`` command is used as follows::

  Usage: plot [arguments]
 
  Required arguments :
    par_file=filename  - the parameter file (e.g. plot.par)
    input=filename     - the input FITS file
    output=filename    - the output plots directory
 
Histogram plotting with ``plot_params_1d``
------------------------------------------

This plotting program can produce histograms of the parameters of a given number of YSO fits. How many fits are shown is controlled as before using the ``oform`` and ``onumb`` parameters.

.. note:
    All the fits which contribute to the histogram are weighted equally.

As well as the common options specified in :ref:`plotcommon` an advanced option is available to specify the number of bins in the histogram. This number should be a multiple of 10 for best results::

    30     = histn = number of bins in histogram 

This tool is used as::

  Usage: plot_param_1d [arguments]
  
  Required arguments :
    par_file=filename  - the parameter file (e.g. plot.par)
    input=filename     - the input fitter FITS file
    output=directory   - the output plots directory
    parameter=name     - the name of the parameter to show (e.g. MDISK)
    log=yes/no         - whether to plot the parameter on a log scale
    zero=yes/no        - if log=yes, whether to show zero values

For a list of possible ``parameter`` values for the R06 or the Kurucz models, see :ref:`parnames`. If you define your own set of models, whatever numerical parameters you include in ``parameters.fits.gz`` will be available here.

Two histograms are shown in each plot. In grey, the distribution of models in the grid is shown, normalized so that the maximum value is 1. The hashed histogram shows the distribution of the fits, also normalized to its maximum value (the normalization factor is not the same for the two histograms as the good fits only represent a small fraction of all the models in the grid).

2-D parameter plotting with ``plot_params_2d``
----------------------------------------------

This plotting program can produce two-dimensional maps of the parameters.  How many fits are shown is controlled as before using the ``oform`` and ``onumb`` parameters.

As well as the common options specified in \S\ref{sec:plotoverview} a few advanced options are available, although the default values should be sensible::

    8       = greyspl2d = greyscale resampling
    2048    = greybig2d = greyscale resolution (major direction)
    136     = greysma2d = greyscale resolution (minor direction)
    13.     = greycli2d = greyscale clipping value
    40.     = greymax2d = greyscale stretch
    0.75    = ch2dpoint = size of points in 2d plots

This tool is used as::

  Usage: plot_param_2d [arguments]
  
  Required arguments :
    par_file=filename  - the parameter file (e.g. plot.par)
    input=filename     - the input fitter FITS file
    output=directory   - the output plots directory
    parameterx=name    - the name of the parameter to show on the x axis
    parametery=name    - the name of the parameter to show on the y axis
    logx=yes/no        - whether to plot the x-axis parameter on a log scale
    logy=yes/no        - whether to plot the y-axis parameter on a log scale
    zerox=yes/no       - if logx=yes, whether to show zero values for the x-axis
    zeroy=yes/no       - if logy=yes, whether to show zero values for the y-axis

The distribution of models in the grid is shown in grey points. All the fits are shown as black points.

