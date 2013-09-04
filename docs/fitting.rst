Fitting
=======

In order to fit the data, the :func:`sedfitter.fit` function should be used as follows::

    from sedfitter import fit
    fit(data_file, filters, apertures, model_dir, output_file,
        extinction_law=..., av_range=..., distance_range=...)

The arguments are:

* ``data_file``: the filename of an ASCII data file containing the fluxes for
  the sources (see :doc:`data` for a description of the format)

* ``filters``: a list of filter names (given as individual strings) for which
  the data is defined. The filter names should be the name of the files in the
  ``convolved`` directory for the models, without the extensions. This is
  typically ``2J``, ``I1``, ``M1``, etc.

* ``apertures``: a list of floating point values giving the radii of the
  circular apertures for which the fluxes are defined. The fluxes may not be
  measured from aperture photometry, but this is meant to give an indication of
  the sizescale of the emission, and can be used to reject models that would
  have been clearly resolved at the distance specified.

* ``model_dir``: the path to the directory containing the models. This
  directory should have the standard layout described `here <broken_link>`_.

* ``output_file``: the name of the output file to write with the results of the
  fitting.

* ``extinction_law``: an :class:`~sedfitter.extinction.Extinction` instance.

* ``av_range``: a tuple of two values giving the lower and upper acceptable Av
  values.

* ``distance_range``: a tuple of two values giving the lower and upper
  acceptable distance values, in kpc.

In addition to the above compulstory arguments, the following optional arguments can be specified:

* ``output_format``: which models to output (see :doc:`select_syntax`)

* ``output_convolved``: whether or not to output convolved fluxes in the output
  file (specified as a boolean - default is ``False``)

* ``remove_resolved``: whether to remove models that are larger than the
  apertures in any of the bands (specified as a boolean - default is ``False``)

For more details, you can also check out the :func:`sedfitter.fit` documentation.

The following is a complete example showing how to use the :func:`sedfitter.fit` function::

    from sedfitter import fit
    from sedfitter.extinction import Extinction

    # Define path to models

    model_dir = '/Users/tom/Models/models_r06'

    # Read in extinction law

    extinction = Extinction.from_file('kmh94.par')

    # Define filters and apertures

    filters = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']

    apertures = [3., 3., 3., 3., 3., 3., 3.]

    # Run the fitting

    fit('data_glimpse', filters, apertures, model_dir,
        'output_glimpse.fitinfo',
        extinction_law=extinction,
        distance_range=[1., 2.],
        av_range=[0., 40.],
        output_format=('F', 3.),
        output_convolved=False,
        remove_resolved=True)


