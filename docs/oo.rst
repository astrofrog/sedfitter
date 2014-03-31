Object-oriented fitting interface
=================================

In addition to fitting via the :func:`~sedfitter.fit` function, it is possible
to make use of the :class:`~sedfitter.fit.Fitter` class to load in all the
models and set the fitting parameters, and then fit sources on-demand (this can
be useful for GUI tools for example). To use the :class:`~sedfitter.fit.Fitter`
class, first instantiate the class by setting up the fit parameters, similarly
to the :func:`~sedfitter.fit` function, except that no data file or output file
is passed to the :class:`~sedfitter.fit.Fitter` class::

    from astropy import units as u
    from sedfitter import Fitter
    from sedfitter.extinction import Extinction

    # Define path to models
    model_dir = '/Volumes/Data/models/models_r06'

    # Read in extinction law)
    extinction = Extinction.from_file('kmh94.par', columns=[0, 3],
                                      wav_unit=u.micron, chi_unit=u.cm**2 / u.g)

    # Define filters and apertures
    filters = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
    apertures = [3., 3., 3., 3., 3., 3., 3.] * u.arcsec

    # Run the fitting
    fitter = Fitter(filters, apertures, model_dir,
                    extinction_law=extinction,
                    distance_range=[1., 2.] * u.kpc,
                    av_range=[0., 40.])
                    
To fit sources, you will then need to set up :class:`~sedfitter.source.Source`
instances. You can either set up a :class:`~sedfitter.source.Source` instance
with a single data file line::

    from sedfitter.source import Source
    s = Source.from_ascii('source_2 0.0 0.0 1 1 1 0.2 0.05 1.2 0.1 1.8 0.3')
    
or by instantiating the class and setting the parameters manually::    

    s = Source()
    s.name = ...
    s.x = ...
    
Once you have a :class:`~sedfitter.source.Source` instance, you can pass it to
the :meth:`~sedfitter.fit.Fitter.fit` method::

    info = fitter.fit(s)
    
The returned object is a :class:`~sedfitter.fit_info.FitInfo` instance which
can be passed instead of a filename to all the post-processing functions (such
as :func:`~sedfitter.plot`).
