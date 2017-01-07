from __future__ import print_function, division

# Still to implement:
# - Performance monitoring
# - Remove resolved models
# - Optional FITS input/output
# - Output convolved fluxes

import numpy as np
from astropy import units as u

from . import timer

from .models import Models
from .source import Source
from .utils import io
from .utils.validator import validate_array
from . import six
from .fit_info import FitInfoFile


class Fitter(object):
    """
    A fitter class that can be used to fit sources.

    This class is initialized using a particular set of models, and with
    specific fit parameters. It can then be used to fit data given by Source
    instances, and returns a FitInfo instance. Once initialized, the fit
    parameters cannot be changed, because changing most of them would require
    re-reading the models from disk.

    Parameters
    ----------
    filter_names : tuple or list
        List of filter names (given as individual strings) for which the data
        is defined. The filter names should be the name of the files in the
        ``convolved`` directory for the models, without the extensions. This is
        typically ``2J``, ``I1``, ``M1``, etc. You can also specify the
        wavelength as a :class:`~astropy.units.quantity.Quantity` instance
        instead of a filter name, and this will indicate that the SED fluxes
        closest to the requested wavelength should be used in the fitting.
    apertures : :class:`~astropy.units.quantity.Quantity` array instance
        The aperture radii that the data is specified in (as an angle). The
        fluxes may not be measured from aperture photometry, but this is meant
        to give an indication of the sizescale of the emission, and can be used
        to reject models that would have been clearly resolved at the distance
        specified.
    models_dir : str
        Name of the directory containing the models to use.
    extinction_law : :class:`~sedfitter.extinction.Extinction` instance
        The extinction law to use.
    av_range : tuple
        Minimum and maximum Av to allow in the fitting.
    distance_range : :class:`~astropy.units.quantity.Quantity` array instance
        Minimum and maximum distance to allow in the fitting in units of length.
    remove_resolved : bool, optional
        If set, then models larger than the aperture are removed. See
        Robitaille et al. (2007) for a discussion of this criterion.
    """

    def __init__(self, filter_names, apertures, model_dir,
                 extinction_law=None, av_range=None, distance_range=None,
                 remove_resolved=False):

        validate_array('apertures', apertures, domain='positive', ndim=1, physical_type='angle')
        validate_array('distance_range', distance_range, domain='positive', ndim=1, shape=(2,), physical_type='length')

        if len(apertures) != len(filter_names):
            raise ValueError("length of apertures list should match length of filter names list")

        # Construct filters dictionary
        self.filters = []
        for i in range(len(apertures)):
            filt = {'aperture_arcsec': apertures[i].to(u.arcsec).value}
            if isinstance(filter_names[i], six.string_types):
                filt['name'] = filter_names[i]
            elif isinstance(filter_names[i], u.Quantity):
                filt['wav'] = filter_names[i]
            else:
                raise ValueError("filter should be a string or a Quantity")

            self.filters.append(filt)

        # Read in models
        self.models = Models.read(model_dir, self.filters, distance_range=distance_range, remove_resolved=remove_resolved)

        # Add wavelength to filters
        for i, f in enumerate(self.filters):
            if 'wav' not in f:
                f['wav'] = self.models.wavelengths[i]

        # Set Av law
        self.av_law = extinction_law.get_av(self.models.wavelengths)

        # Set scale model - make this a scalar
        self.sc_law = -2. * np.ones(self.av_law.shape)

        self.model_dir = model_dir
        self.av_range = av_range
        self.extinction_law = extinction_law

    def fit(self, source):
        """
        Fit the specified source.

        Parameters
        ----------
        source : `~sedfitter.source.Source`
            The source to fit.

        Returns
        -------
        fit_info : `sedfitter.fit_info.FitInfo`
            The results of the fit.
        """

        info = self.models.fit(source, self.av_law, self.sc_law,
                               self.av_range[0], self.av_range[1])

        info.meta.model_dir = self.model_dir
        info.meta.filters = self.filters
        info.meta.extinction_law = self.extinction_law

        return info


def fit(data, filter_names, apertures, model_dir, output, n_data_min=3,
        extinction_law=None, av_range=None, distance_range=None,
        output_format=('F', 6.), output_convolved=False,
        remove_resolved=False):
    """
    Fit a set of sources with models.

    Parameters
    ----------
    data : str
        Filename of the file containing the data, one source per line (see
        documentation for a description of the required format).
    filter_names : tuple or list
        List of filter names (given as individual strings) for which the data
        is defined. The filter names should be the name of the files in the
        ``convolved`` directory for the models, without the extensions. This is
        typically ``2J``, ``I1``, ``M1``, etc. You can also specify the
        wavelength as a :class:`~astropy.units.quantity.Quantity` instance
        instead of a filter name, and this will indicate that the SED fluxes
        closest to the requested wavelength should be used in the fitting.
    apertures : :class:`~astropy.units.quantity.Quantity` array instance
        The aperture radii that the data is specified in (as an angle). The
        fluxes may not be measured from aperture photometry, but this is meant
        to give an indication of the sizescale of the emission, and can be used
        to reject models that would have been clearly resolved at the distance
        specified.
    models_dir : str
        Name of the directory containing the models to use.
    output : str
        Name of the file to output the fit information to (in binary format).
    extinction_law : :class:`~sedfitter.extinction.Extinction` instance
        The extinction law to use.
    av_range : tuple
        Minimum and maximum Av to allow in the fitting.
    distance_range : :class:`~astropy.units.quantity.Quantity` array instance
        Minimum and maximum distance to allow in the fitting in units of length.
    n_data_min : int, optional
        The minimum number of points a source needs to be fit.
    output_format : tuple, optional
        Tuple specifying which fits should be output. See the documentation
        for a description of the tuple syntax.
    output_convolved : bool, optional
        Whether to output the convolved fluxes (necessary if the convolved
        model fluxes are needed for the SED plot).
    remove_resolved : bool, optional
        If set, then models larger than the aperture are removed. See
        Robitaille et al. (2007) for a discussion of this criterion.
    """

    fitter = Fitter(filter_names, apertures, model_dir,
                    extinction_law=extinction_law, av_range=av_range,
                    distance_range=distance_range,
                    remove_resolved=remove_resolved)

    print(" ------------------------------------------------------------")
    print("  => Fitting parameters")
    print(" ------------------------------------------------------------")
    print("")
    print("   Minimum A_V      : %9.3f mag" % av_range[0])
    print("   Maximum A_V      : %9.3f mag" % av_range[1])
    print("   Minimum distance : %9.3f %s" % (distance_range[0].value, distance_range.unit))
    print("   Maximum distance : %9.3f %s" % (distance_range[1].value, distance_range.unit))
    print("")
    print(" ------------------------------------------------------------")
    print("  => Output parameters")
    print(" ------------------------------------------------------------")
    print("")
    print("   File   : %s" % output)
    print("   Format : %s" % output_format[0])
    print("   Number : %g" % output_format[1])
    print("")
    print(" ------------------------------------------------------------")
    print("  => Data format parameters")
    print(" ------------------------------------------------------------")
    print("")
    print("   Number of filters :  %i" % len(filter_names))
    print("")

    # Open datafile
    if isinstance(data, six.string_types):
        data_file = open(data, 'r')
    else:
        data_file = data

    print('')
    print('     Filter    Wavelength    Aperture (")   ')
    print('    ----------------------------------------')
    for f in fitter.filters:
        print('       %5s   %9.2f  %9.2f        ' % (f.get('name', ''), f['wav'].to(u.micron).value, f['aperture_arcsec']))
    print('')

    # Cycle through sources

    io.delete_file(output)

    fout = FitInfoFile(output, 'w')

    s = Source()

    t = timer.Timer()

    while True:

        try:
            s = Source.from_ascii(data_file.readline())
        except EOFError:
            break

        if s.n_data >= n_data_min:

            info = fitter.fit(s)

            if not output_convolved:
                info.model_fluxes = None

            info.keep(output_format)

            fout.write(info)

            t.display()

    t.display(force=True)

    fout.close()
