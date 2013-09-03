from __future__ import print_function, division

# Still to implement:
# - Performance monitoring
# - Remove resolved models
# - Optional FITS input/output
# - Output convolved fluxes

try:
    import cPickle as pickle
except:
    import pickle

import numpy as np

from . import timer

from .models import Models
from .fit_info import FitInfo
from .source import Source
from . import util
from . import six


def fit(data, filter_names, apertures, model_dir, output, n_data_min=3,
        extinction_law=None, av_range=None, distance_range=None,
        output_format=('F', 6.), output_convolved=False,
        remove_resolved=False):
    '''
    Fit a set of sources with models

    Required Arguments:

        *data*: [ string ]
            Filename of the file containing the data, one source per line.
            (see :ref:`dataformat` for details about the format)

        *filter_names*: [ tuple or list ]
            List of the filters that the data is specified in.

        *apertures*: [ tuple or list ]
            List of the aperture radii that the data is specified in.

        *models_dir*: [ string ]
            Name of the directory containing the models to use.

        *output*: [ string ]
            Name of the file to output the fit information to (in binary format).

    Required Keyword Arguments:

        *extinction_law*: [ Extinction instance ]
            The extinction law to use

        *av_range*: [ tuple ]
            Tuple of the minimum and maximum Av to allow in the fitting

    Optional Arguments:

        *n_data_min*

    '''

    print(" ------------------------------------------------------------")
    print("  => Fitting parameters")
    print(" ------------------------------------------------------------")
    print("")
    print("   Minimum A_V      : %9.3f mag" % av_range[0])
    print("   Maximum A_V      : %9.3f mag" % av_range[1])
    print("   Minimum distance : %9.3f kpc" % distance_range[0])
    print("   Maximum distance : %9.3f kpc" % distance_range[1])
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
        data_file = file(data, 'rb')
    else:
        data_file = data

    # Construct filters dictionary
    filters = []
    for i in range(len(apertures)):
        filters.append({'aperture_arcsec': apertures[i], 'name': filter_names[i]})

    # Read in models
    models = Models(model_dir, filters, distance_range=distance_range, remove_resolved=remove_resolved)

    # Add wavelength to filters
    for i, f in enumerate(filters):
        f['wav'] = models.wavelengths[i]

    print('')
    print('     Filter    Wavelength    Aperture (")   ')
    print('    ----------------------------------------')
    for f in filters:
        print('       %5s   %9.2f  %9.2f        ' % (f['name'], f['aperture_arcsec'], f['wav']))
    print('')

    # Set Av law
    av_law = extinction_law.av(models.wavelengths)

    # Set scale model - make this a scalar
    sc_law = -2. * np.ones(av_law.shape)

    # Cycle through sources

    util.delete_file(output)

    fout = file(output, 'wb')
    pickle.dump(model_dir, fout, 2)
    pickle.dump(filters, fout, 2)
    extinction_law.write_binary(fout)

    # NOTE _ CAN USE PROTOCOL 2 IN SOURCE FOR COORINDATES

    s = Source()

    t = timer.Timer()

    while True:

        try:
            s.read_ascii(data_file)
        except:
            break

        if s.n_data >= n_data_min:

            info = FitInfo(source=s)
            info.av, info.sc, info.chi2, info.model_name, model_fluxes = models.fit(s, av_law, sc_law, av_range[0], av_range[1])

            if output_convolved:
                info.model_fluxes = model_fluxes

            info.sort()
            info.keep(output_format[0], output_format[1])

            info.write_binary(fout)

            t.display()

    t.display(force=True)

    fout.close()
