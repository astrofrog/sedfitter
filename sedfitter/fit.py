# Still to implement:
# - Performance monitoring
# - Remove resolved models
# - Optional FITS input/output
# - Output convolved fluxes

import cPickle as pickle

import numpy as np

import timer

from models import Models
from fit_info import FitInfo
from source import Source


def fit(data, filter_names, apertures, model_dir, output, n_data_min=3,
        extinction_law=None, distance_range=None, av_range=None,
        output_format=('F', 6.), output_convolved=False,
        remove_resolved=False):

    # Open datafile
    data_file = file(data, 'rb')

    # Construct filters dictionary
    filters = []
    for i in range(len(apertures)):
        filters.append({'aperture_arcsec': apertures[i], 'name': filter_names[i]})

    # Read in models
    models = Models(model_dir, filters, distance_range=distance_range, remove_resolved=remove_resolved)

    # Add wavelength to filters
    for i, f in enumerate(filters):
        f['wav'] = models.wavelengths[i]

    # Set Av law
    av_law = extinction_law.av(models.wavelengths)

    # Set scale model - make this a scalar
    sc_law = -2. * np.ones(av_law.shape)

    # Cycle through sources

    t = timer.Timer()

    fout = file(output, 'wb')
    pickle.dump(model_dir, fout, 2)
    pickle.dump(filters, fout, 2)
    pickle.dump(models.names, fout, 2)
    extinction_law.write_binary(fout)

    # NOTE _ CAN USE PROTOCOL 2 IN SOURCE FOR COORINDATES

    s = Source()

    while True:

        try:
            s.read_ascii(data_file)
        except:
            break

        if s.n_data > n_data_min:

            info = FitInfo(source=s)
            info.av, info.sc, info.chi2 = models.fit(s, av_law, sc_law, av_range[0], av_range[1])

            info.sort()
            info.keep(output_format[0], output_format[1])

            info.write_binary(fout)

            t.display()

    t.display(force=True)

    fout.close()
