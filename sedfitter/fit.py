# Still to implement:
# - Min/max A_v - yes
# - Min/max scale - won't implement?
# - Minimum number of datapoints - done
# - Performance monitoring
# - YSO version - done
# - Output convolved fluxes
# - Option to removed resolved models
# - Optimize output -> go back to using FITS file, at least guarantee perfect backwards compatibility (even if it means performance decreases)

import numpy as np

import parfile
import source
import timer

from models import Models
from extinction import Extinction
from fit_info import FitInfo

import cPickle as pickle


def fit(parameter_file):

    # Read in fitting parameters
    par = parfile.read(parameter_file, 'par')

    # Read in data
    sources = source.read_sources(par['dfile'], n_min_valid=par['drequ'])

    # Read in data format
    f = file(par['dform'], 'rb')
    f.readline()
    nwav = int(f.readline().split('=')[0])
    filter_names = f.readline().split('=')[0].strip().split()
    apertures = np.array(f.readline().split('=')[0].strip().split(), float)

    # Read in model parameters
    modpar = parfile.read("%s/models.conf" % par['modir'], 'conf')

    if modpar['aperture_dependent']:
        if 'mind' in par and 'maxd' in par:
            n_distances = 1 + (np.log10(par['maxd']) - np.log10(par['mind'])) / modpar['logd_step']
            distances = np.logspace(np.log10(par['mind']), np.log10(par['maxd']), n_distances)
        else:
            raise Exception("For aperture-dependent models, mind/maxd are required")
    else:
        distances = None

    # Construct filters dictionary
    filters = []
    for i in range(nwav):
        filters.append({'aperture_arcsec': apertures[i], 'name': filter_names[i]})

    # Read in models
    models = Models(par['modir'], filters, distances=distances)

    # Construct filters dictionary
    filters = []
    for i in range(nwav):
        filters.append({'wav': models.wavelengths[i], 'ap': apertures[i], 'name': filter_names[i]})

    # Set extinction model
    extinction = Extinction(par['exlaw'])
    av_law = extinction.av(models.wavelengths)

    # Set scale model
    sc_law = -2. * np.ones(av_law.shape)

    # Cycle through sources

    t = timer.Timer()

    fout = file(par['ofile'], 'wb')
    pickle.dump(par['modir'], fout, 2)
    pickle.dump(filters, fout, 2)
    pickle.dump(models.names, fout, 2)
    pickle.dump(par['exlaw'], fout, 2)

    for s in sources:

        if s.n_data > int(par['drequ']):

            info = FitInfo(s)
            info.av, info.sc, info.chi2 = models.fit(s, av_law, sc_law, par['minav'], par['maxav'])

            info.sort()
            info.keep(par['oform'], par['onumb'])

            pickle.dump(info, fout, 2)

            t.display()

    t.display(force=True)

    fout.close()
