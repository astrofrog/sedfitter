# Still to implement:
# - Min/max A_v
# - Minimum number of datapoints
# - Performance monitoring
# - YSO version

import numpy as np

import parfile
from models import Models
import source
from extinction import Extinction
from output import FitInfo, OutputFile
import timer


def fit(parameter_file):

    # Read in fitting parameters
    par = parfile.read(parameter_file)

    # Read in data
    sources = source.read_sources(par['dfile'])

    # Read in data format
    f = file(par['dform'], 'rb')
    f.readline()
    nwav = int(f.readline().split('=')[0])
    filter_names = f.readline().split('=')[0].strip().split()
    apertures = np.array(f.readline().split('=')[0].strip().split(), float)

    # Read in models
    models = Models(par['modir'], filter_names)

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

    o = OutputFile()

    for s in sources:

        if s.n_data > int(par['drequ']):

            info = FitInfo(s)
            info.av, info.sc, info.chi2 = models.fit(s, av_law, sc_law)
            info.sort()
            info.keep(par['oform'], par['onumb'])

            o.append(info)

            t.display()

    t.display(force=True)

    o.write(par, filters, models.names)
