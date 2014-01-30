from __future__ import print_function, division

import os
import glob

import numpy as np
from astropy.io import fits
from astropy.logger import log

from ..source import Source
from .. import six


def list2data(file_in, file_out, models_dir=None, source_name='source'):
    """
    Convert a list-formatted data file into the standard fitting format

    In some cases it can be convenient to specify an SED in a list-based
    format with one source per file, one wavelength per line. This script
    converts the former into the latter.

    Parameters
    ----------
    file_in: str or file
        The input file.
    file_out: str or file
        The output file.
    models_dir: str, optional
        The models directory that will be used to fit the source. This can be
        useful to check that the desired filters exist, and to find out which
        monochromatic filters can be used.
    source_name: str, optional
        The source name to use. Defaults to 'source'.

    Returns
    -------
    filter_names: list of str
        The list of filters for the datafile (can be passed to fit())
    apertures: list of float
        The list of apertures for the datafile (can be passed to fit())
    """

    # Load datafile
    datafile = np.loadtxt(file_in, dtype=[('flux_type', int), ('name', 'S30'),
                                          ('valid', int),
                                          ('flux', float), ('err', float),
                                          ('aperture', float)])

    # Read in monochromatic wavelengths if available
    if models_dir is not None:

        log.info('Reading in monochromatic wavelengths')

        filter_files = glob.glob(os.path.join(models_dir, 'convolved/MO*'))

        mono_wav = np.zeros(len(filter_files), dtype=float)
        mono_filter = np.zeros(len(filter_files), dtype='S30')

        for i, filter_file in enumerate(filter_files):

            log.debug('Reading {}'.format(filter_file))

            mono_wav[i] = fits.getheader(filter_file)['FILTWAV']
            mono_filter[i] = os.path.basename(filter_file).replace('.fits', '').replace('.gz', '')

    elif np.any(datafile['flux_type'] == 2):

        raise Exception("list contains monochromatic fluxes, so need " +
                        "models_dir to be specified")

    # Create source
    s = Source()
    s.name = source_name
    s.x = 0.
    s.y = 0.
    s.n_wav = len(datafile)
    s.valid = datafile['valid']
    s.flux = datafile['flux']
    s.error = datafile['err']
    if isinstance(file_out, six.string_types):
        s.write_ascii(open(file_out, 'wb'))
    else:
        s.write_ascii(file_out)

    # Determine filter names
    for i in range(len(datafile)):
        if datafile['flux_type'][i] == 2:
            wav = float(datafile['name'][i])
            j = np.argmin(np.abs(mono_wav - wav))
            datafile['name'][i] = mono_filter[j]

    return datafile['name'].tolist(), datafile['aperture'].tolist()
