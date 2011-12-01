import os
import glob

import pyfits
import numpy as np

from .convolved_fluxes import ConvolvedFluxes
from ..sed import SED
from ..logger import log


def convolve_model_dir(model_dir, filters):
    '''
    Convolve all the model SEDs in a model directory

    Parameters
    ----------
    model_dir: str
        The path to the model directory
    filters: list
        A list of Filter objects to use for the convolution
    '''

    for f in filters:
        if f.name is None:
            raise Exception("filters need to be named")

    # Find all SED files to convolve
    sed_files = glob.glob(model_dir + '/seds/*.fits.gz') + \
                glob.glob(model_dir + '/seds/*.fits')

    # Find out number of apertures
    n_ap = pyfits.getheader(sed_files[0], memmap=False)['NAP']

    # Set up convolved fluxes
    fluxes = [ConvolvedFluxes(n_models=len(sed_files), n_ap=n_ap) \
                for i in range(len(filters))]

    # Set up list of binned filters
    binned_filters = []
    binned_nu = None

    # Loop over SEDs
    for im, sed_file in enumerate(sed_files):

        log.debug('Convolving {}'.format(os.path.basename(sed_file)))

        # Read in SED
        s = SED()
        s.read(sed_file, unit_freq='Hz', unit_flux='mJy', order='nu')

        # Check if filters need to be re-binned
        if np.any(s.nu != binned_nu):
            log.info('Rebinning filters')
            binned_filters = [f.rebin(s.nu) for f in filters]
            binned_nu = s.nu

        # Convolve
        for i, f in enumerate(binned_filters):
            fluxes[i].model_names[im] = s.name
            if n_ap == 1:
                fluxes[i].flux[im] = np.sum(s.flux * f.r)
            else:
                raise NotImplemented()

    for i, f in enumerate(binned_filters):
        fluxes[i].write(model_dir + '/convolved/' + f.name + '.fits')
