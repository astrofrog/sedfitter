from __future__ import print_function, division

import os
import glob

import numpy as np

from astropy.io import fits
from astropy.logger import log
from astropy import units as u
from astropy.utils.console import ProgressBar

from ..convolved_fluxes import ConvolvedFluxes
from ..sed import SED
from ..models import load_parameter_table
from .. import six


def convolve_model_dir(model_dir, filters, overwrite=False):
    """
    Convolve all the model SEDs in a model directory

    Parameters
    ----------
    model_dir: str
        The path to the model directory
    filters: list
        A list of Filter objects to use for the convolution
    overwrite: bool, optional
        Whether to overwrite the output files
    """

    for f in filters:
        if f.name is None:
            raise Exception("filter name needs to be set")
        if f.wavelength is None:
            raise Exception("filter central wavelength needs to be set")

    # Create 'convolved' sub-directory if needed
    if not os.path.exists(model_dir + '/convolved'):
        os.mkdir(model_dir + '/convolved')

    # Find all SED files to convolve
    sed_files = (glob.glob(model_dir + '/seds/*.fits.gz') +
                 glob.glob(model_dir + '/seds/*/*.fits.gz') +
                 glob.glob(model_dir + '/seds/*.fits') +
                 glob.glob(model_dir + '/seds/*/*.fits'))

    par_table = load_parameter_table(model_dir)

    if len(sed_files) == 0:
        raise Exception("No SEDs found in %s" % model_dir)
    else:
        log.info("{0} SEDs found in {1}".format(len(sed_files), model_dir))

    # Find out apertures
    first_sed = SED.read(sed_files[0])
    n_ap = first_sed.n_ap
    apertures = first_sed.apertures

    # Set up convolved fluxes
    fluxes = [ConvolvedFluxes(model_names=np.zeros(len(sed_files), dtype='U30' if six.PY3 else 'S30'), apertures=apertures, initialize_arrays=True) for i in range(len(filters))]

    # Set up list of binned filters
    binned_filters = []
    binned_nu = None

    # Loop over SEDs

    b = ProgressBar(len(sed_files))

    for im, sed_file in enumerate(sed_files):

        log.debug('Convolving {0}'.format(os.path.basename(sed_file)))

        # Read in SED
        s = SED.read(sed_file, unit_freq=u.Hz, unit_flux=u.mJy, order='nu')

        # Check if filters need to be re-binned
        try:
            assert binned_nu is not None
            np.testing.assert_array_almost_equal_nulp(s.nu.value, binned_nu.value, 100)
        except AssertionError:
            log.info('Rebinning filters')
            binned_filters = [f.rebin(s.nu) for f in filters]
            binned_nu = s.nu

        b.update()

        # Convolve
        for i, f in enumerate(binned_filters):

            fluxes[i].wavelength = f.wavelength
            fluxes[i].apertures = apertures
            fluxes[i].model_names[im] = s.name

            if n_ap == 1:
                fluxes[i].flux[im] = np.sum(s.flux * f.r)
                fluxes[i].error[im] = np.sqrt(np.sum((s.error * f.r) ** 2))
            else:
                fluxes[i].flux[im, :] = np.sum(s.flux * f.r, axis=1)
                fluxes[i].error[im] = np.sqrt(np.sum((s.error * f.r) ** 2, axis=1))

    for i, f in enumerate(binned_filters):
        fluxes[i].sort_to_match(par_table['MODEL_NAME'])
        fluxes[i].write(model_dir + '/convolved/' + f.name + '.fits',
                        overwrite=overwrite)
