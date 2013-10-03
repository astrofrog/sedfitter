from __future__ import print_function, division

import os
import glob

import numpy as np
from astropy.io import fits
from astropy.logger import log
from astropy.table import Table
from astropy import units as u

from ..convolved_fluxes import ConvolvedFluxes
from ..sed import SED


def convolve_model_dir_monochromatic(model_dir, overwrite=False, max_ram=8,
                                     wav_min=-np.inf, wav_max=np.inf):
    '''
    Convolve all the model SEDs in a model directory

    Parameters
    ----------
    model_dir: str
        The path to the model directory
    overwrite: bool, optional
        Whether to overwrite the output files
    max_ram: float, optional
        The maximum amount of RAM that can be used (in Gb)
    wav_min: float, optional
        The minimum wavelength to consider. Only wavelengths above this value
        will be output.
    wav_max: float, optional
        The maximum wavelength to consider. Only wavelengths below this value
        will be output.
    '''

    # Create 'convolved' sub-directory if needed
    if not os.path.exists(model_dir + '/convolved'):
        os.mkdir(model_dir + '/convolved')

    # Find all SED files to convolve
    sed_files = glob.glob(model_dir + '/seds/*.fits.gz') + \
        glob.glob(model_dir + '/seds/*/*.fits.gz') + \
        glob.glob(model_dir + '/seds/*.fits') + \
        glob.glob(model_dir + '/seds/*/*.fits')

    # Find number of models
    n_models = len(sed_files)

    if n_models == 0:
        raise Exception("No SEDs found in %s" % model_dir)
    else:
        log.info("{0} SEDs found in {1}".format(n_models, model_dir))

    # Find out number of apertures
    n_ap = fits.getheader(sed_files[0], memmap=False)['NAP']

    # Find out apertures
    apertures = fits.open(sed_files[0], memmap=False)[2].data['APERTURE']

    # Find out number of wavelengths
    n_wav = fits.getheader(sed_files[0], memmap=False)['NWAV']

    # Find out wavelengths
    wavelengths = fits.open(sed_files[0], memmap=False)[1].data['WAVELENGTH']

    # For model grids that are very large, it is not possible to compute all
    # fluxes in one go, so we need to process in chunks in wavelength space.
    chunk_size = min(n_wav, int(np.floor(max_ram * 1024. ** 3 / (4. * 2. * n_models * n_ap))))

    if chunk_size == n_wav:
        log.info("Producing all monochromatic files in one go")
    else:
        log.info("Producing monochromatic files in chunks of {0}".format(chunk_size))

    filters = Table()
    filters['wav'] = wavelengths
    filters['filter'] = np.zeros(wavelengths.shape, dtype='S10')

    # Figure out range of wavelength indices to use
    # (wavelengths array is sorted in reverse order)
    jlo = n_wav - wavelengths[::-1].searchsorted(wav_max)
    jhi = n_wav - wavelengths[::-1].searchsorted(wav_min)
    chunk_size = min(chunk_size, jhi - jlo + 1)

    # Loop over wavelength chunks
    for jmin in range(jlo, jhi, chunk_size):

        # Find upper wavelength to compute
        jmax = min(jmin + chunk_size - 1, jhi)

        log.info('Processing wavelengths {0} to {1}'.format(jmin, jmax))

        # Set up convolved fluxes
        fluxes = [ConvolvedFluxes(model_names=np.zeros(n_models, dtype='S30'), apertures=apertures, initialize_arrays=True) for i in range(chunk_size)]

        # Loop over SEDs
        for im, sed_file in enumerate(sed_files):

            log.debug('Processing {0}'.format(os.path.basename(sed_file)))

            # Read in SED
            s = SED.read(sed_file, unit_freq=u.Hz, unit_flux=u.mJy, order='nu')

            # Convolve
            for j in range(chunk_size):

                fluxes[j].wavelength = wavelengths[j + jmin]
                fluxes[j].apertures = apertures
                fluxes[j].model_names[im] = s.name

                if n_ap == 1:
                    fluxes[j].flux[im] = s.flux[0, j + jmin]
                    fluxes[j].error[im] = s.error[0, j + jmin]
                else:
                    fluxes[j].flux[im, :] = s.flux[:, j + jmin]
                    fluxes[j].error[im, :] = s.error[:, j + jmin]

        for j in range(chunk_size):
            fluxes[j].write('{:s}/convolved/MO{:03d}.fits'.format(model_dir, j + jmin + 1),
                            overwrite=overwrite)
            filters['filter'][j + jmin] = "MO{:03d}".format(j + jmin + 1)

    return filters
