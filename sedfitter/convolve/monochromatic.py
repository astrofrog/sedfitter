from __future__ import print_function, division

import os
import glob

import numpy as np
from astropy.logger import log
from astropy.table import Table
from astropy import units as u
from astropy.utils.console import ProgressBar

from ..convolved_fluxes import ConvolvedFluxes
from ..sed import SED
from ..models import load_parameter_table
from .. import six
from ..utils import parfile

__all__ = ['convolve_model_dir_monochromatic', ]


def convolve_model_dir_monochromatic(model_dir, overwrite=False, max_ram=8,
                                     wav_min=-np.inf * u.micron, wav_max=np.inf * u.micron):
    """
    Convolve all the model SEDs in a model directory

    Parameters
    ----------
    model_dir : str
        The path to the model directory
    overwrite : bool, optional
        Whether to overwrite the output files
    max_ram : float, optional
        The maximum amount of RAM that can be used (in Gb)
    wav_min : float, optional
        The minimum wavelength to consider. Only wavelengths above this value
        will be output.
    wav_max : float, optional
        The maximum wavelength to consider. Only wavelengths below this value
        will be output.
    """

    modpar = parfile.read(os.path.join(model_dir, 'models.conf'), 'conf')
    if modpar.get('version', 1) > 1:
        raise ValueError("monochromatic filters are no longer used for new-style model directories")

    # Create 'convolved' sub-directory if needed
    if not os.path.exists(model_dir + '/convolved'):
        os.mkdir(model_dir + '/convolved')

    # Find all SED files to convolve
    sed_files = sorted(glob.glob(model_dir + '/seds/*.fits.gz') +
                       glob.glob(model_dir + '/seds/*/*.fits.gz') +
                       glob.glob(model_dir + '/seds/*.fits') +
                       glob.glob(model_dir + '/seds/*/*.fits'))

    par_table = load_parameter_table(model_dir)

    # Find number of models
    n_models = len(sed_files)

    if n_models == 0:
        raise Exception("No SEDs found in %s" % model_dir)
    else:
        log.info("{0} SEDs found in {1}".format(n_models, model_dir))

    # Find out apertures and wavelengths
    first_sed = SED.read(sed_files[0])
    n_ap = first_sed.n_ap
    apertures = first_sed.apertures
    n_wav = first_sed.n_wav
    wavelengths = first_sed.wav

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
    jlo = n_wav - 1 - (wavelengths[::-1].searchsorted(wav_max) - 1)
    jhi = n_wav - 1 - wavelengths[::-1].searchsorted(wav_min)
    chunk_size = min(chunk_size, jhi - jlo + 1)

    # Loop over wavelength chunks
    for jmin in range(jlo, jhi, chunk_size):

        # Find upper wavelength to compute
        jmax = min(jmin + chunk_size - 1, jhi)

        log.info('Processing wavelengths {0} to {1}'.format(jmin, jmax))

        # Set up convolved fluxes
        fluxes = [ConvolvedFluxes(model_names=np.zeros(n_models, dtype='U30' if six.PY3 else 'S30'), apertures=apertures, initialize_arrays=True) for i in range(chunk_size)]

        b = ProgressBar(len(sed_files))

        # Loop over SEDs
        for im, sed_file in enumerate(sed_files):

            b.update()

            log.debug('Processing {0}'.format(os.path.basename(sed_file)))

            # Read in SED
            s = SED.read(sed_file, unit_freq=u.Hz, unit_flux=u.mJy, order='nu')

            # Convolve
            for j in range(chunk_size):

                fluxes[j].central_wavelength = wavelengths[j + jmin]
                fluxes[j].apertures = apertures
                fluxes[j].model_names[im] = s.name

                if n_ap == 1:
                    fluxes[j].flux[im] = s.flux[0, j + jmin]
                    fluxes[j].error[im] = s.error[0, j + jmin]
                else:
                    fluxes[j].flux[im, :] = s.flux[:, j + jmin]
                    fluxes[j].error[im, :] = s.error[:, j + jmin]

        for j in range(chunk_size):
            fluxes[j].sort_to_match(par_table['MODEL_NAME'])
            fluxes[j].write('{0:s}/convolved/MO{1:03d}.fits'.format(model_dir, j + jmin + 1),
                            overwrite=overwrite)
            filters['filter'][j + jmin] = "MO{0:03d}".format(j + jmin + 1)

    return filters
