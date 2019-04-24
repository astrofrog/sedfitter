from __future__ import print_function, division

import os
import glob

import numpy as np

from astropy.logger import log
from astropy import units as u
from astropy.utils.console import ProgressBar

from ..convolved_fluxes import ConvolvedFluxes
from ..sed import SED, SEDCube
from ..models import load_parameter_table
from .. import six
from ..utils import parfile


def convolve_model_dir(model_dir, filters, overwrite=False, memmap=True):
    """
    Convolve all the model SEDs in a model directory

    Parameters
    ----------
    model_dir : str
        The path to the model directory
    filters : list
        A list of :class:`~sedfitter.filter.Filter` objects to use for the
        convolution
    overwrite : bool, optional
        Whether to overwrite the output files
    memmap : bool, optional
        Whether to use memory mapping when using the SED cubes. If you have
        enough memory, the convolution will be much faster to set this to
        False, since the whole cube will need to be read in, so it's faster to
        do it in one go than many small reads. This option is ignored if
        not using SED cubes.
    """
    modpar = parfile.read(os.path.join(model_dir, 'models.conf'), 'conf')
    if modpar.get('version', 1) == 1:
        return _convolve_model_dir_1(model_dir, filters, overwrite=overwrite)
    else:
        return _convolve_model_dir_2(model_dir, filters, overwrite=overwrite, memmap=memmap)


def _convolve_model_dir_1(model_dir, filters, overwrite=False):

    for f in filters:
        if f.name is None:
            raise Exception("filter name needs to be set")
        if f.central_wavelength is None:
            raise Exception("filter central wavelength needs to be set")

    # Create 'convolved' sub-directory if needed
    if not os.path.exists(model_dir + '/convolved'):
        os.mkdir(model_dir + '/convolved')

    # Find all SED files to convolve
    sed_files = sorted(glob.glob(model_dir + '/seds/*.fits.gz') +
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
        except (ValueError, AssertionError):
            log.info('Rebinning filters')
            binned_filters = [f.rebin(s.nu) for f in filters]
            binned_nu = s.nu

        b.update()

        # Convolve
        for i, f in enumerate(binned_filters):

            fluxes[i].central_wavelength = f.central_wavelength
            fluxes[i].apertures = apertures
            fluxes[i].model_names[im] = s.name

            if n_ap == 1:
                fluxes[i].flux[im] = np.sum(s.flux * f.response)
                fluxes[i].error[im] = np.sqrt(np.sum((s.error * f.response) ** 2))
            else:
                fluxes[i].flux[im, :] = np.sum(s.flux * f.response, axis=1)
                fluxes[i].error[im] = np.sqrt(np.sum((s.error * f.response) ** 2, axis=1))

    for i, f in enumerate(binned_filters):
        fluxes[i].sort_to_match(par_table['MODEL_NAME'])
        fluxes[i].write(model_dir + '/convolved/' + f.name + '.fits',
                        overwrite=overwrite)


def _convolve_model_dir_2(model_dir, filters, overwrite=False, memmap=True):

    for f in filters:
        if f.name is None:
            raise Exception("filter name needs to be set")
        if f.central_wavelength is None:
            raise Exception("filter central wavelength needs to be set")

    # Create 'convolved' sub-directory if needed
    if not os.path.exists(model_dir + '/convolved'):
        os.mkdir(model_dir + '/convolved')

    # Find all SED files to convolve
    sed_cube = SEDCube.read(os.path.join(model_dir, 'flux.fits'), order='nu',
                            memmap=memmap)

    par_table = load_parameter_table(model_dir)

    if not np.all(par_table['MODEL_NAME'] == sed_cube.names):
        raise ValueError("Model names in SED cube and parameter file do not match")

    log.info("{0} SEDs found in {1}".format(sed_cube.n_models, model_dir))

    # Set up convolved fluxes
    fluxes = [ConvolvedFluxes(model_names=sed_cube.names,
                              apertures=sed_cube.apertures,
                              initialize_arrays=True) for i in range(len(filters))]

    # Set up list of binned filters
    binned_filters = [f.rebin(sed_cube.nu) for f in filters]

    # We do the unit conversion - if needed - at the last minute
    val_factor = sed_cube.val.unit.to(u.mJy)
    unc_factor = sed_cube.unc.unit.to(u.mJy)

    # Loop over apertures
    for i_ap in ProgressBar(range(sed_cube.n_ap)):

        sed_val = sed_cube.val[:, i_ap, :]
        sed_unc = sed_cube.val[:, i_ap, :]

        for i, f in enumerate(binned_filters):

            response = f.response.astype(sed_val.dtype)

            fluxes[i].flux[:, i_ap] = np.sum(sed_val * response, axis=1) * val_factor
            fluxes[i].error[:, i_ap] = np.sqrt(np.sum((sed_unc * response) ** 2, axis=1)) * unc_factor

    for i, f in enumerate(binned_filters):

        fluxes[i].central_wavelength = f.central_wavelength
        fluxes[i].apertures = sed_cube.apertures
        fluxes[i].model_names = sed_cube.names

        fluxes[i].write(model_dir + '/convolved/' + f.name + '.fits',
                        overwrite=overwrite)
