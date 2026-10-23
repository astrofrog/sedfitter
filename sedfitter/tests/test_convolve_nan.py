import os

import numpy as np
from astropy import units as u
from astropy.table import Table
from numpy.testing import assert_allclose

from ..convolve import convolve_model_dir
from ..convolved_fluxes import ConvolvedFluxes
from ..filter import Filter
from ..sed import SEDCube


def test_convolve_nan(tmpdir):

    # Models may contain NaN values (for example at long wavelengths in the
    # Richardson et al. 2024 models). NaN values outside a filter should be
    # ignored, while models with NaN values inside the filter should get zero
    # flux so that they are treated as invalid when fitting.

    models_dir = tmpdir.strpath

    cube = SEDCube()
    cube.names = np.array(['model_{0:04d}'.format(i) for i in range(3)])
    cube.distance = 1 * u.kpc
    cube.wav = np.logspace(-2., 3., 100) * u.micron
    cube.apertures = None
    cube.val = np.ones((3, 1, 100)) * u.mJy
    cube.unc = cube.val * 0.1

    # Model 1 has NaN values outside the filter, model 2 inside the filter
    long_wav = cube.wav > 100 * u.micron
    in_filter = (cube.wav > 10 * u.micron) & (cube.wav < 15 * u.micron)
    cube.val[1, 0, long_wav] = np.nan
    cube.val[2, 0, in_filter] = np.nan
    cube.unc[2, 0, in_filter] = np.nan

    cube.write(os.path.join(models_dir, 'flux.fits'))

    with open(os.path.join(models_dir, 'models.conf'), 'w') as f:
        f.write("name = test\nlength_subdir = 0\naperture_dependent = no\n"
                "logd_step = 0.02\nversion = 2\n")

    t = Table()
    t['MODEL_NAME'] = np.array(cube.names, dtype='S')
    t.write(os.path.join(models_dir, 'parameters.fits'))

    filt = Filter()
    filt.name = 'bob'
    filt.central_wavelength = 12. * u.micron
    filt.nu = (np.linspace(15., 10., 100) * u.micron).to(u.Hz, equivalencies=u.spectral())
    filt.response = np.ones(100)
    filt.normalize()

    convolve_model_dir(models_dir, filters=[filt])

    conv = ConvolvedFluxes.read(os.path.join(models_dir, 'convolved', 'bob.fits'))

    assert conv.flux[0, 0] > 0
    assert_allclose(conv.flux[1, 0], conv.flux[0, 0])
    assert_allclose(conv.error[1, 0], conv.error[0, 0])
    assert conv.flux[2, 0] == 0
    assert conv.error[2, 0] == 0
