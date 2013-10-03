import os
import shutil
import tempfile

import numpy as np
from astropy import units as u


def generate_random_models(models_dir, aperture_dependent=False):

    # Create fake SEDs
    np.random.seed(12345)
    os.mkdir(os.path.join(models_dir, 'seds'))

    for i in range(5):

        from ..sed import SED

        sed = SED()

        sed.name = 'model_{0:04d}'.format(i)
        sed.distance = 1 * u.kpc
        sed.wav = np.logspace(-2., 3., 100) * u.micron
        sed.nu = sed.wav.to(u.Hz, equivalencies=u.spectral())

        if aperture_dependent:
            sed.apertures = np.logspace(1., 6., 10) * u.au
            sed.flux = np.cumsum(np.random.random((10, 100)), axis=0) * u.mJy
        else:
            sed.apertures = None
            sed.flux = np.random.random((1,100)) * u.mJy

        sed.error = sed.flux * np.random.random(100)

        sed.write(os.path.join(models_dir, 'seds', sed.name + '.fits'))


class BasePipelineTest(object):

    def setup_class(self):

        self.tmpdir = tempfile.mkdtemp()
        generate_random_models(self.tmpdir, aperture_dependent=self.aperture_dependent)

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_full_pipeline(self):

        from ..convolve import convolve_model_dir
        from ..filter import Filter

        np.random.seed(12345)

        f1 = Filter()
        f1.name = 'alice'
        f1.wav = np.linspace(1., 5., 100) * u.micron
        f1.nu = f1.wav.to(u.Hz, equivalencies=u.spectral())
        f1.r = np.random.random(100)

        f2 = Filter()
        f2.name = 'bob'
        f2.wav = np.linspace(10., 15., 100) * u.micron
        f2.nu = f2.wav.to(u.Hz, equivalencies=u.spectral())
        f2.r = np.random.random(100)

        convolve_model_dir(self.tmpdir, filters=[f1, f2], overwrite=True)


class TestApertureIndependentPipeline(BasePipelineTest):
    aperture_dependent = False


class TestApertureDependentPipeline(BasePipelineTest):
    aperture_dependent = True
