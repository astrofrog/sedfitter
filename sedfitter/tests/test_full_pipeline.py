import os
import shutil
import tempfile
import pytest

import numpy as np
from astropy.table import Table
from astropy import units as u


DATA = """
source_1 0.0 0.0 1 1 1 0.2 0.1 1.3 0.2 1.5 0.3
source_2 0.0 0.0 1 1 1 0.2 0.05 1.2 0.1 1.8 0.3
"""


def generate_random_models_1(models_dir, aperture_dependent=False):
    """
    Generate random models in original grid format
    """

    # Create fake SEDs
    np.random.seed(12345)
    os.mkdir(os.path.join(models_dir, 'seds'))

    names = []

    for i in range(5):

        from ..sed import SED

        sed = SED()

        sed.name = 'model_{0:04d}'.format(i)
        sed.distance = 1 * u.kpc

        # Make at least one SED have a different size as a regression test for a
        # bug that caused the convolution to crash if SEDs did not all have the
        # same size. Note that technically speaking the monochromatic
        # 'convolution' will be wrong in this case, but this doesn't matter
        # since the monochromatic convolution is not really going to be
        # supported going forward.
        n_wav = 100 if i < 4 else 105

        sed.wav = np.logspace(-2., 3., n_wav) * u.micron
        sed.nu = sed.wav.to(u.Hz, equivalencies=u.spectral())

        if aperture_dependent:
            sed.apertures = np.logspace(1., 6., 10) * u.au
            sed.flux = np.cumsum(np.random.random((10, n_wav)), axis=0) * u.mJy
        else:
            sed.apertures = None
            sed.flux = (1 + np.random.random((1, n_wav))) * u.mJy

        sed.error = sed.flux * np.random.random(n_wav) / 100.

        sed.write(os.path.join(models_dir, 'seds', sed.name + '_sed.fits'))

        names.append(sed.name)

    # Generate model conf file
    f = open(os.path.join(models_dir, 'models.conf'), 'w')
    f.write("name = test\n")
    f.write("length_subdir = 0\n")
    f.write("aperture_dependent = {0}\n".format('yes' if aperture_dependent else 'no'))
    f.write("logd_step = 0.02\n")
    f.close()

    # Generate model parameter file
    t = Table()
    t['MODEL_NAME'] = np.array(names, dtype='S30')
    t['par1'] = np.random.random(5)
    t['par2'] = np.random.random(5)
    t = t[[0, 4, 1, 3, 2]]
    t.write(os.path.join(models_dir, 'parameters.fits'))


def generate_random_models_2(models_dir, aperture_dependent=False):
    """
    Generate random models in original grid format
    """

    # Create fake SEDs
    np.random.seed(12345)

    from ..sed import SEDCube

    cube = SEDCube()

    cube.names = np.array(['model_{0:04d}'.format(i) for i in range(5)])
    cube.distance = 1 * u.kpc

    cube.wav = np.logspace(-2., 3., 100) * u.micron

    if aperture_dependent:
        cube.apertures = np.logspace(1., 6., 10) * u.au
        cube.val = np.cumsum(np.random.random((5, 10, 100)), axis=0) * u.mJy
    else:
        cube.apertures = None
        cube.val = (1 + np.random.random((5, 1, 100))) * u.mJy

    cube.unc = cube.val * 0.01 * np.random.random(cube.val.shape)

    cube.write(os.path.join(models_dir, 'flux.fits'))

    # Generate model conf file
    f = open(os.path.join(models_dir, 'models.conf'), 'w')
    f.write("name = test\n")
    f.write("length_subdir = 0\n")
    f.write("aperture_dependent = {0}\n".format('yes' if aperture_dependent else 'no'))
    f.write("logd_step = 0.02\n")
    f.write("version = 2\n")
    f.close()

    # Generate model parameter file
    t = Table()
    t['MODEL_NAME'] = np.array(cube.names, dtype='S')
    t['par1'] = np.random.random(5)
    t['par2'] = np.random.random(5)
    t.write(os.path.join(models_dir, 'parameters.fits'))


class BasePipelineTest(object):

    def setup_class(self):

        self.tmpdir = tempfile.mkdtemp()

        if self.model_version == 1:
            model_generator = generate_random_models_1
        else:
            model_generator = generate_random_models_2

        model_generator(self.tmpdir, aperture_dependent=self.aperture_dependent)

        from ..extinction import Extinction

        self.extinction = Extinction()
        self.extinction.wav = np.logspace(-2., 3.) * u.micron
        self.extinction.chi = self.extinction.wav.value ** -2 * u.cm ** 2 / u.g

    def teardown_class(self):
        shutil.rmtree(self.tmpdir)

    def test_broadband(self, tmpdir):

        np.random.seed(12345)

        from ..filter import Filter

        f1_wav = np.linspace(5., 1., 100) * u.micron

        f1 = Filter()
        f1.name = 'alice'
        f1.central_wavelength = 7. * u.micron
        f1.nu = f1_wav.to(u.Hz, equivalencies=u.spectral())
        f1.response = np.random.random(100)
        f1.normalize()

        f2_wav = np.linspace(15., 10., 100) * u.micron

        f2 = Filter()
        f2.name = 'bob'
        f2.central_wavelength = 12. * u.micron
        f2.nu = f2_wav.to(u.Hz, equivalencies=u.spectral())
        f2.response = np.random.random(100)
        f2.normalize()

        f3_wav = np.linspace(25., 15., 100) * u.micron

        f3 = Filter()
        f3.name = 'eve'
        f3.central_wavelength = 20. * u.micron
        f3.nu = f3_wav.to(u.Hz, equivalencies=u.spectral())
        f3.response = np.random.random(100)
        f3.normalize()

        from ..convolve import convolve_model_dir

        convolve_model_dir(self.tmpdir, filters=[f1, f2, f3])

        from ..fit import fit

        data_file = tmpdir.join('data').strpath
        open(data_file, 'w').write(DATA.strip())

        output_file = tmpdir.join('output').strpath

        fit(data_file, ['bob', 'alice', 'eve'], [1., 3., 3.] * u.arcsec, self.tmpdir, output_file,
            extinction_law=self.extinction,
            distance_range=[1., 2.] * u.kpc,
            av_range=[0., 0.1],
            output_format=('F', 3.),
            output_convolved=False)

        self._postprocess(output_file, tmpdir)

    def test_monochromatic(self, tmpdir):

        from ..convolve import convolve_model_dir_monochromatic
        if self.model_version == 1:
            convolve_model_dir_monochromatic(self.tmpdir)
            filters = ['MO001', 'MO002', 'MO020']
        else:
            with pytest.raises(ValueError) as exc:
                convolve_model_dir_monochromatic(self.tmpdir)
            assert exc.value.args[0] == "monochromatic filters are no longer used for new-style model directories"
            filters = [3.4 * u.micron, 8.0 * u.micron, 15. * u.micron]

        from ..fit import fit

        data_file = tmpdir.join('data').strpath
        open(data_file, 'w').write(DATA.strip())

        output_file = tmpdir.join('output').strpath

        fit(data_file, filters, [1., 3., 3.] * u.arcsec, self.tmpdir, output_file,
            extinction_law=self.extinction,
            distance_range=[1., 2.] * u.kpc,
            av_range=[0., 0.1],
            output_format=('F', 3.),
            output_convolved=False)

        self._postprocess(output_file, tmpdir)

    def _postprocess(self, output_file, tmpdir):

        from ..plot import plot

        plots_dir = tmpdir.join('plots_sed_auto').strpath

        plot(output_file,
             plot_mode='A',
             output_dir=plots_dir,
             select_format=('F', 3.),
             format='png')

        plots_dir = tmpdir.join('plots_sed_manual').strpath

        plot(output_file,
             plot_mode='A',
             output_dir=plots_dir,
             select_format=('F', 3.),
             format='png',
             show_convolved=False, show_sed=True,
             x_mode='M', x_range=(0.1, 2000),
             y_mode='M', y_range=(1.e-14, 2e-8),
             plot_max=100)

        from ..plot_params_1d import plot_params_1d

        plots_dir = tmpdir.join('plots_1d').strpath

        plot_params_1d(output_file, 'par1',
                       log_x=True,
                       output_dir=plots_dir,
                       select_format=('F', 2.), format='png')

        from ..plot_params_2d import plot_params_2d

        plots_dir = tmpdir.join('plots_2d').strpath

        plot_params_2d(output_file, 'par1', 'par2',
                       log_x=True, log_y=True,
                       output_dir=plots_dir,
                       select_format=('F', 2.), format='png')

        from ..extract_parameters import extract_parameters

        output_prefix = tmpdir.join('extract_parameters').strpath

        extract_parameters(output_file, output_prefix)

        additional = {}
        additional['add'] = {}
        for i in range(5):
            additional['add']['model_{0:04d}'.format(i)] = np.random.random()

        from ..write_parameters import write_parameters

        output_ascii_file = tmpdir.join('write_parameters').strpath

        write_parameters(output_file, output_ascii_file)

        write_parameters(output_file, output_ascii_file, additional=additional)

        from ..write_parameter_ranges import write_parameter_ranges

        output_ascii_file = tmpdir.join('write_parameters_ranges').strpath

        write_parameter_ranges(output_file, output_ascii_file)

        write_parameter_ranges(output_file, output_ascii_file, additional=additional)

        from ..filter_output import filter_output

        filter_output(output_file, output_good='auto', output_bad='auto', cpd=3.)


class TestApertureIndependentPipeline1(BasePipelineTest):

    aperture_dependent = False
    model_version = 1


class TestApertureIndependentPipeline2(BasePipelineTest):

    aperture_dependent = False
    model_version = 2


class TestApertureDependentPipeline1(BasePipelineTest):

    aperture_dependent = True
    model_version = 1


class TestApertureDependentPipeline2(BasePipelineTest):

    aperture_dependent = True
    model_version = 2
