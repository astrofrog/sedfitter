import pytest
import numpy as np

from .. import ConvolvedFluxes


def test_fluxes_init():
    ConvolvedFluxes()


def test_fluxes_wavelength():
    c = ConvolvedFluxes()
    c.wavelength = 3.14


@pytest.mark.parametrize('value', ['string', object(), np.array([1, 2, 3])])
def test_fluxes_wavelength_invalid(value):
    c = ConvolvedFluxes()
    with pytest.raises(TypeError) as exc:
        c.wavelength = value
    assert exc.value.args[0] == 'wavelength should be a scalar floating point value'


def test_fluxes_single_aperture():
    c = ConvolvedFluxes(n_models=5)
    c.model_names = ['a', 'b', 'c', 'd', 'e']
    c.apertures = [1.]
    c.flux = [1., 2., 3., 4., 5.]
    c.error = [0.1, 0.2, 0.3, 0.4, 0.5]


def test_fluxes_single_aperture_invalid_lengths():

    c = ConvolvedFluxes(n_models=5)

    with pytest.raises(ValueError) as exc:
        c.apertures = [1., 2.]
    assert exc.value.args[0] == 'Expected 1 apertures, but got 2'

    with pytest.raises(ValueError) as exc:
        c.model_names = ['a', 'b', 'c', 'd']
    assert exc.value.args[0] == 'Expected 5 model names, but got 4'

    with pytest.raises(ValueError) as exc:
        c.flux = [1., 2., 3.]
    assert exc.value.args[0] == 'Expected 5 model fluxes, but got 3'

    with pytest.raises(ValueError) as exc:
        c.error = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    assert exc.value.args[0] == 'Expected 5 model flux errors, but got 6'


@pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
def test_fluxes_single_aperture_invalid_types(value):

    c = ConvolvedFluxes(n_models=5)

    with pytest.raises(ValueError) as exc:
        c.apertures = value
    assert exc.value.args[0] == 'apertures should be a 1-d sequence'

    with pytest.raises(ValueError) as exc:
        c.model_names = value
    assert exc.value.args[0] == 'model_names should be a 1-d sequence'

    with pytest.raises(ValueError) as exc:
        c.flux = value
    assert exc.value.args[0] == 'flux should be a 1-d sequence'

    with pytest.raises(ValueError) as exc:
        c.error = value
    assert exc.value.args[0] == 'error should be a 1-d sequence'


def test_fluxes_single_aperture_io_roundtrip(tmpdir):

    c1 = ConvolvedFluxes(n_models=5)
    c1.model_names = ['a', 'b', 'c', 'd', 'e']
    c1.apertures = [1.]
    c1.flux = [1., 2., 3., 4., 5.]
    c1.error = [0.1, 0.2, 0.3, 0.4, 0.5]

    filename = str(tmpdir.join('test_single.fits'))

    c1.write(filename)
    c2 = ConvolvedFluxes.read(filename)

    assert c1 == c2


def test_fluxes_single_aperture_interpolate(tmpdir):

    c1 = ConvolvedFluxes(n_models=5)
    c1.model_names = ['a', 'b', 'c', 'd', 'e']
    c1.apertures = [1.]
    c1.flux = [1., 2., 3., 4., 5.]
    c1.error = [0.1, 0.2, 0.3, 0.4, 0.5]

    c2 = c1.interpolate([0.5, 1.0, 1.5])

    assert np.all(c1.model_names == c2.model_names)

    assert np.all(c2.apertures == [0.5, 1.0, 1.5])

    assert np.all(c2.flux[:, 0] == c1.flux[:])
    assert np.all(c2.flux[:, 1] == c1.flux[:])
    assert np.all(c2.flux[:, 2] == c1.flux[:])

    assert np.all(c2.error[:, 0] == c1.error[:])
    assert np.all(c2.error[:, 1] == c1.error[:])
    assert np.all(c2.error[:, 2] == c1.error[:])


def test_fluxes_multiple_apertures():
    c = ConvolvedFluxes(n_models=3, n_ap=2)
    c.model_names = ['a', 'b', 'c']
    c.apertures = [1., 2.]
    c.flux = [[1., 2.], [3., 4.], [5., 6.]]
    c.error = [[0.1, 0.2], [0.1, 0.2], [0.1, 0.2]]


def test_fluxes_multiple_apertures_invalid_lengths():

    c = ConvolvedFluxes(n_models=3, n_ap=2)

    with pytest.raises(ValueError) as exc:
        c.apertures = [1.]
    assert exc.value.args[0] == 'Expected 2 apertures, but got 1'

    with pytest.raises(ValueError) as exc:
        c.model_names = ['a', 'b']
    assert exc.value.args[0] == 'Expected 3 model names, but got 2'

    with pytest.raises(ValueError) as exc:
        c.flux = np.ones((2, 2, 2))
    assert exc.value.args[0] == 'Expected (3, 2) model fluxes, but got (2, 2, 2)'

    with pytest.raises(ValueError) as exc:
        c.error = np.ones((3, 3))
    assert exc.value.args[0] == 'Expected (3, 2) model flux errors, but got (3, 3)'


@pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
def test_fluxes_multiple_apertures_invalid_types(value):

    c = ConvolvedFluxes(n_models=3, n_ap=2)

    with pytest.raises(ValueError) as exc:
        c.apertures = value
    assert exc.value.args[0] == 'apertures should be a 1-d sequence'

    with pytest.raises(ValueError) as exc:
        c.model_names = value
    assert exc.value.args[0] == 'model_names should be a 2-d array'

    with pytest.raises(ValueError) as exc:
        c.flux = value
    assert exc.value.args[0] == 'flux should be a 2-d array'

    with pytest.raises(ValueError) as exc:
        c.error = value
    assert exc.value.args[0] == 'error should be a 2-d array'


def test_fluxes_multiple_apertures_roundtrip():

    c1 = ConvolvedFluxes(n_models=3, n_ap=2)
    c1.model_names = ['a', 'b', 'c']
    c1.apertures = [1., 2.]
    c1.flux = [[1., 2.], [3., 4.], [5., 6.]]
    c1.error = [[0.1, 0.2], [0.1, 0.2], [0.1, 0.2]]

    filename = str(tmpdir.join('test_multiple.fits'))

    c1.write(filename)
    c2 = ConvolvedFluxes.read(filename)

    assert c1 == c2


def test_fluxes_multiple_apertures_interpolate():

    c1 = ConvolvedFluxes(n_models=3, n_ap=2)
    c1.model_names = ['a', 'b', 'c']
    c1.apertures = [1., 2., 3.]
    c1.flux = [[1., 2., 3.], [4., 5., 6.]]
    c1.error = [[0.1, 0.2, 0.4], [0.5, 0.3, 0.1]]

    c2 = c1.interpolate([1.5, 2.5])

    assert np.all(c1.model_names == c2.model_names)

    assert np.all(c2.apertures == [1.5, 2.5])

    assert np.all(c2.flux[:, 0] == np.array([1.5, 4.5]))
    assert np.all(c2.flux[:, 1] == np.array([2.5, 5.5]))

    assert np.all(c2.error[:, 0] == np.array([0.15, 0.4]))
    assert np.all(c2.error[:, 1] == np.array([0.3, 0.2]))
