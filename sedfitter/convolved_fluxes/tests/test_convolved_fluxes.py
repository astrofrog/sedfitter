import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal_nulp

from astropy import units as u

from .. import ConvolvedFluxes

def assert_allclose_quantity(q1, q2):
    np.testing.assert_allclose(q1.value, q2.to(q1.unit).value)


def test_init():
    ConvolvedFluxes()


def test_wavelength():
    c = ConvolvedFluxes()
    c.central_wavelength = 3.14 * u.micron


@pytest.mark.parametrize('value', ['string', object(), np.array([1, 2, 3])])
def test_wavelength_invalid_type(value):
    c = ConvolvedFluxes()
    with pytest.raises(TypeError) as exc:
        c.central_wavelength = value
    assert exc.value.args[0] == 'central_wavelength should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.Jy])
def test_wavelength_invalid_unit(unit):
    c = ConvolvedFluxes()
    with pytest.raises(TypeError) as exc:
        c.central_wavelength = [1., 2., 3.] * unit
    assert exc.value.args[0] == 'central_wavelength should be given in units of length'


def test_aperture_none():
    c = ConvolvedFluxes()
    c.apertures = None


def test_aperture_list():
    c = ConvolvedFluxes()
    c.apertures = [1., 2., 3.] * u.au


@pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
def test_aperture_invalid_type(value):

    c = ConvolvedFluxes()

    with pytest.raises(TypeError) as exc:
        c.apertures = value
    assert exc.value.args[0] == 'apertures should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.Jy])
def test_aperture_invalid_unit(unit):

    c = ConvolvedFluxes()

    with pytest.raises(TypeError) as exc:
        c.apertures = [1., 2., 3.] * unit
    assert exc.value.args[0] == 'apertures should be given in units of length'


def test_model_name_list():
    c = ConvolvedFluxes()
    c.model_names = ['a', 'b', 'c']


@pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
def test_model_name_invalid_type(value):

    c = ConvolvedFluxes()

    with pytest.raises(TypeError) as exc:
        c.model_names = value
    assert exc.value.args[0] == 'model_names should be a 1-d sequence'


def test_model_name_unset():

    c = ConvolvedFluxes()

    with pytest.raises(ValueError) as exc:
        c.flux = [1., 2., 3.] * u.mJy
    assert exc.value.args[0] == "model_names has not been set"

    with pytest.raises(ValueError) as exc:
        c.error = [1., 2., 3.] * u.mJy
    assert exc.value.args[0] == "model_names has not been set"


def test_single():
    c = ConvolvedFluxes()
    c.model_names = ['a', 'b', 'c', 'd', 'e']  # this sets the number of models
    c.apertures = None  # this sets the number of apertures
    c.flux = [(1.,), (2.,), (3.,), (4.,), (5.,)] * u.mJy
    c.error = [(0.1,), (0.2,), (0.3,), (0.4,), (0.5,)] * u.mJy


def test_single_invalid_lengths():

    c = ConvolvedFluxes()
    c.apertures = None
    c.model_names = ['a', 'b', 'c', 'd']

    with pytest.raises(ValueError) as exc:
        c.flux = [(1.,), (2.,), (3.,)] * u.mJy
    assert exc.value.args[0] == 'flux has incorrect shape (expected (4, 1) but found (3, 1))'

    with pytest.raises(ValueError) as exc:
        c.error = [(0.1,), (0.2,), (0.3,), (0.4,), (0.5,), (0.6,)] * u.mJy
    assert exc.value.args[0] == 'error has incorrect shape (expected (4, 1) but found (6, 1))'


@pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
def test_single_invalid_types(value):

    c = ConvolvedFluxes()

    c.apertures = None
    c.model_names = ['a', 'b', 'c', 'd']

    with pytest.raises(TypeError) as exc:
        c.flux = value
    assert exc.value.args[0] == 'flux should be given as a Quantity object'

    with pytest.raises(TypeError) as exc:
        c.error = value
    assert exc.value.args[0] == 'error should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.m])
def test_single_invalid_units(unit):

    c = ConvolvedFluxes()

    c.apertures = None
    c.model_names = ['a', 'b', 'c', 'd']

    with pytest.raises(TypeError) as exc:
        c.flux = [(1.,), (2.,), (3.,), (4.,)] * unit
    assert exc.value.args[0] == 'flux should be given in units of power, flux, spectral flux density'

    with pytest.raises(TypeError) as exc:
        c.error = [(1.,), (2.,), (3.,), (4.,)] * unit
    assert exc.value.args[0] == 'error should be given in units of power, flux, spectral flux density'

def test_single_io_roundtrip(tmpdir):

    c1 = ConvolvedFluxes()
    c1.model_names = ['a', 'b', 'c', 'd', 'e']
    c1.apertures = None
    c1.flux = [(1.,), (2.,), (3.,), (4.,), (5.,)] * u.mJy
    c1.error = [(0.1,), (0.2,), (0.3,), (0.4,), (0.5,)] * u.mJy

    filename = str(tmpdir.join('test_single.fits'))

    c1.write(filename)
    c2 = ConvolvedFluxes.read(filename)

    assert c1 == c2


def test_single_interpolate(tmpdir):

    c1 = ConvolvedFluxes()
    c1.model_names = ['a', 'b', 'c', 'd', 'e']
    c1.apertures = None
    c1.flux = [(1.,), (2.,), (3.,), (4.,), (5.,)] * u.mJy
    c1.error = [(0.1,), (0.2,), (0.3,), (0.4,), (0.5,)] * u.mJy

    c2 = c1.interpolate([0.5, 1.0, 1.5] * u.au)

    assert np.all(c1.model_names == c2.model_names)

    assert np.all(c2.apertures == [0.5, 1.0, 1.5] * u.au)

    assert np.all(c2.flux[:, 0] == c1.flux[:, 0])
    assert np.all(c2.flux[:, 1] == c1.flux[:, 0])
    assert np.all(c2.flux[:, 2] == c1.flux[:, 0])

    assert np.all(c2.error[:, 0] == c1.error[:, 0])
    assert np.all(c2.error[:, 1] == c1.error[:, 0])
    assert np.all(c2.error[:, 2] == c1.error[:, 0])


def test_multiple():
    c = ConvolvedFluxes()
    c.model_names = ['a', 'b', 'c']
    c.apertures = [1., 2.] * u.au
    c.flux = [[1., 2.], [3., 4.], [5., 6.]] * u.mJy
    c.error = [[0.1, 0.2], [0.1, 0.2], [0.1, 0.2]] * u.mJy


def test_multiple_invalid_lengths():

    c = ConvolvedFluxes()

    c.apertures = [1., 2.] * u.au
    c.model_names = ['a', 'b', 'c']

    with pytest.raises(ValueError) as exc:
        c.flux = np.ones((5, 2)) * u.mJy
    assert exc.value.args[0] == 'flux has incorrect shape (expected (3, 2) but found (5, 2))'

    with pytest.raises(ValueError) as exc:
        c.error = np.ones((3, 3)) * u.mJy
    assert exc.value.args[0] == 'error has incorrect shape (expected (3, 2) but found (3, 3))'


@pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
def test_multiple_invalid_types(value):

    c = ConvolvedFluxes()

    c.apertures = [1., 2.] * u.au
    c.model_names = ['a', 'b', 'c']

    with pytest.raises(TypeError) as exc:
        c.flux = value
    assert exc.value.args[0] == 'flux should be given as a Quantity object'

    with pytest.raises(TypeError) as exc:
        c.error = value
    assert exc.value.args[0] == 'error should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.m])
def test_multiple_invalid_unit(unit):

    c = ConvolvedFluxes()

    c.apertures = [1., 2.] * u.au
    c.model_names = ['a', 'b', 'c']

    with pytest.raises(TypeError) as exc:
        c.flux = [[1., 2.], [3., 4.], [5., 6.]] * unit
    assert exc.value.args[0] == 'flux should be given in units of power, flux, spectral flux density'

    with pytest.raises(TypeError) as exc:
        c.error = [[1., 2.], [3., 4.], [5., 6.]] * unit
    assert exc.value.args[0] == 'error should be given in units of power, flux, spectral flux density'


def test_multiple_roundtrip(tmpdir):

    c1 = ConvolvedFluxes()
    c1.model_names = ['a', 'b', 'c']
    c1.apertures = [1., 2.] * u.au
    c1.flux = [[1., 2.], [3., 4.], [5., 6.]] * u.mJy
    c1.error = [[0.1, 0.2], [0.1, 0.2], [0.1, 0.2]] * u.mJy

    filename = str(tmpdir.join('test_multiple.fits'))

    c1.write(filename)
    c2 = ConvolvedFluxes.read(filename)

    assert c1 == c2


def test_multiple_interpolate():

    c1 = ConvolvedFluxes()
    c1.model_names = ['a', 'b']
    c1.apertures = [1., 2., 3.] * u.au
    c1.flux = [[1., 2., 3.], [4., 5., 6.]] * u.mJy
    c1.error = [[0.1, 0.2, 0.4], [0.5, 0.3, 0.1]] * u.mJy

    c2 = c1.interpolate([1.5, 2.5] * u.au)

    assert np.all(c1.model_names == c2.model_names)

    assert_allclose_quantity(c2.apertures, [1.5, 2.5] * u.au)

    assert_allclose_quantity(c2.flux[:, 0], np.array([1.5, 4.5]) * u.mJy)
    assert_allclose_quantity(c2.flux[:, 1], np.array([2.5, 5.5]) * u.mJy)

    assert_allclose_quantity(c2.error[:, 0], np.array([0.15, 0.4]) * u.mJy)
    assert_allclose_quantity(c2.error[:, 1], np.array([0.3, 0.2]) * u.mJy)
