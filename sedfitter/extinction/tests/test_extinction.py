import pytest
import numpy as np
from numpy.testing import assert_allclose

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u

from .. import Extinction



def test_extinction_init():
    Extinction()


@pytest.mark.parametrize('value', [[1, 0, 1, 1, 2], (0, 1, 2, 3, 4), np.array([1., 2., 3.])])
def test_extinction_wav(value):
    e = Extinction()
    e.wav = value * u.micron


@pytest.mark.parametrize('value', [1, 1., 'a', object()])
def test_extinction_wav_invalid_type(value):
    e = Extinction()
    with pytest.raises(TypeError) as exc:
        e.wav = value
    assert exc.value.args[0] == 'wav should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.Jy])
def test_extinction_wav_invalid_unit(unit):
    e = Extinction()
    with pytest.raises(TypeError) as exc:
        e.wav = [1., 2., 3.] * unit
    assert exc.value.args[0] == 'wav should be given in units of length'


@pytest.mark.parametrize('value', [np.array([[1, 2]]), np.array(3)])
def test_extinction_wav_invalid_shape(value):
    e = Extinction()
    with pytest.raises(TypeError) as exc:
        e.wav = value * u.micron
    assert exc.value.args[0] == 'wav should be a 1-d sequence'


@pytest.mark.parametrize('value', [[1, 0, 1, 1, 2], (0, 1, 2, 3, 4), np.array([1., 2., 3.])])
def test_extinction_chi(value):
    e = Extinction()
    e.chi = value * u.cm ** 2 / u.g


@pytest.mark.parametrize('value', [1, 1., 'a', object()])
def test_extinction_chi_invalid_type(value):
    e = Extinction()
    with pytest.raises(TypeError) as exc:
        e.chi = value
    assert exc.value.args[0] == 'chi should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.Jy])
def test_extinction_chi_invalid_unit(unit):
    e = Extinction()
    with pytest.raises(TypeError) as exc:
        e.chi = [1., 2., 3.] * unit
    assert exc.value.args[0] == 'chi should be given in units of area per unit mass'


@pytest.mark.parametrize('value', [np.array([[1, 2]]), np.array(3)])
def test_extinction_chi_invalid_shape(value):
    e = Extinction()
    with pytest.raises(TypeError) as exc:
        e.chi = value * u.cm ** 2 / u.g
    assert exc.value.args[0] == 'chi should be a 1-d sequence'


@pytest.mark.parametrize(('attribute'), ['wav', 'chi'])
def test_dimensions_match(attribute):
    units = {}
    units['wav'] = u.micron
    units['chi'] = u.cm ** 2 / u.g
    e = Extinction()
    setattr(e, attribute, [1, 2, 3] * units[attribute])
    for other in ['wav', 'chi']:
        if other != attribute:
            with pytest.raises(ValueError) as exc:
                setattr(e, other, [1, 2, 3, 4] * units[other])
            assert exc.value.args[0] == "{0} has incorrect length (expected 3 but found 4)".format(other)


def test_extinction_table():

    e = Extinction()
    e.wav = [0., 1., 2., 3.3, 5.5] * u.micron
    e.chi = [1., 3., 2., 3., 5.] * u.cm ** 2 / u.g

    table = e.to_table()

    e2 = Extinction.from_table(table)

    assert_quantity_allclose(e.wav, e2.wav)
    assert_quantity_allclose(e.chi, e2.chi)


def test_extinction_numpy_io(tmpdir):

    wav = [0., 1., 2., 3.3, 5.5]
    chi = [1., 3., 2., 3., 5.]

    filename = str(tmpdir.join('test_extinction.fits'))
    np.savetxt(filename, list(zip(wav, chi)))

    e2 = Extinction.from_file(filename)

    assert_allclose(wav, e2.wav.to(u.micron).value)
    assert_allclose(chi, e2.chi.to(u.cm**2 / u.g).value)
