import os

import pytest
import numpy as np

from astropy import units as u

from .. import SED

ROOT = os.path.dirname(__file__)


def test_init():
    SED()


def test_wav():
    s = SED()
    s.wav = [1,2,3] * u.m


@pytest.mark.parametrize('value', [1. * u.m, [[1,2,3],[4,5,6]] * u.m])
def test_wav_invalid_shape(value):
    s = SED()
    with pytest.raises(TypeError) as exc:
        s.wav = value
    assert exc.value.args[0] == 'wav should be a 1-d sequence'


@pytest.mark.parametrize('value', ['string', object(), 1., [1,2,3]])
def test_wav_invalid_type(value):
    s = SED()
    with pytest.raises(TypeError) as exc:
        s.wav = value
    assert exc.value.args[0] == 'wav should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.Jy])
def test_wav_invalid_unit(unit):
    s = SED()
    with pytest.raises(TypeError) as exc:
        s.wav = [1., 2., 3.] * unit
    assert exc.value.args[0] == 'wav should be given in units of length'


# def test_aperture_none():
#     s = SED()
#     s.apertures = None
#
#
# def test_aperture_list():
#     s = SED()
#     s.apertures = [1., 2., 3.]
#
#
# @pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
# def test_aperture_invalid_type(value):
#
#     s = SED()
#
#     with pytest.raises(TypeError) as exc:
#         s.apertures = value
#     assert exc.value.args[0] == 'apertures should be a 1-d sequence'
#
#
# def test_model_name_list():
#     s = SED()
#     s.model_names = ['a', 'b', 'c']
#
#
# @pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
# def test_model_name_invalid_type(value):
#
#     s = SED()
#
#     with pytest.raises(TypeError) as exc:
#         s.model_names = value
#     assert exc.value.args[0] == 'model_names should be a 1-d sequence'
#
#
# def test_model_name_unset():
#
#     s = SED()
#
#     with pytest.raises(ValueError) as exc:
#         s.flux = [1., 2., 3.]
#     assert exc.value.args[0] == "model_names has not been set"
#
#     with pytest.raises(ValueError) as exc:
#         s.error = [1., 2., 3.]
#     assert exc.value.args[0] == "model_names has not been set"
#
#
# def test_single():
#     s = SED()
#     s.model_names = ['a', 'b', 'c', 'd', 'e']  # this sets the number of models
#     s.apertures = None  # this sets the number of apertures
#     s.flux = [1., 2., 3., 4., 5.]
#     s.error = [0.1, 0.2, 0.3, 0.4, 0.5]
#
#
# def test_single_invalid_lengths():
#
#     s = SED()
#     s.apertures = None
#     s.model_names = ['a', 'b', 'c', 'd']
#
#     with pytest.raises(ValueError) as exc:
#         s.flux = [1., 2., 3.]
#     assert exc.value.args[0] == 'flux has incorrect length (expected 4 but found 3)'
#
#     with pytest.raises(ValueError) as exc:
#         s.error = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
#     assert exc.value.args[0] == 'error has incorrect length (expected 4 but found 6)'
#
#
# @pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
# def test_single_invalid_types(value):
#
#     s = SED()
#
#     s.apertures = None
#     s.model_names = ['a', 'b', 'c', 'd']
#
#     with pytest.raises(TypeError) as exc:
#         s.flux = value
#     assert exc.value.args[0] == 'flux should be a 1-d sequence'
#
#     with pytest.raises(TypeError) as exc:
#         s.error = value
#     assert exc.value.args[0] == 'error should be a 1-d sequence'
#
#
# def test_single_io_roundtrip(tmpdir):
#
#     c1 = SED()
#     c1.model_names = ['a', 'b', 'c', 'd', 'e']
#     c1.apertures = None
#     c1.flux = [1., 2., 3., 4., 5.]
#     c1.error = [0.1, 0.2, 0.3, 0.4, 0.5]
#
#     filename = str(tmpdir.join('test_single.fits'))
#
#     c1.write(filename)
#     c2 = SED.read(filename)
#
#     assert c1 == c2
#
#
# def test_single_interpolate(tmpdir):
#
#     c1 = SED()
#     c1.model_names = ['a', 'b', 'c', 'd', 'e']
#     c1.apertures = None
#     c1.flux = [1., 2., 3., 4., 5.]
#     c1.error = [0.1, 0.2, 0.3, 0.4, 0.5]
#
#     c2 = c1.interpolate([0.5, 1.0, 1.5])
#
#     assert np.all(c1.model_names == c2.model_names)
#
#     assert np.all(c2.apertures == [0.5, 1.0, 1.5])
#
#     assert np.all(c2.flux[:, 0] == c1.flux[:])
#     assert np.all(c2.flux[:, 1] == c1.flux[:])
#     assert np.all(c2.flux[:, 2] == c1.flux[:])
#
#     assert np.all(c2.error[:, 0] == c1.error[:])
#     assert np.all(c2.error[:, 1] == c1.error[:])
#     assert np.all(c2.error[:, 2] == c1.error[:])
#
#
# def test_multiple():
#     s = SED()
#     s.model_names = ['a', 'b', 'c']
#     s.apertures = [1., 2.]
#     s.flux = [[1., 2.], [3., 4.], [5., 6.]]
#     s.error = [[0.1, 0.2], [0.1, 0.2], [0.1, 0.2]]
#
#
# def test_multiple_invalid_lengths():
#
#     s = SED()
#
#     s.apertures = [1., 2.]
#     s.model_names = ['a', 'b', 'c']
#
#     with pytest.raises(ValueError) as exc:
#         s.flux = np.ones((5, 2))
#     assert exc.value.args[0] == 'flux has incorrect shape (expected (3, 2) but found (5, 2))'
#
#     with pytest.raises(ValueError) as exc:
#         s.error = np.ones((3, 3))
#     assert exc.value.args[0] == 'error has incorrect shape (expected (3, 2) but found (3, 3))'
#
#
# @pytest.mark.parametrize('value', ['string', 1, 0.5, np.zeros((2, 2, 3))])
# def test_multiple_invalid_types(value):
#
#     s = SED()
#
#     s.apertures = [1., 2.]
#     s.model_names = ['a', 'b', 'c']
#
#     with pytest.raises(TypeError) as exc:
#         s.flux = value
#     assert exc.value.args[0] == 'flux should be a 2-d array'
#
#     with pytest.raises(TypeError) as exc:
#         s.error = value
#     assert exc.value.args[0] == 'error should be a 2-d array'
#
#
# def test_multiple_roundtrip(tmpdir):
#
#     c1 = SED()
#     c1.model_names = ['a', 'b', 'c']
#     c1.apertures = [1., 2.]
#     c1.flux = [[1., 2.], [3., 4.], [5., 6.]]
#     c1.error = [[0.1, 0.2], [0.1, 0.2], [0.1, 0.2]]
#
#     filename = str(tmpdir.join('test_multiple.fits'))
#
#     c1.write(filename)
#     c2 = SED.read(filename)
#
#     assert c1 == c2
#
#
# def test_multiple_interpolate():
#
#     c1 = SED()
#     c1.model_names = ['a', 'b']
#     c1.apertures = [1., 2., 3.]
#     c1.flux = [[1., 2., 3.], [4., 5., 6.]]
#     c1.error = [[0.1, 0.2, 0.4], [0.5, 0.3, 0.1]]
#
#     c2 = c1.interpolate([1.5, 2.5])
#
#     assert np.all(c1.model_names == c2.model_names)
#
#     assert_array_almost_equal_nulp(c2.apertures, [1.5, 2.5], 10)
#
#     assert_array_almost_equal_nulp(c2.flux[:, 0], np.array([1.5, 4.5]), 10)
#     assert_array_almost_equal_nulp(c2.flux[:, 1], np.array([2.5, 5.5]), 10)
#
#     assert_array_almost_equal_nulp(c2.error[:, 0], np.array([0.15, 0.4]), 10)
#     assert_array_almost_equal_nulp(c2.error[:, 1], np.array([0.3, 0.2]), 10)




def test_read_single_aperture():
    # This actually tests a lot of internal type/dimension checking for the
    # properties.
    s = SED.read(os.path.join(ROOT,'data','kt09250g+4.0z-0.5_sed.fits.gz'))
    assert s.n_wav == 1321
    assert s.n_ap == 1


def test_read_multiple_aperture():
    # This actually tests a lot of internal type/dimension checking for the
    # properties.
    s = SED.read(os.path.join(ROOT,'data','3000030_1_sed.fits.gz'))
    assert s.n_wav == 250
    assert s.n_ap == 50


def test_roundtrip_single_aperture(tmpdir):
    s = SED.read(os.path.join(ROOT,'data','kt09250g+4.0z-0.5_sed.fits.gz'))
    temp_file = tmpdir.join('test_roundtrip_single_aperture').strpath
    s.write(temp_file)
    s2 = SED.read(temp_file)
    assert s == s2


def test_roundtrip_multiple_aperture(tmpdir):
    s = SED.read(os.path.join(ROOT,'data','3000030_1_sed.fits.gz'))
    temp_file = tmpdir.join('test_roundtrip_multipleaperture').strpath
    s.write(temp_file)
    s2 = SED.read(temp_file)
    assert s == s2

def test_full_roundtrip(tmpdir):

    n_ap = 10
    n_wav = 30

    s = SED()

    s.name = 'test'
    s.distance = 1 * u.kpc
    s.apertures = np.linspace(10, 100, n_ap) * u.au
    s.wav = np.linspace(0.01, 5000, n_wav)[::-1] * u.micron
    s.flux = np.random.random((n_ap, n_wav))[::-1] * u.mJy
    s.error = np.random.random((n_ap, n_wav)) * u.mJy

    temp_file = tmpdir.join('test_roundtrip_multipleaperture').strpath

    s.write(temp_file)
    s2 = SED.read(temp_file, unit_flux=u.mJy)

    assert s == s2
