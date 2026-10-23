import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table
from numpy.testing import assert_allclose

from ..parameter_columns import (aperture_to_au, interpolate_to_aperture,
                                 parameter_names, parameter_values,
                                 select_parameter, validate_aperture)

APERTURES = np.logspace(2., 4., 3) * u.au


def make_table():
    t = Table()
    t['MODEL_NAME'] = ['a', 'b']
    t['scalar'] = [1., 2.]
    t['per_ap'] = [[1., 2., 3.], [4., 5., 6.]]
    t['other'] = [[1., 2.], [3., 4.]]
    return t


def test_names_split():
    assert parameter_names(make_table(), APERTURES) == ['scalar', 'per_ap_100au', 'per_ap_1000au',
                                                        'per_ap_10000au', 'other_0', 'other_1']


def test_names_interpolate():
    assert parameter_names(make_table(), APERTURES, interpolate=True) == ['scalar', 'per_ap',
                                                                          'other_0', 'other_1']


def test_names_no_apertures():
    # Without known model apertures, vector columns are split by index
    assert parameter_names(make_table(), None) == ['scalar', 'per_ap_0', 'per_ap_1', 'per_ap_2',
                                                   'other_0', 'other_1']


def test_values_split():
    values = parameter_values(make_table(), APERTURES)
    assert_allclose(values[0], [1., 2.])
    assert_allclose(values[1], [1., 4.])
    assert_allclose(values[3], [3., 6.])
    assert_allclose(values[5], [2., 4.])


def test_values_interpolate():
    values = parameter_values(make_table(), APERTURES, aperture_au=[100., 1000.] * u.au)
    assert len(values) == 4
    assert_allclose(values[1], [1., 5.])


def test_interpolate_log_spacing():
    # Interpolation is linear in log(aperture), so the geometric mean of two
    # apertures is halfway between them
    values = np.array([[1., 2., 3.]])
    assert_allclose(interpolate_to_aperture(values, APERTURES, np.sqrt(10.) * 100. * u.au), [1.5])


def test_interpolate_bounds():
    values = np.array([[1., 2., 3.]])
    assert_allclose(interpolate_to_aperture(values, APERTURES, 1. * u.pc), [3.])
    assert_allclose(interpolate_to_aperture(values, APERTURES, 100. * u.au), [1.])
    with pytest.raises(ValueError, match='smaller than the smallest'):
        interpolate_to_aperture(values, APERTURES, 50. * u.au)


def test_interpolate_ignores_nan_neighbour():
    values = np.array([[1., 2., np.nan]])
    assert_allclose(interpolate_to_aperture(values, APERTURES, 1000. * u.au), [2.])


def test_interpolate_masked():
    values = np.ma.array([[1., 2., 3.]], mask=[[False, True, False]])
    assert np.isnan(interpolate_to_aperture(values, APERTURES, 1000. * u.au)[0])


def test_select_parameter():
    t = make_table()
    assert_allclose(select_parameter(t, 'scalar', APERTURES), [1., 2.])
    assert_allclose(select_parameter(t, 'per_ap', APERTURES, 1000. * u.au), [2., 5.])
    with pytest.raises(ValueError, match='an aperture needs to be specified'):
        select_parameter(t, 'per_ap', APERTURES)
    with pytest.raises(ValueError, match='cannot be reduced'):
        select_parameter(t, 'other', APERTURES, 1000. * u.au)


def test_aperture_to_au():
    log_distance_kpc = np.array([0., 1.])
    assert_allclose(aperture_to_au(2. * u.arcsec, log_distance_kpc).to_value(u.au), [2000., 20000.])
    assert_allclose(aperture_to_au(1. * u.pc, log_distance_kpc).to_value(u.au), [206264.806] * 2)


def test_validate_aperture():
    assert validate_aperture(None) is None
    validate_aperture(1. * u.arcsec)
    with pytest.raises(TypeError, match='length or angle'):
        validate_aperture(3.)
    with pytest.raises(TypeError, match='length$'):
        validate_aperture(1. * u.arcsec, allow_angle=False)
    with pytest.raises(TypeError):
        validate_aperture([1., 2.] * u.au)
