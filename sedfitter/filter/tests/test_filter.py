import os

import pytest
import numpy as np

from astropy import units as u

from .. import Filter

ROOT = os.path.dirname(__file__)


def test_init():
    Filter()


def test_nu():
    f = Filter()
    f.nu = [1,2,3] * u.Hz


@pytest.mark.parametrize('value', [1. * u.Hz, [[1,2,3],[4,5,6]] * u.Hz])
def test_nu_invalid_shape(value):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.nu = value
    assert exc.value.args[0] == 'nu should be a 1-d sequence'


@pytest.mark.parametrize('value', ['string', object(), 1., [1,2,3]])
def test_nu_invalid_type(value):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.nu = value
    assert exc.value.args[0] == 'nu should be given as a Quantity object'


@pytest.mark.parametrize('unit', [u.arcsec, u.kg, u.Jy])
def test_nu_invalid_unit(unit):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.nu = [1., 2., 3.] * unit
    assert exc.value.args[0] == 'nu should be given in units of frequency'


@pytest.mark.parametrize('filename', ['2J.txt', 'I1.txt', 'S3.txt'])
def test_read(filename):
    f = Filter.read(os.path.join(ROOT, 'data', filename))
    if filename == '2J.txt':
        assert len(f.nu) == 107
    elif filename == 'I1.txt':
        assert len(f.nu) == 391
    else:
        assert len(f.nu) == 21
        
        
def test_norm():
    f = Filter()
    f.nu = [1.,2.,3.] * u.Hz
    f.response = [1.,2.,1.]
    f.normalize()
    np.testing.assert_allclose(f.response,[1/3., 2./3., 1./3.])

    
def test_norm_int():
    f = Filter()
    f.nu = [1,2,3] * u.Hz
    f.response = [1,2,1]
    f.normalize()
    np.testing.assert_allclose(f.response,[1/3., 2./3., 1./3.])
    
    
def test_norm_inverted():
    f = Filter()
    f.nu = [3.,2.,1.] * u.Hz
    f.response = [1.,2.,1.]
    f.normalize()
    np.testing.assert_allclose(f.response,[1/3., 2./3., 1./3.])
