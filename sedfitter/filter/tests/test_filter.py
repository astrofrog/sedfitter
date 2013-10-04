import os

import pytest
import numpy as np

from astropy import units as u

from .. import Filter

ROOT = os.path.dirname(__file__)


def test_init():
    Filter()


def test_wav():
    f = Filter()
    f.wav = [1,2,3] * u.m


@pytest.mark.parametrize('value', [1. * u.m, [[1,2,3],[4,5,6]] * u.m])
def test_wav_invalid_shape(value):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.wav = value
    assert exc.value.args[0] == 'wav should be a 1-d sequence'


@pytest.mark.parametrize('value', ['string', object(), 1., [1,2,3]])
def test_wav_invalid(value):
    f = Filter()
    with pytest.raises(TypeError) as exc:
        f.wav = value
    assert exc.value.args[0] == 'wavelengths should be given as a Quantity object with units of distance'


@pytest.mark.parametrize('filename', ['2J.txt', 'I1.txt', 'S3.txt'])
def test_read(filename):
    f = Filter.read(os.path.join(ROOT, 'data', filename))
    if filename == '2J.txt':
        assert len(f.wav) == 107
    elif filename == 'I1.txt':
        assert len(f.wav) == 391
    else:
        assert len(f.wav) == 21