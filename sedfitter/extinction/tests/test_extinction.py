import pytest

import pytest
import numpy as np

from .. import Extinction


def test_extinction_init():
    e = Extinction()

@pytest.mark.parametrize('value', [[1, 0, 1, 1, 2], (0, 1, 2, 3, 4), np.array([1., 2., 3.])])
def test_extinction_wav(value):
    e = Extinction()
    e.wav = value


@pytest.mark.parametrize('value', [1, 1., 'a', object()])
def test_extinction_wav_invalid(value):
    e = Extinction()
    with pytest.raises(TypeError) as exc:
        e.wav = value
    assert exc.value.args[0] == 'wav should be a 1-d sequence'

@pytest.mark.parametrize('value', [[1, 0, 1, 1, 2], (0, 1, 2, 3, 4), np.array([1., 2., 3.])])
def test_extinction_chi(value):
    e = Extinction()
    e.chi = value


@pytest.mark.parametrize('value', [1, 1., 'a', object()])
def test_extinction_chi_invalid(value):
    e = Extinction()
    with pytest.raises(TypeError) as exc:
        e.chi = value
    assert exc.value.args[0] == 'chi should be a 1-d sequence'


@pytest.mark.parametrize(('attribute'), ['wav', 'chi'])
def test_dimensions_match(attribute):
    e = Extinction()
    setattr(e, attribute, [1, 2, 3])
    for other in ['wav', 'chi']:
        if other != attribute:
            with pytest.raises(ValueError) as exc:
                setattr(e, other, [1, 2, 3, 4])
            assert exc.value.args[0] == "{0} has incorrect length (expected 3 but found 4)".format(other)

