import pytest
import numpy as np

from .. import Source


def test_source_init():
    s = Source()


def test_source_name():
    s = Source()
    s.name = "Source Name"


@pytest.mark.parametrize('value', [1, 1., object(), np.array([1, 2, 3])])
def test_source_name_invalid(value):
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.name = value
    assert exc.value.args[0] == 'name should be a string'


def test_source_position():
    s = Source()
    s.x = 1.
    s.y = 2.


@pytest.mark.parametrize('value', ['a', object(), np.array([1, 2, 3])])
def test_source_position_invalid(value):

    s = Source()

    with pytest.raises(TypeError) as exc:
        s.x = value
    assert exc.value.args[0] == "x should be a scalar floating point value"

    with pytest.raises(TypeError) as exc:
        s.y = value
    assert exc.value.args[0] == "y should be a scalar floating point value"


@pytest.mark.parametrize('value', [[1, 0, 1, 1, 2], (0, 1, 2, 3, 4), np.array([1., 2., 3.])])
def test_source_valid(value):
    s = Source()
    s.valid = value


@pytest.mark.parametrize('value', [1, 1., 'a', object()])
def test_source_valid_invalid(value):
    s = Source()
    with pytest.raises(TypeError) as exc:
        s.valid = value
    assert exc.value.args[0] == 'valid should be a 1-d sequence'


def test_source_valid_ints():
    s = Source()
    with pytest.raises(ValueError) as exc:
        s.valid = [1., 2.3, 4.]
    assert exc.value.args[0] == "valid values should be integers"


def test_source_valid_range():
    s = Source()
    with pytest.raises(ValueError) as exc:
        s.valid = [-1, 2, 3]
    assert exc.value.args[0] == "valid values should be in the range [0,4]"


@pytest.mark.parametrize(('attribute', 'value'), zip(['flux', 'error'], [[1., 2., 3., 4., 5], (1, 2, 3), np.array([1, 2, 3, 4])]))
def test_source_flux(attribute, value):
    s = Source()
    setattr(s, attribute, value)


@pytest.mark.parametrize(('attribute', 'value'), zip(['flux', 'error'], [1, 1., 'a', object()]))
def test_source_flux_invalid(attribute, value):
    s = Source()
    with pytest.raises(TypeError) as exc:
        setattr(s, attribute, value)
    assert exc.value.args[0] == '{0} should be a 1-d sequence'.format(attribute)


@pytest.mark.parametrize(('attribute'), ['valid', 'flux', 'error'])
def test_dimensions_match(attribute):
    s = Source()
    setattr(s, attribute, [1, 2, 3])
    for other in ['valid', 'flux', 'error']:
        if other != attribute:
            with pytest.raises(ValueError) as exc:
                setattr(s, other, [1, 2, 3, 4])
            assert exc.value.args[0] == "{0} has incorrect length (expected 3 but found 4)".format(other)


def test_dict():

    # Populate source object
    s = Source()
    s.name = "Source name"
    s.x = 1.
    s.y = 2.
    s.valid = [1, 0, 1]
    s.flux = [3., 4.4, 5.5]
    s.error = [0.1, 0.4, 0.3]
    d = s.to_dict()

    # Make sure the dictionary has the values we expect
    assert d['name'] == "Source name"
    assert d['x'] == 1.
    assert d['y'] == 2.
    assert np.all(d['valid'] == np.array([1, 0, 1]))
    assert np.all(d['flux'] == np.array([3., 4.4, 5.5]))
    assert np.all(d['error'] == np.array([0.1, 0.4, 0.3]))

    # Make sure that to_dict/from_dict round-trip
    assert Source.from_dict(d) == s
