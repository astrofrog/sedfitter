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


@pytest.mark.parametrize('value', [[1, 0, 1, 1, 2, 9], (0, 1, 2, 3, 4), np.array([1., 2., 3.])])
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
    assert exc.value.args[0] == "valid values should be in the range [0:4] or set to 9"


@pytest.mark.parametrize(('attribute', 'value'), list(zip(['flux', 'error'], [[1., 2., 3., 4., 5], (1, 2, 3), np.array([1, 2, 3, 4])])))
def test_source_flux(attribute, value):
    s = Source()
    setattr(s, attribute, value)


@pytest.mark.parametrize(('attribute', 'value'), list(zip(['flux', 'error'], [1, 1., 'a', object()])))
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


def test_get_log_fluxes():

    # Populate source object
    s = Source()
    s.name = "Source name"
    s.x = 1.
    s.y = 2.
    s.valid = [0, 1, 2, 3, 4]
    s.flux = [10., 10., 10., 10., 10.]
    s.error = [1., 1., 1., 1., 1.]

    weight, log_flux, log_error = s.get_log_fluxes()

    assert log_flux[0] == 0.
    assert log_flux[1] == 1. - 5.e-3 / np.log(10.)
    assert log_flux[2] == 1.
    assert log_flux[3] == 1.
    assert log_flux[4] == 10.

    assert log_error[0] == 0.
    assert log_error[1] == 0.1 / np.log(10.)
    assert log_error[2] == 1.
    assert log_error[3] == 1.
    assert log_error[4] == 1.

    assert weight[0] == 0.
    assert weight[1] == 1. / log_error[1] ** 2
    assert weight[2] == 0.
    assert weight[3] == 0.
    assert weight[4] == 1. / log_error[4] ** 2


def test_source_from_ascii():

    line = "a 1. 2. 0 1 1 2 0 1. 0.1 2. 0.2 3. 0.3 4. 0.4 5. 0.5"

    s = Source.from_ascii(line)

    assert s.name == "a"
    assert s.x == 1.
    assert s.y == 2.
    assert np.all(s.valid == np.array([0, 1, 1, 2, 0]))
    assert np.all(s.flux == np.array([1., 2., 3., 4., 5.]))
    assert np.all(s.error == np.array([0.1, 0.2, 0.3, 0.4, 0.5]))


def test_source_ascii_roundtrip():

    s = Source()
    s.name = "a"
    s.x = 1.
    s.y = 2.
    s.valid = [0, 1, 2, 3, 4]
    s.flux = [10., 20., 30., 40., 50.]
    s.error = [1., 0.4, 1.2, 1.3, 1.1]

    # Check that round-tripping works
    s2 = Source.from_ascii(s.to_ascii())
    assert s == s2
