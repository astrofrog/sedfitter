from __future__ import print_function, division

import numpy as np

from .. import six


def is_numpy_array(variable):
    return type(variable) in [np.ndarray,
                              np.core.records.recarray,
                              np.ma.core.MaskedArray]


def validate_1d_array(name, value):

    if type(value) in [list, tuple]:
        value = np.array(value)
    if not is_numpy_array(value) or value.ndim != 1:
        raise ValueError(name + " should be a 1-D sequence")

    return value


class Source(object):

    def __init__(self):

        self.name = ""
        self.x = 0.
        self.y = 0.
        self.valid = None
        self.flux = None
        self.error = None

    def __getstate__(self):
        return {
            'name': self.name,
            'x': self.x,
            'y': self.y,
            'valid': self.valid,
            'flux': self.flux,
            'error': self.error
        }

    def __setstate__(self, d):
        self.__init__()
        self.name = d['name']
        self.x = d['x']
        self.y = d['y']
        self.valid = d['valid']
        self.flux = d['flux']
        self.error = d['error']

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value is None:
            self._name = value
        elif isinstance(value, six.string_types):
            self._name = value
        else:
            raise TypeError("name should be a string")

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        if value is None:
            self._x = value
        elif np.isscalar(value) and np.isreal(value):
            self._x = value
        else:
            raise TypeError("x should be a scalar floating point value")

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        if value is None:
            self._y = value
        elif np.isscalar(value) and np.isreal(value):
            self._y = value
        else:
            raise TypeError("y should be a scalar floating point value")

    @property
    def valid(self):
        return self._valid

    @valid.setter
    def valid(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._valid = value
        elif isinstance(value, np.ndarray) and value.ndim == 1:
            if self.n_wav is not None and len(value) != self.n_wav:
                raise ValueError("valid has incorrect length (expected {0} but found {1})".format(self.n_wav, len(value)))
            else:
                if np.any(value.astype(int) != value):
                    raise ValueError("valid values should be integers")
                elif np.any((value < 0) | ((value > 4) & (value != 9))):
                    raise ValueError("valid values should be in the range [0:4] or set to 9")
                else:
                    self._valid = value
        else:
            raise TypeError("valid should be a 1-d sequence")

    @property
    def flux(self):
        return self._flux

    @flux.setter
    def flux(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._flux = value
        elif isinstance(value, np.ndarray) and value.ndim == 1:
            if self.n_wav is not None and len(value) != self.n_wav:
                raise ValueError("flux has incorrect length (expected {0} but found {1})".format(self.n_wav, len(value)))
            else:
                self._flux = value
        else:
            raise TypeError("flux should be a 1-d sequence")

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._error = value
        elif isinstance(value, np.ndarray) and value.ndim == 1:
            if self.n_wav is not None and len(value) != self.n_wav:
                raise ValueError("error has incorrect length (expected {0} but found {1})".format(self.n_wav, len(value)))
            else:
                self._error = value
        else:
            raise TypeError("error should be a 1-d sequence")

    @property
    def n_wav(self):
        if self.valid is not None:
            return len(self.valid)
        elif self.flux is not None:
            return len(self.flux)
        elif self.error is not None:
            return len(self.error)
        else:
            return None

    @property
    def n_data(self):
        return np.sum((self.valid == 1) | (self.valid == 4))

    def get_log_fluxes(self):

        # Initialize arrays
        log_flux = np.zeros(self.flux.shape, dtype=np.float64)
        log_error = np.zeros(self.error.shape, dtype=np.float64)
        weight = np.zeros(self.valid.shape, dtype=np.float64)

        # Fluxes
        r = self.valid == 1
        log_flux[r] = np.log10(self.flux[r]) - 0.5 * (self.error[r] / self.flux[r]) ** 2. / np.log(10.)
        log_error[r] = np.abs(self.error[r] / self.flux[r]) / np.log(10.)
        weight[r] = 1. / log_error[r] ** 2.

        # Lower and upper limits
        r = (self.valid == 2) | (self.valid == 3)
        log_flux[r] = np.log10(self.flux[r])
        log_error[r] = self.error[r]

        # Log10[Fluxes]
        r = self.valid == 4
        log_flux[r] = self.flux[r]
        log_error[r] = self.error[r]
        weight[r] = 1. / log_error[r] ** 2.

        # Ignored points
        r = self.valid == 9
        log_flux[r] = np.log10(self.flux[r]) - 0.5 * (self.error[r] / self.flux[r]) ** 2. / np.log(10.)
        log_error[r] = np.abs(self.error[r] / self.flux[r]) / np.log(10.)

        return weight, log_flux, log_error

    def __str__(self):

        weight, log_flux, log_error = self.get_log_fluxes()

        string = "Source name : %s\n" % self.name
        string += "x           : %9.5f\n" % self.x
        string += "y           : %9.5f\n" % self.y
        for j in range(self.n_wav):
            string += "F = %12.4e + / - %12.4e mJy (%1i)  Log[F] = %8.5f+ / -%8.5f\n" % \
                      (self.flux[j], self.error[j], self.valid[j], log_flux[j], log_error[j])

        return string

    @classmethod
    def from_ascii(cls, line):

        s = cls()

        cols = line.split()

        if len(cols) < 3:
            raise EOFError()

        s.name = cols[0]
        s.x = np.float64(cols[1])
        s.y = np.float64(cols[2])
        n_wav = np.int32((len(cols) - 3) / 3)
        s.valid = np.array(cols[3:3 + n_wav], dtype=int)
        flux_and_error = np.array(cols[3 + n_wav:], dtype=float)
        s.flux = flux_and_error[::2]
        s.error = flux_and_error[1::2]

        return s

    def to_ascii(self):
        line = ""
        line += "{0:30s} {1:9.5f} {2:9.5f} ".format(self.name, self.x, self.y)
        for v in self.valid:
            line += "{0:1d} ".format(v)
        for j in range(self.n_wav):
            line += "{0:11.3e} {1:11.3e} ".format(self.flux[j], self.error[j])
        return line

    @classmethod
    def from_dict(cls, source_dict):
        s = cls()
        s.name = source_dict['name']
        s.x = source_dict['x']
        s.y = source_dict['y']
        s.valid = source_dict['valid']
        s.flux = source_dict['flux']
        s.error = source_dict['error']
        return s

    def to_dict(self):
        return {
            'name': self.name,
            'x': self.x,
            'y': self.y,
            'valid': self.valid,
            'flux': self.flux,
            'error': self.error
        }

    def __eq__(self, other):
        return self.name == other.name and \
            self.x == other.x and \
            self.y == other.y and \
            np.all(self.valid == other.valid) and \
            np.all(self.flux == other.flux) and \
            np.all(self.error == other.error)
