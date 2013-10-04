from __future__ import print_function, division

import os

import numpy as np
from astropy import units as u

from ..utils.integrate import integrate_subset, integrate
from ..utils.validator import validate_array, validate_scalar

c = 299792458.


class Filter(object):

    def __init__(self):
        self.name = None
        self.wavelength = None
        self.wav = None
        self.nu = None
        self.r = None

    @classmethod
    def read(cls, filename):
        """
        Read a filter from a file

        Parameters
        ----------
        filename: str
            The name of the file containing the filter
        """

        self = cls()

        # Read in central wavelength
        self.wavelength = float(open(filename, 'r').readline().split('=')[1]) * u.micron

        # Read in spectral response curve
        self.wav = np.loadtxt(filename, usecols=[0], dtype=float) * u.micron
        self.r = np.loadtxt(filename, usecols=[1], dtype=float)

        # Compute frequency
        self.nu = self.wav.to(u.Hz, equivalencies=u.spectral())

        # Set name
        if self.name is None:
            self.name = os.path.basename(filename).split('.')[0]

        return self

    def normalize(self):
        """
        Normalize so the integral over nu is 1
        """
        self.r /= integrate(self.nu.to(u.Hz).value, self.r)

    def rebin(self, nu_new):
        """
        Re-bin the filter onto a new frequency grid

        Parameters
        ----------
        nu_new: np.ndarray
            The new frequency grid

        Returns
        -------
        filter: Filter
            The binned filter
        """

        # Define new filter
        f = Filter()
        f.name = self.name
        f.wavelength = self.wavelength
        f.nu = nu_new
        f.wav = f.nu.to(u.micron, equivalencies=u.spectral())

        self_nu_hz = self.nu.to(u.Hz).value
        nu_new_hz = f.nu.to(u.Hz).value

        # Compute re-binned transmission

        f.r = np.zeros(nu_new.shape)

        for i in range(len(f.r)):

            if i == 0:
                nu1 = nu_new_hz[0]
            else:
                nu1 = 0.5 * (nu_new_hz[i - 1] + nu_new_hz[i])

            if i == len(nu_new_hz) - 1:
                nu2 = nu_new_hz[-1]
            else:
                nu2 = 0.5 * (nu_new_hz[i] + nu_new_hz[i + 1])

            nu1 = min(max(nu1, self_nu_hz[0]), self_nu_hz[-1])
            nu2 = min(max(nu2, self_nu_hz[0]), self_nu_hz[-1])

            if nu2 != nu1:
                f.r[i] = integrate_subset(self_nu_hz, self.r, nu1, nu2)

        return f

    @property
    def wavelength(self):
        """
        The central or characteristic wavelength of the filter
        """
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        if value is None:
            self._wavelength = None
        else:
            if isinstance(value, u.Quantity) and value.unit.is_equivalent(u.m):
                if not value.isscalar:
                    raise TypeError("wavelength should be a scalar Quantity")
                if not value > 0 * u.micron:
                    raise ValueError("wavelength should be strictly positive")
                self._wavelength = value
            else:
                raise TypeError("central wavelength should be given as a Quantity object with units of distance")

    @property
    def wav(self):
        """
        The wavelengths at which the filter is defined
        """
        return self._wav

    @wav.setter
    def wav(self, value):
        if value is None:
            self._wav = None
        else:
            if isinstance(value, u.Quantity) and value.unit.is_equivalent(u.m):
                self._wav = validate_array('wav', value, domain='strictly-positive', ndim=1, shape=None if self.nu is None else (len(self.nu),))
            else:
                raise TypeError("wavelengths should be given as a Quantity object with units of distance")

    @property
    def nu(self):
        """
        The frequencies at which the filter is defined
        """
        return self._nu

    @nu.setter
    def nu(self, value):
        if value is None:
            self._nu = None
        else:
            if isinstance(value, u.Quantity) and value.unit.is_equivalent(u.Hz):
                self._nu = validate_array('nu', value, domain='strictly-positive', ndim=1, shape=None if self.wav is None else (len(self.wav),))
            else:
                raise TypeError("frequencies should be given as a Quantity object with units of frequency")

    @property
    def r(self):
        """
        The filter response
        """
        return self._r

    @r.setter
    def r(self, value):
        if value is None:
            self._r = None
        else:
            self._r = validate_array('r', value, domain='positive', ndim=1, shape=None if (self.wav or self.nu) is None else (len(self.wav or self.nu),))
