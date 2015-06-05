from __future__ import print_function, division

import os

import numpy as np
from astropy import units as u

from ..utils.integrate import integrate_subset, integrate
from ..utils.validator import validate_array, validate_scalar

c = 299792458.


class Filter(object):
    """
    A filter used to convolve SED models.

    Parameters
    ----------
    name : str
        The name of the filter (typically short and with no spaces, such as
        ``2K`` or 2MASS K-band or ``I1`` for IRAC 3.6 microns).
    central_wavelength : :class:`~astropy.units.quantity.Quantity`
        The central wavelength (in units of length)
    nu : :class:`~astropy.units.quantity.Quantity`
        The frequencies at which the filter is defined.
    response : `numpy.ndarray`
        The relative response at the different frequencies. Note that this
        should already take into account the assumptions about whether this is
        the response per photon or energy, and any other assumptions about
        converting the SED to a monochromatic flux. See the documentation for
        more details.
    """
    def __init__(self, name=None, central_wavelength=None, nu=None, response=None):
        self.name = name
        self.central_wavelength = central_wavelength
        self.nu = nu
        self.response = response

    @classmethod
    def read(cls, filename):
        """
        Read a filter from a file.

        This method assumes that the file contains a column with wavelengths
        (in micron) and the other column contains the response. It also
        assumes that the first line contains a header with the central
        wavelength (in microns), using the following syntax::

            # wav = 1.25

        Parameters
        ----------
        filename: str
            The name of the file containing the filter
        """

        self = cls()

        # Read in central wavelength
        self.central_wavelength = float(open(filename, 'r').readline().split('=')[1]) * u.micron

        # Read in spectral response curve
        self.wav = np.loadtxt(filename, usecols=[0], dtype=float) * u.micron
        self.response = np.loadtxt(filename, usecols=[1], dtype=float)

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
        if self.nu is None:
            raise ValueError("nu has not been set")
        if self.response is None:
            raise ValueError("response has not been set")
        self.response = self.response / np.abs(integrate(self.nu.to(u.Hz).value, self.response))

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
        f.central_wavelength = self.central_wavelength
        f.nu = nu_new

        self_nu_hz = self.nu.to(u.Hz).value
        nu_new_hz = f.nu.to(u.Hz).value

        # Compute re-binned transmission

        f.response = np.zeros(nu_new.shape)

        for i in range(len(f.response)):

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
                f.response[i] = integrate_subset(self_nu_hz, self.response, nu1, nu2)

        return f

    @property
    def central_wavelength(self):
        """
        The central or characteristic wavelength of the filter
        """
        return self._wavelength

    @central_wavelength.setter
    def central_wavelength(self, value):
        if value is None:
            self._wavelength = None
        else:
            self._wavelength = validate_scalar('central_wavelength', value,
                                               domain='strictly-positive',
                                               physical_type='length')

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
            self._nu = validate_array('nu', value, domain='strictly-positive', ndim=1,
                                      physical_type='frequency')

    @property
    def response(self):
        """
        The filter response
        """
        return self._r

    @response.setter
    def response(self, value):
        if value is None:
            self._r = None
        else:
            self._r = validate_array('r', value, domain='positive', ndim=1,
                                     shape=None if self.nu is None else (len(self.nu),))
