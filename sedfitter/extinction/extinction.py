from __future__ import print_function, division

import numpy as np

from astropy import units as u
from astropy.table import Table

from ..utils.validator import validate_array

__all__ = ['Extinction']


class Extinction(object):

    def __init__(self):
        self.wav = None
        self.chi = None

    @property
    def wav(self):
        return self._wav

    @wav.setter
    def wav(self, value):
        if value is None:
            self._wav = None
        else:
            self._wav = validate_array('wav', value, ndim=1,
                                       shape=None if self.chi is None else self.chi.shape,
                                       physical_type='length')

    @property
    def chi(self):
        return self._chi

    @chi.setter
    def chi(self, value):
        if value is None:
            self._chi = None
        else:
            self._chi = validate_array('chi', value, ndim=1,
                                       shape=None if self.wav is None else self.wav.shape,
                                       physical_type='area per unit mass')

    @classmethod
    def from_file(cls, filename, columns=(0, 1),
                  wav_unit=u.micron, chi_unit=u.cm ** 2 / u.g):
        """
        Read an extinction law from an ASCII file.

        This reads in two columns: the wavelength, and the opacity (in units
        of area per unit mass).

        Parameters
        ----------
        filename : str, optional
            The name of the file to read the extinction law from
        columns : tuple or list, optional
            The columns to use for the wavelength and opacity respectively
        wav_unit : :class:`~astropy.units.core.Unit`
            The units to assume for the wavelength
        chi_unit : :class:`~astropy.units.core.Unit`
            The units to assume for the opacity
        """

        self = cls()

        f = np.loadtxt(filename, dtype=[('wav', float), ('chi', float)],
                       usecols=columns)

        self.wav = f['wav'] * wav_unit
        self.chi = f['chi'] * chi_unit

        return self

    def get_av(self, wav):
        """
        Interpolate the Av at given wavelengths

        Parameters
        ----------
        wav : :class:`~astropy.units.quantity.Quantity`
            The wavelengths at which to interpolate the visual extinction.
        """
        if isinstance(wav, u.Quantity) and wav.unit.is_equivalent(u.m):
            return (-0.4 * np.interp(wav.to(self.wav.unit), self.wav, self.chi, left=0., right=0.)
                    / np.interp(([0.55] * u.micron).to(self.wav.unit), self.wav, self.chi))
        else:
            raise TypeError("wav should be given as a Quantity object with units of length")

    @classmethod
    def from_table(cls, table):
        self = cls()
        self.wav = table['wav'].data * table['wav'].unit
        self.chi = table['chi'].data * table['chi'].unit
        return self

    def to_table(self):
        t = Table()
        t['wav'] = self.wav
        t['chi'] = self.chi
        return t

    def __getstate__(self):
        return {
            'wav': self.wav,
            'chi': self.chi,
        }

    def __setstate__(self, d):
        self.__init__()
        self.wav = d['wav']
        self.chi = d['chi']
