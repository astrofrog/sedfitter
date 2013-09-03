from __future__ import print_function, division

try:
    import cPickle as pickle
except ImportError:
    import pickle

import numpy as np
from ..utils.validator import validate_array


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
            self._wav = validate_array('wav', value, ndim=1, shape=None if self.chi is None else self.chi.shape)

    @property
    def chi(self):
        return self._chi

    @chi.setter
    def chi(self, value):
        if value is None:
            self._chi = None
        else:
            self._chi = validate_array('chi', value, ndim=1, shape=None if self.wav is None else self.wav.shape)

    @classmethod
    def from_file(cls, filename, columns=[0, 1]):

        e = cls()

        f = np.loadtxt(filename, dtype=[('wav', float), ('chi', float)],
                       usecols=columns)

        e.wav = f['wav']
        e.chi = f['chi']

        return e

    def get_av(self, wav):
        '''
        Interpolate the Av at given wavelengths

        Parameters
        ----------
        wav : sequence
            The wavelengths at which to interpolate the visual extinction
        '''
        return -0.4 * np.interp(wav, self.wav, self.chi, left=0., right=0.) \
            / np.interp(0.55, self.wav, self.chi)

    def from_table(self, table):
        self.wav = table['wav']
        self.chi = table['chi']

    def to_table(self):
        from astropy.table import Table, Column
        from astropy import units as u
        t = Table()
        t.add_column(Column('wav', self.wav, units=u.micron))
        t.add_column(Column('chi', self.chi, units=u.cm ** 2 / u.g))
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
