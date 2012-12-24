from __future__ import print_function, division

import numpy as np


class Extinction(object):

    def __init__(self):
        self.wav = None
        self.chi = None

    @property
    def wav(self):
        return self._wav

    @wav.setter
    def wav(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._wav = value
        elif isinstance(value, np.ndarray) and value.ndim == 1:
            if self.chi is not None and len(value) != len(self.chi):
                raise ValueError("wav has incorrect length (expected {0} but found {1})".format(len(self.chi), len(value)))
            else:
                self._wav = value
        else:
            raise TypeError("wav should be a 1-d sequence")

    @property
    def chi(self):
        return self._chi

    @chi.setter
    def chi(self, value):
        if type(value) in [list, tuple]:
            value = np.array(value)
        if value is None:
            self._chi = value
        elif isinstance(value, np.ndarray) and value.ndim == 1:
            if self.wav is not None and len(value) != len(self.wav):
                raise ValueError("chi has incorrect length (expected {0} but found {1})".format(len(self.wav), len(value)))
            else:
                self._chi = value
        else:
            raise TypeError("chi should be a 1-d sequence")

    @classmethod
    def from_file(self, filename, columns=[0, 1]):

        f = np.loadtxt(filename, dtype=[('wav', float), ('chi', float)],
                       usecols=columns)

        self.wav = f['wav']
        self.chi = f['chi']

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
