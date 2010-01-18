import numpy as np
from scipy.interpolate import interp1d


class Extinction(object):

    def __init__(self, filename):
        '''
        Initialize an extinction object from a file
        '''

        f = np.loadtxt(filename, dtype=[('wav', float), ('kappa', float)],
                       usecols=[0, 3])

        self._kappa = interp1d(f['wav'], f['kappa'], bounds_error=False,
                               fill_value=0)

    def kappa(self, wavelengths):
        '''
        Interpolate the opacity at given wavelengths
        '''
        return self._kappa(wavelengths)

    def av(self, wavelengths):
        '''
        Interpolate the Av at given wavelengths
        '''
        return -0.4 * self.kappa(wavelengths) / self.kappa(0.55000)
