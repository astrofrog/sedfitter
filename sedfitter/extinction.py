import numpy as np
from scipy.interpolate import interp1d

import cPickle as pickle

class Extinction(object):

    def __init__(self, filename=None):
        '''
        Initialize an extinction object from a file
        '''
        
        if filename:

            f = np.loadtxt(filename, dtype=[('wav', float), ('kappa', float)],
                           usecols=[0, 3])

            self._wav = f['wav']
            self._kappa = f['kappa']
        
            self._update_interp()
        
    def _update_interp(self):
        self.kappa = interp1d(self._wav, self._kappa, bounds_error=False,
                               fill_value=0)

    def av(self, wavelengths):
        '''
        Interpolate the Av at given wavelengths
        '''
        return -0.4 * self.kappa(wavelengths) / self.kappa(0.55000)
        
    def write_binary(self, file_handle):
        pickle.dump(len(self._wav), file_handle)
        file_handle.write(self._wav.astype(np.float32).tostring())
        file_handle.write(self._kappa.astype(np.float32).tostring())
        
    def read_binary(self, file_handle):
        n_wav = pickle.load(file_handle)
        self._wav = np.fromstring(file_handle.read(n_wav*4), dtype=np.float32)
        self._kappa = np.fromstring(file_handle.read(n_wav*4), dtype=np.float32)
        self._update_interp()