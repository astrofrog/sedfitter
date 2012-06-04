'''
This file defines the base Filter class
'''

import os

import numpy as np

from ..utils.integrate import integrate_subset

c = 299792458


class Filter(object):

    def __init__(self, *args, **kwargs):
        self.name = None
        self.wavelength = None
        self.wav = None
        self.nu = None
        self.r = None
        if len(args) > 0:
            self.read(*args, **kwargs)

    def read(self, filename):
        '''
        Read a filter from a file

        Parameters
        ----------
        filename: str
            The name of the file containing the filter
        '''

        # Read in central wavelength
        self.wavelength = float(open(filename, 'rb').readline().split('=')[1])

        # Read in spectral response curve
        self.wav = np.loadtxt(filename, usecols=[0], dtype=float)
        self.r = np.loadtxt(filename, usecols=[1], dtype=float)

        # Compute frequency
        self.nu = c / (self.wav * 1.e-6)

        # Set name
        if self.name is None:
            self.name = os.path.basename(filename).split('.')[0]

    def rebin(self, nu_new):
        '''
        Re-bin the filter onto a new frequency grid

        Parameters
        ----------
        nu_new: np.ndarray
            The new frequency grid

        Returns
        -------
        filter: Filter
            The binned filter
        '''

        # Define new filter
        f = Filter()
        f.name = self.name
        f.wavelength = self.wavelength
        f.nu = nu_new
        f.wav = c / self.nu * 1.e6

        # Compute re-binned transmission

        f.r = np.zeros(nu_new.shape)

        for i in range(len(f.r)):

            if i == 0:
                nu1 = nu_new[0]
            else:
                nu1 = 0.5 * (nu_new[i - 1] + nu_new[i])

            if i == len(nu_new) - 1:
                nu2 = nu_new[-1]
            else:
                nu2 = 0.5 * (nu_new[i] + nu_new[i + 1])

            nu1 = min(max(nu1, self.nu[0]), self.nu[-1])
            nu2 = min(max(nu2, self.nu[0]), self.nu[-1])

            if nu2 != nu1:
                f.r[i] = integrate_subset(self.nu, self.r, nu1, nu2)

        return f
