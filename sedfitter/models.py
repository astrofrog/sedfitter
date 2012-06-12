from __future__ import print_function, division

import os

import numpy as np

from .convolve import ConvolvedFluxes
from . import fitting_routines as f
from . import parfile


class Models(object):

    def __init__(self, *args, **kwargs):

        self.names = None
        self.fluxes = None
        self.wavelengths = None
        self.apertures = None
        self.logd = None
        self.extended = []

        if args:
            self.read(*args, **kwargs)

        return

    def read(self, directory, filters, distance_range=None, remove_resolved=False):

        # Read in model parameters
        modpar = parfile.read("%s/models.conf" % directory, 'conf')

        print(" ------------------------------------------------------------")
        print("  => Model parameters")
        print(" ------------------------------------------------------------")
        print("")
        print("   Models              :  %s" % modpar['name'])
        print("   Log[d] stepping     :  %g" % modpar['logd_step'])

        if modpar['aperture_dependent']:
            if distance_range:
                if distance_range[0] == distance_range[1]:
                    self.n_distances = 1
                    self.distances = np.array([distance_range[0]])
                else:
                    self.n_distances = 1 + (np.log10(distance_range[1]) - np.log10(distance_range[0])) / modpar['logd_step']
                    self.distances = np.logspace(np.log10(distance_range[0]), np.log10(distance_range[1]), self.n_distances)
                print("   Number of distances :  %i" % self.n_distances)
            else:
                raise Exception("For aperture-dependent models, a distange range is required")
        else:
            self.n_distances = None
            self.distances = None

        print("")
        print(" ------------------------------------------------------------")
        print("  => Reading in convolved fluxes")
        print(" ------------------------------------------------------------")
        print("")

        model_fluxes = []
        self.wavelengths = []

        for filt in filters:

            filename = '%s/convolved/%s.fits' % (directory, filt['name'])

            if not os.path.exists(filename):
                if os.path.exists(filename + '.gz'):
                    filename += '.gz'
                else:
                    raise Exception("File not found: " + filename)

            print("   Reading " + filename)

            conv = ConvolvedFluxes()
            conv.read(filename)

            self.wavelengths.append(conv.wavelength)

            if self.n_distances is not None:
                apertures_au = filt['aperture_arcsec'] * self.distances * 1.e3
                conv = conv.interpolate(apertures_au)
                conv.flux = conv.flux / self.distances ** 2
                self.logd = np.log10(self.distances)
                if remove_resolved:
                    self.extended.append(apertures_au[np.newaxis, :] < conv.radius_sigma_50[:, np.newaxis])

            model_fluxes.append(conv.flux)

        if self.n_distances is not None:
            self.fluxes = np.column_stack(model_fluxes).reshape(conv.n_models, len(filters), self.n_distances)
            self.fluxes = self.fluxes.swapaxes(1, 2)
            if remove_resolved:
                self.extended = np.column_stack(self.extended).reshape(conv.n_models, len(filters), self.n_distances)
                self.extended = self.extended.swapaxes(1, 2)
        else:
            self.fluxes = np.column_stack(model_fluxes)

        try:
            self.names = np.char.strip(conv.model_names)
        except:
            self.names = np.array([x.strip() for x in conv.model_names], dtype=conv.model_names.dtype)

        self.n_models = conv.n_models

        self.valid = self.fluxes != 0

        self.fluxes[~self.valid] = -np.inf
        self.fluxes[self.valid] = np.log10(self.fluxes[self.valid])

        return

    def fit(self, source, av_law, sc_law, av_min, av_max):

        if self.fluxes.ndim == 2:  # Aperture-independent fitting

            # Use 2-parameter linear regression to find the best-fit av and scale for each model
            residual = source.logflux - self.fluxes
            av_best, sc_best = f.linear_regression(residual, source.weight, av_law, sc_law)

            # Use optimal scaling for Avs that are outside range
            reset1 = (av_best < av_min)
            reset2 = (av_best > av_max)
            av_best[reset1] = av_min
            av_best[reset2] = av_max
            reset = reset1 | reset2
            sc_best[reset] = f.optimal_scaling(residual[reset] - av_best[reset][:, np.newaxis] * av_law[np.newaxis, :], source.weight, sc_law)

            # Compute best-fit model in each case
            model = av_best[:, np.newaxis] * av_law[np.newaxis, :] + sc_best[:, np.newaxis] * sc_law[np.newaxis, :]

            # Calculate the chi-squared value
            ch_best = f.chi_squared(source.valid, residual, source.logerror, source.weight, model)

            # Extract convolved model fluxes for best-fit
            model_fluxes = model + self.fluxes

        elif self.fluxes.ndim == 3:  # Aperture dependent fitting

            # Use optimal scaling to fit the Av
            residual = source.logflux - self.fluxes
            av_best = f.optimal_scaling(residual, source.weight, av_law)

            # Reset to valid range
            av_best[av_best < av_min] = av_min
            av_best[av_best > av_max] = av_max

            # Compute best-fit model in each case
            model = av_best[:, :, np.newaxis] * av_law[np.newaxis, np.newaxis, :]

            # Calculate the chi-squared value
            ch_best = f.chi_squared(source.valid, residual, source.logerror, source.weight, model)

            # Remove extended objects
            if type(self.extended) == np.ndarray:
                reset = np.any(self.extended[:, :, source.valid > 0], axis=2)
                ch_best[reset] = np.inf

            # Find best-fit distance in each case
            best = np.argmin(ch_best, axis=1)

            sc_best = self.logd[best]

            ch_best = ch_best[np.arange(self.n_models), best]
            av_best = av_best[np.arange(self.n_models), best]

            # Extract convolved model fluxes for best-fit
            model_fluxes = (model + self.fluxes)[np.arange(self.n_models), best, :]

        else:

            raise Exception("Unexpected number of dimensions in flux array")

        return av_best, sc_best, ch_best, self.names, model_fluxes
