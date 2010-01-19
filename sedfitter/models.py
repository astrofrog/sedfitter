import os

import numpy as np

from convolved_fluxes import ConvolvedFluxes
from import fitting routines as f

class Models(object):

    def __init__(self, *args):

        self.names = None
        self.fluxes = None
        self.wavelengths = None

        if args:
            self.read(*args)

        return

    def read(self, directory, filters, distances=None):

        self.model_fluxes = []
        self.wavelengths = []

        for filt in filters:

            filename = '%s/convolved/%s.fits' % (directory, filt)

            if not os.path.exists(filename):
                if os.path.exists(filename+'.gz'):
                    filename += '.gz'
                else:
                    raise Exception("File not found: "+filename)

            print "Reading "+filename

            conv = ConvolvedFluxes(filename)

            self.wavelengths.append(conv.wavelength)
            
            if distances:
                apertures_au = filt.ap * distances
                conv.interpolate(apertures_au)
                model_fluxes.append(conv.interp_fluxes)
            else:
                model_fluxes.append(conv.fluxes)

        self.fluxes = np.column_stack(model_fluxes)

        self.names = conv.model_names
        self.n_models = conv.n_models

        self.valid = self.fluxes <> 0

        self.fluxes[~self.valid] = -np.inf
        self.fluxes[self.valid] = np.log10(self.fluxes[self.valid])

        return

    def fit(self, source, av_law, sc_law):

        # Tile source info to match model shape
        valid = np.tile(source.valid, self.n_models).reshape(self.fluxes.shape)
        flux = np.tile(source.logflux, self.n_models).reshape(self.fluxes.shape)
        flux_error = np.tile(source.logerror, self.n_models).reshape(self.fluxes.shape)
        weight = np.tile(source.weight, self.n_models).reshape(self.fluxes.shape)

        # Use 2-parameter linear regression to find the best-fit av and scale for each model
        residual = flux - self.fluxes
        av_best, sc_best = f.linear_regression(residual, weight, av_law, sc_law)

        # Compute best-fit model in each case
        model = av_best[:, np.newaxis] * av_law[np.newaxis, :] + sc_best[:, np.newaxis] * sc_law[np.newaxis, :]

        # Calculate the chi-squared value
        ch_best = f.chi_squared(valid, residual, flux_error, weight, model)

        return av_best, sc_best, ch_best
