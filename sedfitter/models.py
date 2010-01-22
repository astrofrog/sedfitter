import os

import numpy as np

from convolved_fluxes import ConvolvedFluxes
import fitting_routines as f

class Models(object):

    def __init__(self, *args, **kwargs):

        self.names = None
        self.fluxes = None
        self.wavelengths = None
        self.apertures = None
        self.logd = None

        if args:
            self.read(*args, **kwargs)

        return

    def read(self, directory, filters, distances=None):

        if type(distances) == np.ndarray:
            self.n_distances = len(distances)
        else:
            self.n_distances = None

        model_fluxes = []
        self.wavelengths = []

        for filt in filters:

            filename = '%s/convolved/%s.fits' % (directory, filt['name'])

            if not os.path.exists(filename):
                if os.path.exists(filename+'.gz'):
                    filename += '.gz'
                else:
                    raise Exception("File not found: "+filename)

            print "Reading "+filename

            conv = ConvolvedFluxes(filename)

            self.wavelengths.append(conv.wavelength)

            if self.n_distances:
                apertures_au = filt['aperture_arcsec'] * distances * 1.e3
                conv.interpolate(apertures_au)
                conv.fluxes = conv.fluxes / distances**2
                self.logd = np.log10(distances)

            model_fluxes.append(conv.fluxes)

        if self.n_distances:
            self.fluxes = np.column_stack(model_fluxes).reshape(conv.n_models, len(filters), len(distances))
            self.fluxes = self.fluxes.swapaxes(1,2)
            print self.fluxes.shape

        else:
            self.fluxes = np.column_stack(model_fluxes)
            print self.fluxes.shape

        self.names = conv.model_names
        self.n_models = conv.n_models

        self.valid = self.fluxes <> 0

        self.fluxes[~self.valid] = -np.inf
        self.fluxes[self.valid] = np.log10(self.fluxes[self.valid])

        return

    def fit(self, source, av_law, sc_law):

        if self.fluxes.ndim == 2: # Aperture-independent fitting

            # Use 2-parameter linear regression to find the best-fit av and scale for each model
            residual = source.logflux - self.fluxes
            av_best, sc_best = f.linear_regression(residual, source.weight, av_law, sc_law)

            # Compute best-fit model in each case
            model = av_best[:, np.newaxis] * av_law[np.newaxis, :] + sc_best[:, np.newaxis] * sc_law[np.newaxis, :]

            # Calculate the chi-squared value
            ch_best = f.chi_squared(source.valid, residual, source.logerror, source.weight, model)

        elif self.fluxes.ndim == 3: # Aperture dependent fitting

            # print "Log10(flux) = ", source.logflux
            # print "Model fluxes = ", self.fluxes[0,0,:]

            # Use optimal scaling to fit the Av
            residual = source.logflux - self.fluxes
            
            # print "Residual = ", residual[0,0,:]
            # print "Av law = ", av_law
            # print "Weights = ",source.weight
            
            av_best = f.optimal_scaling(residual, source.weight, av_law)

            # print "Av = ", av_best[0,0]

            # Compute best-fit model in each case
            model = av_best[:, :, np.newaxis] * av_law[np.newaxis, np.newaxis, :]

            # Calculate the chi-squared value
            ch_best = f.chi_squared(source.valid, residual, source.logerror, source.weight, model)

            # print "Chi^2 = ", ch_best[0,0]

            # Find best-fit distance in each case
            best = np.argmin(ch_best, axis=1)

            sc_best = self.logd[best]

            ch_best = ch_best[np.arange(self.n_models),best]
            av_best = av_best[np.arange(self.n_models),best]

        else:

            raise Exception("Unexpected number of dimensions in flux array")

        return av_best, sc_best, ch_best
