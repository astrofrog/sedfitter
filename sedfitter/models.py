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

            model_fluxes.append(conv.fluxes)

        if self.n_distances:
            self.fluxes = np.vstack(model_fluxes).reshape(len(filters), len(distances), conv.n_models)
        else:
            self.fluxes = np.column_stack(model_fluxes)
            
        self.names = conv.model_names
        self.n_models = conv.n_models

        self.valid = self.fluxes <> 0

        self.fluxes[~self.valid] = -np.inf
        self.fluxes[self.valid] = np.log10(self.fluxes[self.valid])

        return

    def fit(self, source, av_law, sc_law):

        if self.fluxes.ndim == 2: # Aperture-independent fitting

            # Tile source info to match model shape
            # valid = np.tile(source.valid, self.n_models).reshape(self.fluxes.shape)
            # flux = np.tile(source.logflux, self.n_models).reshape(self.fluxes.shape)
            # flux_error = np.tile(source.logerror, self.n_models).reshape(self.fluxes.shape)
            # weight = np.tile(source.weight, self.n_models).reshape(self.fluxes.shape)

            # Use 2-parameter linear regression to find the best-fit av and scale for each model
            residual = source.logflux - self.fluxes
            av_best, sc_best = f.linear_regression(residual, source.weight, av_law, sc_law)

            # Compute best-fit model in each case
            model = av_best[:, np.newaxis] * av_law[np.newaxis, :] + sc_best[:, np.newaxis] * sc_law[np.newaxis, :]

            # Calculate the chi-squared value
            ch_best = f.chi_squared(source.valid, residual, source.logerror, source.weight, model)

        elif self.fluxes.ndim == 3: # Aperture dependent fitting

            # Tile source info to match model shape (n_wav, n_distances, n_models)
            valid = np.tile(source.valid, self.n_models*self.n_distances).reshape(self.fluxes.shape)
            flux = np.tile(source.logflux, self.n_models*self.n_distances).reshape(self.fluxes.shape)
            flux_error = np.tile(source.logerror, self.n_models*self.n_distances).reshape(self.fluxes.shape)
            weight = np.tile(source.weight, self.n_models*self.n_distances).reshape(self.fluxes.shape)

            print valid[:,3,3]
        
        else:
        
            raise Exception("Unexpected number of dimensions in flux array")



        return av_best, sc_best, ch_best
