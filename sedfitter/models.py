import os

import numpy as np

from convolved_fluxes import ConvolvedFluxes


def linear_regression(data, weights, pattern1, pattern2):

    c1 = np.sum(data*pattern1*weights, axis=1)
    c2 = np.sum(data*pattern2*weights, axis=1)
    m11 = np.sum(pattern1*pattern1*weights, axis=1)
    m12 = np.sum(pattern1*pattern2*weights, axis=1)
    m22 = np.sum(pattern2*pattern2*weights, axis=1)

    inv_det = 1./(m11*m22-m12*m12)

    p1 = (m22*c1-m12*c2)*inv_det
    p2 = (m11*c2-m12*c1)*inv_det

    return p1, p2


def optimal_scaling(data, weights, pattern1):

    return np.sum(data*pattern1*weights) / np.sum(pattern1*pattern1*weights)


def chi_squared(valid, data, error, weight, model):

    chi2_array = np.zeros(data.shape, dtype=np.float32)

    elem = (valid==1) | (valid==4)
    chi2_array[elem] = (data[elem] - model[elem])**2 * weight[elem]

    elem = ((valid==2) & (data > model)) | ((valid==3) & (data > model))
    hard = error==1.
    chi2_array[elem & hard] = 1.e30
    chi2_array[elem & ~hard] = -2. * np.log10(1.-error[elem & ~hard])

    return np.sum(chi2_array, axis=1)


class Models(object):

    def __init__(self, *args):

        self.names = None
        self.fluxes = None
        self.wavelengths = None

        if args:
            self.read(*args)

        return

    def read(self, directory, filters):

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
            self.model_fluxes.append(conv.fluxes)

        self.fluxes = np.column_stack(self.model_fluxes)
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
        av_best, sc_best = linear_regression(residual, weight, av_law, sc_law)

        # Compute best-fit model in each case
        model = av_best[:, np.newaxis] * av_law[np.newaxis, :] + sc_best[:, np.newaxis] * sc_law[np.newaxis, :]

        # Calculate the chi-squared value
        ch_best = chi_squared(valid, residual, flux_error, weight, model)

        return av_best, sc_best, ch_best
