from __future__ import print_function, division

import os

import numpy as np

from .convolved_fluxes import ConvolvedFluxes
from . import fitting_routines as f
from .utils import parfile
from astropy.table import Table
from astropy import units as u


class Models(object):

    def __init__(self):

        self.names = None
        self.fluxes = None
        self.wavelengths = None
        self.apertures = None
        self.logd = None
        self.extended = []

    @classmethod
    def read(cls, directory, filters, distance_range=None, remove_resolved=False):

        m = cls()

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
                    m.n_distances = 1
                    m.distances = np.array([distance_range[0]])
                else:
                    m.n_distances = 1 + (np.log10(distance_range[1]) - np.log10(distance_range[0])) / modpar['logd_step']
                    m.distances = np.logspace(np.log10(distance_range[0]), np.log10(distance_range[1]), m.n_distances)
                print("   Number of distances :  %i" % m.n_distances)
            else:
                raise Exception("For aperture-dependent models, a distange range is required")
        else:
            m.n_distances = None
            m.distances = None

        print("")
        print(" ------------------------------------------------------------")
        print("  => Reading in convolved fluxes")
        print(" ------------------------------------------------------------")
        print("")

        model_fluxes = []
        m.wavelengths = []

        for filt in filters:

            filename = '%s/convolved/%s.fits' % (directory, filt['name'])

            if not os.path.exists(filename):
                if os.path.exists(filename + '.gz'):
                    filename += '.gz'
                else:
                    raise Exception("File not found: " + filename)

            print("   Reading " + filename)

            conv = ConvolvedFluxes.read(filename)

            m.wavelengths.append(conv.wavelength.to(u.micron).value)

            if m.n_distances is not None:
                apertures_au = filt['aperture_arcsec'] * m.distances * 1.e3
                conv = conv.interpolate(apertures_au)
                conv.flux = conv.flux / m.distances ** 2
                m.logd = np.log10(m.distances)
                if remove_resolved:
                    m.extended.append(apertures_au[np.newaxis, :] < conv.radius_sigma_50[:, np.newaxis])

            model_fluxes.append(conv.flux.to(u.mJy).value)

        if m.n_distances is not None:
            m.fluxes = np.column_stack(model_fluxes).reshape(conv.n_models, len(filters), m.n_distances)
            m.fluxes = m.fluxes.swapaxes(1, 2)
            if remove_resolved:
                m.extended = np.column_stack(m.extended).reshape(conv.n_models, len(filters), m.n_distances)
                m.extended = m.extended.swapaxes(1, 2)
        else:
            m.fluxes = np.column_stack(model_fluxes)

        try:
            m.names = np.char.strip(conv.model_names)
        except:
            m.names = np.array([x.strip() for x in conv.model_names], dtype=conv.model_names.dtype)

        m.n_models = conv.n_models

        m.valid = m.fluxes != 0

        m.fluxes[~m.valid] = -np.inf
        m.fluxes[m.valid] = np.log10(m.fluxes[m.valid])

        return m

    def fit(self, source, av_law, sc_law, av_min, av_max):

        weight, log_flux, log_error = source.get_log_fluxes()

        if self.fluxes.ndim == 2:  # Aperture-independent fitting

            # Use 2-parameter linear regression to find the best-fit av and scale for each model
            residual = log_flux - self.fluxes
            av_best, sc_best = f.linear_regression(residual, weight, av_law, sc_law)

            # Use optimal scaling for Avs that are outside range
            reset1 = (av_best < av_min)
            reset2 = (av_best > av_max)
            av_best[reset1] = av_min
            av_best[reset2] = av_max
            reset = reset1 | reset2
            sc_best[reset] = f.optimal_scaling(residual[reset] - av_best[reset][:, np.newaxis] * av_law[np.newaxis, :], weight, sc_law)

            # Compute best-fit model in each case
            model = av_best[:, np.newaxis] * av_law[np.newaxis, :] + sc_best[:, np.newaxis] * sc_law[np.newaxis, :]

            # Calculate the chi-squared value
            ch_best = f.chi_squared(source.valid, residual, log_error, weight, model)

            # Extract convolved model fluxes for best-fit
            model_fluxes = model + self.fluxes

        elif self.fluxes.ndim == 3:  # Aperture dependent fitting

            # Use optimal scaling to fit the Av
            residual = log_flux - self.fluxes
            av_best = f.optimal_scaling(residual, weight, av_law)

            # Reset to valid range
            av_best[av_best < av_min] = av_min
            av_best[av_best > av_max] = av_max

            # Compute best-fit model in each case
            model = av_best[:, :, np.newaxis] * av_law[np.newaxis, np.newaxis, :]

            # Calculate the chi-squared value
            ch_best = f.chi_squared(source.valid, residual, log_error, weight, model)

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


def load_parameter_table(model_dir):

    if os.path.exists(model_dir + '/parameters.fits'):
        t = Table.read(model_dir + '/parameters.fits')
    elif os.path.exists(model_dir + '/parameters.fits.gz'):
        t = Table.read(model_dir + '/parameters.fits.gz')
    else:
        raise Exception("Parameter file not found in %s" % model_dir)

    return t