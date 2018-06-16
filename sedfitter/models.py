from __future__ import print_function, division

import os

import numpy as np
from astropy.table import Table
from astropy import units as u

from .convolved_fluxes import ConvolvedFluxes, MonochromaticFluxes
from . import fitting_routines as f
from .utils import parfile
from .utils.validator import validate_array
from .fit_info import FitInfo
from .utils.io import read_table

__all__ = ['Models']


class Models(object):

    def __init__(self):

        self.names = None
        self.fluxes = None
        self.distances = None
        self.apertures = None
        self.logd = None
        self.wavelengths = None
        self.distances = None
        self.extended = []

    @property
    def wavelengths(self):
        """
        The wavelengths at which the models are defined
        """
        return self._wavelengths

    @wavelengths.setter
    def wavelengths(self, value):
        if value is None:
            self._wavelengths = None
        else:
            self._wavelengths = validate_array('wavelengths', value, domain='positive', ndim=1, physical_type='length')

    @property
    def distances(self):
        """
        The distances at which the models are defined
        """
        return self._distances

    @distances.setter
    def distances(self, value):
        if value is None:
            self._distances = None
        else:
            self._distances = validate_array('distances', value, domain='positive', ndim=1, physical_type='length')

    @property
    def apertures(self):
        """
        The apertures at which the fluxes are defined
        """
        return self._apertures

    @apertures.setter
    def apertures(self, value):
        if value is None:
            self._apertures = None
        else:
            self._apertures = validate_array('apertures', value, domain='positive', ndim=1, physical_type='length')

    @property
    def fluxes(self):
        """
        The model fluxes
        """
        return self._fluxes

    @fluxes.setter
    def fluxes(self, value):
        if value is None:
            self._fluxes = value
        else:
            if self.n_distances is None:
                self._fluxes = validate_array('fluxes', value, ndim=2,
                                              shape=(self.n_models, self.n_wav),
                                              physical_type=('power', 'flux', 'spectral flux density'))
            else:
                self._fluxes = validate_array('fluxes', value, ndim=3,
                                              shape=(self.n_models, self.n_distances, self.n_wav),
                                              physical_type=('power', 'flux', 'spectral flux density'))

    @property
    def n_ap(self):
        if self.apertures is None:
            return 1
        else:
            return len(self.apertures)

    @property
    def n_wav(self):
        if self.wavelengths is None:
            return None
        else:
            return len(self.wavelengths)

    @property
    def n_distances(self):
        if self.distances is None:
            return None
        else:
            return len(self.distances)

    @property
    def n_models(self):
        if self.names is None:
            return None
        else:
            return len(self.names)

    @property
    def valid(self):
        if self.fluxes is None:
            return None
        else:
            return self.fluxes != 0

    @property
    def log_fluxes_mJy(self):
        values = np.zeros(self.fluxes.shape)
        values[~self.valid] = -np.inf
        values[self.valid] = np.log10(self.fluxes[self.valid].to(u.mJy).value)
        return values

    @classmethod
    def read(cls, directory, filters, distance_range=None, remove_resolved=False):
        modpar = parfile.read("%s/models.conf" % directory, 'conf')
        if modpar.get('version', 1) == 1:
            return cls._read_version_1(directory, filters,
                                       distance_range=distance_range,
                                       remove_resolved=remove_resolved)
        else:
            return cls._read_version_2(directory, filters,
                                       distance_range=distance_range,
                                       remove_resolved=remove_resolved)

    @classmethod
    def _read_version_1(cls, directory, filters, distance_range=None, remove_resolved=None):

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

            distance_range_kpc = distance_range.to(u.kpc).value

            if distance_range:
                if distance_range_kpc[0] == distance_range_kpc[1]:
                    n_distances = 1
                    m.distances = np.array([distance_range_kpc[0]]) * u.kpc
                else:
                    n_distances = 1 + (np.log10(distance_range_kpc[1]) - np.log10(distance_range_kpc[0])) / modpar['logd_step']
                    m.distances = np.logspace(np.log10(distance_range_kpc[0]), np.log10(distance_range_kpc[1]), n_distances) * u.kpc
                print("   Number of distances :  %i" % m.n_distances)
            else:
                raise Exception("For aperture-dependent models, a distange range is required")

        print("")
        print(" ------------------------------------------------------------")
        print("  => Reading in convolved fluxes")
        print(" ------------------------------------------------------------")
        print("")

        m.wavelengths = np.zeros(len(filters)) * u.micron

        for ifilt, filt in enumerate(filters):

            filename = '%s/convolved/%s.fits' % (directory, filt['name'])

            if not os.path.exists(filename):
                if os.path.exists(filename + '.gz'):
                    filename += '.gz'
                else:
                    raise Exception("File not found: " + filename)

            print("   Reading " + filename)

            conv = ConvolvedFluxes.read(filename)

            if ifilt == 0:
                if m.n_distances is None:
                    model_fluxes = np.zeros((conv.n_models, len(filters))) * u.mJy
                    extended = None
                else:
                    model_fluxes = np.zeros((conv.n_models, m.n_distances, len(filters))) * u.mJy
                    extended = np.zeros((conv.n_models, m.n_distances, len(filters)), dtype=bool)

            m.wavelengths[ifilt] = conv.central_wavelength

            if m.n_distances is not None:
                apertures_au = filt['aperture_arcsec'] * m.distances.to(u.pc).value * u.au
                conv = conv.interpolate(apertures_au)
                conv.flux = conv.flux * (u.kpc / m.distances) ** 2
                m.logd = np.log10(m.distances.to(u.kpc).value)
                if remove_resolved:
                    extended[:, :, ifilt] = apertures_au[np.newaxis,:] < conv.find_radius_sigma(0.5)[:, np.newaxis]
                model_fluxes[:, :, ifilt] = conv.flux
            else:
                model_fluxes[:, ifilt] = conv.flux[:, 0]

        try:
            m.names = np.char.strip(conv.model_names)
        except:
            m.names = np.array([x.strip() for x in conv.model_names], dtype=conv.model_names.dtype)

        m.fluxes = model_fluxes

        if extended is not None:
            m.extended = extended

        return m

    @classmethod
    def _read_version_2(cls, directory, filters, distance_range=None, remove_resolved=None):

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

            distance_range_kpc = distance_range.to(u.kpc).value

            if distance_range:
                if distance_range_kpc[0] == distance_range_kpc[1]:
                    n_distances = 1
                    m.distances = np.array([distance_range_kpc[0]]) * u.kpc
                else:
                    n_distances = 1 + (np.log10(distance_range_kpc[1]) - np.log10(distance_range_kpc[0])) / modpar['logd_step']
                    m.distances = np.logspace(np.log10(distance_range_kpc[0]), np.log10(distance_range_kpc[1]), n_distances) * u.kpc
                print("   Number of distances :  %i" % m.n_distances)
            else:
                raise Exception("For aperture-dependent models, a distange range is required")

        print("")
        print(" ------------------------------------------------------------")
        print("  => Reading in convolved fluxes")
        print(" ------------------------------------------------------------")
        print("")

        # Start off by reading in main flux cube
        from .sed.cube import SEDCube
        cube = SEDCube.read(os.path.join(directory, 'flux.fits'))

        # Initialize model flux array and array to indicate whether models are
        # extended
        if m.n_distances is None:
            model_fluxes = np.zeros((cube.n_models, len(filters))) * u.mJy
            extended = None
        else:
            model_fluxes = np.zeros((cube.n_models, m.n_distances, len(filters))) * u.mJy
            extended = np.zeros((cube.n_models, m.n_distances, len(filters)), dtype=bool)

        # Define empty wavelength array
        m.wavelengths = np.zeros(len(filters)) * u.micron

        for ifilt, filt in enumerate(filters):

            if 'name' in filt:

                filename = '%s/convolved/%s.fits' % (directory, filt['name'])

                if not os.path.exists(filename):
                    if os.path.exists(filename + '.gz'):
                        filename += '.gz'
                    else:
                        raise Exception("File not found: " + filename)

                print("   Reading " + filename)

                conv = ConvolvedFluxes.read(filename)

                m.wavelengths[ifilt] = conv.central_wavelength

            elif 'wav' in filt:

                # Find wavelength index
                wavelength_index = np.argmin(np.abs(cube.wav - filt['wav']))

                print("   Reading fluxes at {0}".format(filt['wav']))

                conv = MonochromaticFluxes.from_sed_cube(cube, wavelength_index)

                m.wavelengths[ifilt] = filt['wav']

            if m.n_distances is not None:
                apertures_au = filt['aperture_arcsec'] * m.distances.to(u.pc).value * u.au
                conv = conv.interpolate(apertures_au)
                conv.flux = conv.flux * (u.kpc / m.distances) ** 2
                m.logd = np.log10(m.distances.to(u.kpc).value)
                # TODO: rather than compute the radius for each model, just
                # check directly the condition.
                if remove_resolved:
                    extended[:, :, ifilt] = apertures_au[np.newaxis,:] < conv.find_radius_sigma(0.5)[:, np.newaxis]
                model_fluxes[:, :, ifilt] = conv.flux
            else:
                model_fluxes[:, ifilt] = conv.flux[:, 0]

        try:
            m.names = np.char.strip(conv.model_names)
        except:
            m.names = np.array([x.strip() for x in conv.model_names], dtype=conv.model_names.dtype)

        m.fluxes = model_fluxes

        if extended is not None:
            m.extended = extended

        return m

    def fit(self, source, av_law, sc_law, av_min, av_max, output_convolved=False):

        weight, log_flux, log_error = source.get_log_fluxes()

        model_fluxes = self.log_fluxes_mJy

        if model_fluxes.ndim == 2:  # Aperture-independent fitting

            # Use 2-parameter linear regression to find the best-fit av and scale for each model
            residual = log_flux - model_fluxes
            av_best, sc_best = f.linear_regression(residual, weight, av_law, sc_law)

            # Use optimal scaling for Avs that are outside range
            reset1 = (av_best < av_min)
            reset2 = (av_best > av_max)
            av_best[reset1] = av_min
            av_best[reset2] = av_max
            reset = reset1 | reset2
            sc_best[reset] = f.optimal_scaling(residual[reset] - av_best[reset][:, np.newaxis] * av_law[np.newaxis, :], weight, sc_law)

            # Compute best-fit model in each case
            model = av_best[:, np.newaxis] * av_law[np.newaxis, :] + sc_best[:, np.newaxis] * sc_law[np.newaxis,:]

            # Calculate the chi-squared value
            ch_best = f.chi_squared(source.valid, residual, log_error, weight, model)

            # Extract convolved model fluxes for best-fit
            model_fluxes = model + model_fluxes

        elif model_fluxes.ndim == 3:  # Aperture dependent fitting

            # Use optimal scaling to fit the Av
            residual = log_flux - model_fluxes
            av_best = f.optimal_scaling(residual, weight, av_law)

            # Reset to valid range
            av_best[av_best < av_min] = av_min
            av_best[av_best > av_max] = av_max

            # Compute best-fit model in each case
            model = av_best[:, :, np.newaxis] * av_law[np.newaxis, np.newaxis,:]

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
            model_fluxes = (model + model_fluxes)[np.arange(self.n_models), best, :]

        else:

            raise Exception("Unexpected number of dimensions in flux array")

        info = FitInfo()
        info.source = source
        info.av = av_best
        info.sc = sc_best
        info.chi2 = ch_best
        info.model_name = self.names
        info.model_fluxes = model_fluxes
        info.sort()

        return info


def load_parameter_table(model_dir):

    if os.path.exists(model_dir + '/parameters.fits'):
        t = read_table(model_dir + '/parameters.fits')
    elif os.path.exists(model_dir + '/parameters.fits.gz'):
        t = read_table(model_dir + '/parameters.fits.gz')
    else:
        raise Exception("Parameter file not found in %s" % model_dir)

    return t
