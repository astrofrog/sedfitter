from scipy.interpolate import interp1d
import atpy
import pyfits
import numpy as np

from ..logger import log

from hyperion.util.interpolate import interp1d_fast


class ConvolvedFluxes(object):

    def __init__(self, n_models=None, n_ap=1, wavelength=None):

        self.wavelength = wavelength
        self.n_models = n_models
        self.n_ap = n_ap

        if n_models is not None:

            self.model_names = np.zeros(n_models, dtype='S30')
            self.radius_sigma_50 = np.zeros(n_models, dtype=float)
            self.radius_cumul_99 = np.zeros(n_models, dtype=float)
            self._apertures = np.zeros(n_ap)

            if n_ap == 1:
                self._flux = np.zeros(n_models, dtype=float)
                self._flux_errors = np.zeros(n_models, dtype=float)
            else:
                self._flux = np.zeros((n_models, n_ap), dtype=float)
                self._flux_errors = np.zeros((n_models, n_ap), dtype=float)

            if self._flux.ndim > 1:
                self._flux_interp = interp1d(self._apertures, self._flux[:])

            # Make the default fluxes and apertures public - this is then
            # overwritten when the user interpolates to custom apertures.
            self.flux = self._flux
            self.flux_err = self._flux_errors
            self.apertures = self._apertures

        else:

            self.flux = None
            self.flux_err = None
            self.apertures = None

    def read(self, filename):
        '''
        Read convolved flux from a FITS file
        '''

        # Open the convolved flux FITS file
        hdulist = pyfits.open(filename)

        # Try and read in the wavelength of the filter
        if 'FILTWAV' in hdulist[0].header:
            self.wavelength = hdulist[0].header['FILTWAV']
        else:
            self.wavelength = None

        # Read in number of models and apertures
        self.n_models = hdulist[0].header['NMODELS']
        self.n_ap = hdulist[0].header['NAP']

        # Read in model names
        self.model_names = hdulist[1].data.field('MODEL_NAME')

        # Read in flux and flux errors
        self._flux = hdulist[1].data.field('TOTAL_FLUX')
        self._flux_errors = hdulist[1].data.field('TOTAL_FLUX_ERR')

        # Read in 99% cumulative and 50% surface brightness radii
        self.radius_sigma_50 = hdulist[1].data.field('RADIUS_SIGMA_50')
        self.radius_cumul_99 = hdulist[1].data.field('RADIUS_CUMUL_99')

        # Read in apertures
        self._apertures = hdulist[2].data.field('APERTURE')

        # Open the convolved flux FITS file
        # t = atpy.TableSet(filename, verbose=False)

        # Try and read in the wavelength of the filter
        # if 'FILTWAV' in t.keywords:
        # if 'FILTWAV' in hdulist[0].header:
        #     self.wavelength = t.keywords['FILTWAV']
        # else:
        #     self.wavelength = None

        # Read in number of models and apertures
        # self.n_models = t.keywords['NMODELS']
        # self.n_ap = t.keywords['NAP']

        # Read in model names
        # self.model_names = t[0].MODEL_NAME

        # Read in flux and flux errors
        # self._flux = t[0].TOTAL_FLUX
        # self._flux_errors = t[0].TOTAL_FLUX_ERR

        # Read in 99% cumulative and 50% surface brightness radii
        # self.radius_sigma_50 = t[0].RADIUS_SIGMA_50
        # self.radius_cumul_99 = t[0].RADIUS_CUMUL_99

        # Read in apertures
        # self._apertures = t[1].APERTURE

        # Create an interpolating function for the flux vs aperture
        if self._flux.ndim > 1:
            self._flux_interp = interp1d(self._apertures, self._flux[:])

        # Make the default fluxes and apertures public - this is then
        # overwritten when the user interpolates to custom apertures.
        self.flux = self._flux
        self.flux_err = self._flux_errors
        self.apertures = self._apertures

        return

    def write(self, filename):
        '''
        Write convolved flux to a FITS file
        '''

        ts = atpy.TableSet()

        ts.add_keyword('FILTWAV', self.wavelength)
        ts.add_keyword('NMODELS', self.n_models)
        ts.add_keyword('NAP', self.n_ap)

        ts.append(atpy.Table(name='CONVOLVED flux'))
        ts[0].add_column('MODEL_NAME', self.model_names)
        ts[0].add_column('TOTAL_FLUX', self.flux)
        ts[0].add_column('TOTAL_FLUX_ERR', self.flux_err)
        ts[0].add_column('RADIUS_SIGMA_50', self.radius_sigma_50)
        ts[0].add_column('RADIUS_CUMUL_99', self.radius_cumul_99)

        ts.append(atpy.Table(name='APERTURES'))
        ts[1].add_column("APERTURE", self.apertures)

        ts.write(filename, verbose=False)

    def interpolate(self, apertures):
        '''
        Interpolate the flux to the apertures specified
        '''

        # If any apertures are larger than the defined max, reset to max
        apertures[apertures > self._apertures.max()] = self._apertures.max()

        # If any apertures are smaller than the defined min, raise Exception
        if np.any(apertures < self._apertures.min()):
            raise Exception("Aperture(s) requested too small")

        # Interpolate to requested apertures
        if self._flux.ndim > 1:
            self.flux = self._flux_interp(apertures)
        else:
            self.flux = np.repeat(self._flux, len(apertures)).reshape(len(self._flux), len(apertures))

        # Save the apertures
        self.apertures = apertures

    def find_radius_cumul(self, fraction):

        log.info("Finding radius containing %g%s of the flux" % (fraction * 100., '%'))

        radius = np.zeros(self.n_models)

        if self.apertures is None:

            return radius

        else:

            required = fraction * self.flux[:, -1]

            # Linear interpolation - need to loop over apertures for vectorization
            for ia in range(len(self.apertures) - 1):
                calc = (required >= self.flux[:, ia]) & (required < self.flux[:, ia + 1])
                radius[calc] = (required[calc] - self.flux[calc, ia]) / \
                               (self.flux[calc, ia + 1] - self.flux[calc, ia]) * \
                               (self.apertures[ia + 1] - self.apertures[ia]) + \
                               self.apertures[ia]

            calc = (required < self.flux[:, 0])
            radius[calc] = self.apertures[0]

            calc = (required >= self.flux[:, -1])
            radius[calc] = self.apertures[-1]

            return radius

    def find_radius_sigma(self, fraction):

        log.info("Finding radius containing %g%s of the flux" % (fraction * 100., '%'))

        sigma = np.zeros(self.flux.shape)
        sigma[:, 0] = self.flux[:, 0] / self.apertures[0] ** 2
        sigma[:, 1:] = (self.flux[:, 1:] - self.flux[:, :-1]) / \
                       (self.apertures[1:] ** 2 - self.apertures[:-1] ** 2)

        maximum = np.max(sigma, axis=1)

        radius = np.zeros(self.n_models)

        # Linear interpolation - need to loop over apertures backwards for vectorization
        for ia in range(len(self.apertures) - 2, -1, -1):
            calc = (sigma[:, ia] > fraction * maximum) & (radius == 0.)
            radius[calc] = (sigma[calc, ia] - fraction * maximum[calc]) / \
                           (sigma[calc, ia] - sigma[calc, ia + 1]) * \
                           (self.apertures[ia + 1] - self.apertures[ia]) + \
                           self.apertures[ia]

        calc = sigma[:, -1] > fraction * maximum
        radius[calc] = self.apertures[-1]

        return radius
