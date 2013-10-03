from __future__ import print_function, division

import numpy as np
np.seterr(all='ignore')

from scipy.interpolate import interp1d

from astropy.logger import log
from astropy.utils.misc import isiterable
from ..utils.validator import validate_scalar, validate_array

# TODO: get rid of use of interp1d


def is_numpy_array(variable):
    return issubclass(variable.__class__, (np.ndarray,
                                           np.core.records.recarray,
                                           np.ma.core.MaskedArray))


class ConvolvedFluxes(object):

    def __init__(self, wavelength=None, model_names=None, apertures=None, flux=None, error=None, initialize_arrays=False):

        self.model_names = model_names
        self.apertures = apertures
        self.wavelength = wavelength

        if initialize_arrays:

            if model_names is None:
                raise ValueError("model_names is required when using initialize_arrays=True")

            if apertures is None:
                raise ValueError("apertures is required when using initialize_arrays=True")

            if flux is None:
                self.flux = np.zeros((self.n_models, self.n_ap))
            else:
                self.flux = flux

            if error is None:
                self.error = np.zeros((self.n_models, self.n_ap))
            else:
                self.error = error

        else:

            self.flux = flux
            self.error = error

    @property
    def wavelength(self):
        """
        The effective wavelength at which the fluxes are defined
        """
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        if value is None:
            self._wavelength = None
        else:
            self._wavelength = validate_scalar('wavelength', value, domain='positive')

    @property
    def model_names(self):
        """
        The names of the models
        """
        return self._model_names

    @model_names.setter
    def model_names(self, value):
        if value is None:
            self._model_names = value
        else:
            self._model_names = validate_array('model_names', value, ndim=1)

    @property
    def apertures(self):
        """
        The apertures for which the models are defined
        """
        return self._apertures

    @apertures.setter
    def apertures(self, value):
        if value is None:
            self._apertures = value
        else:
            self._apertures = validate_array('apertures', value, ndim=1)

    @property
    def flux(self):
        """
        The convolved fluxes
        """
        return self._flux

    @flux.setter
    def flux(self, value):
        if value is None:
            self._flux = value
        else:

            if self.model_names is None:
                raise ValueError("model_names has not been set")

            self._flux = validate_array('flux', value, ndim=2, shape=(self.n_models, self.n_ap))

    @property
    def error(self):
        """
        The convolved flux errors
        """
        return self._error

    @error.setter
    def error(self, value):
        if value is None:
            self._error = value
        else:

            if self.model_names is None:
                raise ValueError("model_names has not been set")

            self._error = validate_array('error', value, ndim=2, shape=(self.n_models, self.n_ap))


    @property
    def n_models(self):
        if self.model_names is None:
            return None
        else:
            return self.model_names.shape[0]

    @property
    def n_ap(self):
        if self.apertures is None:
            return 1
        else:
            return len(self.apertures)

    def __eq__(self, other):
        return self.wavelength == other.wavelength \
            and np.all(self.model_names == other.model_names) \
            and np.all(self.apertures == other.apertures) \
            and np.all(self.flux == other.flux) \
            and np.all(self.error == other.error)

    @classmethod
    def read(cls, filename):
        '''
        Read convolved flux from a FITS file

        Parameters
        ----------
        filename : str
            The name of the FITS file to read the convolved fluxes from
        '''

        from astropy.io import fits
        from astropy.table import Table

        conv = cls()

        # Open the convolved flux FITS file
        convolved = fits.open(filename)
        keywords = convolved[0].header

        # Try and read in the wavelength of the filter
        if 'FILTWAV' in keywords:
            conv.wavelength = keywords['FILTWAV']
        else:
            conv.wavelength = None

        # Read in apertures, if present
        try:
            ta = Table(convolved['APERTURES'].data)
            conv.apertures = ta['APERTURE']
        except KeyError:
            pass

        # Create shortcuts to table
        tc = Table(convolved['CONVOLVED FLUXES'].data)

        # Read in model names
        conv.model_names = tc['MODEL_NAME']

        # Read in flux and flux errors
        if tc['TOTAL_FLUX'].ndim == 1 and conv.n_ap == 1:
            conv.flux = tc['TOTAL_FLUX'].reshape(tc['TOTAL_FLUX'].shape[0], 1)
        else:
            conv.flux = tc['TOTAL_FLUX']
        if tc['TOTAL_FLUX_ERR'].ndim == 1 and conv.n_ap == 1:
            conv.error = tc['TOTAL_FLUX_ERR'].reshape(tc['TOTAL_FLUX_ERR'].shape[0], 1)
        else:
            conv.error = tc['TOTAL_FLUX_ERR']

        # Read in 99% cumulative and 50% surface brightness radii
        try:
            conv.radius_sigma_50 = tc['RADIUS_SIGMA_50']
            conv.radius_cumul_99 = tc['RADIUS_CUMUL_99']
        except KeyError:
            pass

        # Create an interpolating function for the flux vs aperture
        if conv.n_ap > 1:
            conv.flux_interp = interp1d(conv.apertures, conv.flux[:])

        return conv

    def write(self, filename, overwrite=False):
        '''
        Write convolved flux to a FITS file.

        Parameters
        ----------
        filename: str
            The name of the file to output the convolved fluxes to.
        overwrite: bool, optional
            Whether to overwrite the output file
        '''

        from astropy.io import fits
        from astropy.table import Table, Column

        tc = Table()
        tc['MODEL_NAME'] = self.model_names
        tc['TOTAL_FLUX'] = self.flux
        tc['TOTAL_FLUX_ERR'] = self.error

        if self.apertures is not None:
            tc['RADIUS_SIGMA_50'] = self.find_radius_sigma(0.50)
            tc['RADIUS_CUMUL_99'] = self.find_radius_cumul(0.99)

        if self.apertures is not None:
            ta = Table()
            ta['APERTURE'] = self.apertures

        # Convert to FITS
        hdulist = fits.HDUList()
        hdulist.append(fits.PrimaryHDU())
        if self.wavelength is not None:
            hdulist[0].header['FILTWAV'] = self.wavelength
        hdulist[0].header['NMODELS'] = self.n_models
        hdulist[0].header['NAP'] = self.n_ap

        hdulist.append(fits.BinTableHDU(np.array(tc), name='CONVOLVED FLUXES'))
        if self.apertures is not None:
            hdulist.append(fits.BinTableHDU(np.array(ta), name='APERTURES'))
        hdulist.writeto(filename, clobber=overwrite)

    def interpolate(self, apertures):
        '''
        Interpolate the flux to the apertures specified.
        '''

        # Initalize new ConvolvedFluxes object to return
        c = ConvolvedFluxes()

        # Set wavelength
        c.wavelength = self.wavelength

        # Save requested apertures
        c.apertures = apertures[:]

        # Transfer model names
        c.model_names = self.model_names

        # Interpolate to requested apertures
        if self.n_ap > 1:

            # If any apertures are larger than the defined max, reset to max
            if np.any(apertures > self.apertures.max()):
                apertures[apertures > self.apertures.max()] = self.apertures.max()

            # If any apertures are smaller than the defined min, raise error
            if np.any(apertures < self.apertures.min()):
                raise Exception("Aperture(s) requested too small")

            flux_interp = interp1d(self.apertures, self.flux)
            c.flux = flux_interp(apertures)

            # The following is not strictly correct - errors from interpolation is not interpolation of errors
            error_interp = interp1d(self.apertures, self.error)
            c.error = error_interp(apertures)

        else:

            c.flux = np.repeat(self.flux, len(apertures)).reshape(c.n_models, len(apertures))
            c.error = np.repeat(self.error, len(apertures)).reshape(c.n_models, len(apertures))

        try:
            c.radius_sigma_50 = self.radius_sigma_50[:]
            c.radius_cumul_99 = self.radius_cumul_99[:]
        except AttributeError:
            pass

        return c

    def find_radius_cumul(self, fraction):
        '''
        Find for each model the radius containing a fraction of the flux.

        Parameters
        ----------
        fraction: float
            The fraction to use when determining the radius
        '''

        log.info("Calculating radii containing %g%s of the flux" % (fraction * 100., '%'))

        radius = np.zeros(self.n_models, dtype=self.flux.dtype)

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
        '''
        Find for each model a fractional surface brightness radius

        This is the outermost radius where the surface brightness is larger
        than a fraction of the peak surface brightness.

        Parameters
        ----------
        fraction: float
            The fraction to use when determining the radius
        '''

        log.info("Calculating %g%s peak surface brightness radii" % (fraction * 100., '%'))

        sigma = np.zeros(self.flux.shape, dtype=self.flux.dtype)
        sigma[:, 0] = self.flux[:, 0] / self.apertures[0] ** 2
        sigma[:, 1:] = (self.flux[:, 1:] - self.flux[:, :-1]) / \
                       (self.apertures[1:] ** 2 - self.apertures[:-1] ** 2)

        maximum = np.max(sigma, axis=1)

        radius = np.zeros(self.n_models, dtype=self.flux.dtype)

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
