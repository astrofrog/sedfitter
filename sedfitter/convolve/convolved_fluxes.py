from __future__ import print_function, division

import numpy as np
np.seterr(all='ignore')

from scipy.interpolate import interp1d

from astropy.logger import log

# TODO: get rid of use of interp1d

class ConvolvedFluxes(object):

    def __init__(self, n_models=None, n_ap=1, wavelength=None, dtype=np.float32):

        self.wavelength = wavelength

        if n_models is not None:

            self.model_names = np.zeros(n_models, dtype='S30')
            self.apertures = np.zeros(n_ap, dtype=dtype)

            if n_ap == 1:
                self.flux = np.zeros(n_models, dtype=dtype)
                self.err = np.zeros(n_models, dtype=dtype)
                self.flux_interp = None
            else:
                self.flux = np.zeros((n_models, n_ap), dtype=dtype)
                self.err = np.zeros((n_models, n_ap), dtype=dtype)
                self.flux_interp = interp1d(self.apertures, self.flux[:])

        else:

            self.model_names = None
            self.apertures = None
            self.flux = None
            self.err = None
            self.flux_interp = None

    @property
    def n_models(self):
        if self.model_names is None:
            return None
        else:
            return self.model_names.shape[0]

    @property
    def n_ap(self):
        if self.model_names is None:
            return None
        elif self.flux.ndim == 1:
            return 1
        else:
            return self.flux.shape

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

        # Read in number of models and apertures
        conv.n_models = keywords['NMODELS']
        conv.n_ap = keywords['NAP']

        # Create shortcuts to tables
        tc = Table(convolved['CONVOLVED FLUXES'].data)
        ta = Table(convolved['APERTURES'].data)

        # Read in model names
        conv.model_names = tc['MODEL_NAME']

        # Read in flux and flux errors
        conv.flux = tc['TOTAL_FLUX']
        conv.err = tc['TOTAL_FLUX_ERR']

        # Read in 99% cumulative and 50% surface brightness radii
        try:
            conv.radius_sigma_50 = tc['RADIUS_SIGMA_50']
            conv.radius_cumul_99 = tc['RADIUS_CUMUL_99']
        except KeyError:
            pass

        # Read in apertures
        conv.apertures = ta['APERTURE']

        # Create an interpolating function for the flux vs aperture
        if conv.flux.ndim > 1:
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

        keywords = {}
        keywords['FILTWAV'] = self.wavelength
        keywords['NMODELS'] = self.n_models
        keywords['NAP'] = self.n_ap

        tc = Table()
        tc.add_column(Column('MODEL_NAME', self.model_names))
        tc.add_column(Column('TOTAL_FLUX', self.flux))
        tc.add_column(Column('TOTAL_FLUX_ERR', self.err))

        if self.n_ap > 1:
            tc.add_column(Column('RADIUS_SIGMA_50', self.find_radius_sigma(0.50)))
            tc.add_column(Column('RADIUS_CUMUL_99', self.find_radius_cumul(0.99)))

        ta = Table()
        ta.add_column(Column("APERTURE", self.apertures))

        # Convert to FITS
        hdulist = fits.HDUList()
        hdulist.append(fits.PrimaryHDU(header=fits.Header(keywords)))
        hdulist.append(fits.BinTableHDU(np.array(tc), name='CONVOLVED FLUXES'))
        hdulist.append(fits.BinTableHDU(np.array(ta), name='APERTURES'))
        hdulist.writeto(filename, clobber=overwrite)

    def interpolate(self, apertures):
        '''
        Interpolate the flux to the apertures specified.
        '''

        # Initalize new ConvolvedFluxes object to return
        c = ConvolvedFluxes(n_models=self.n_models,
                            n_ap=len(apertures),
                            wavelength=self.wavelength)

        # Save requested apertures
        c.apertures = apertures[:]

        # Transfer model names
        c.model_names = self.model_names

        # If any apertures are larger than the defined max, reset to max
        apertures[apertures > self.apertures.max()] = self.apertures.max()

        # If any apertures are smaller than the defined min, raise Exception
        if np.any(apertures < self.apertures.min()):
            raise Exception("Aperture(s) requested too small")

        # Interpolate to requested apertures
        if self.flux.ndim > 1:
            c.flux = self.flux_interp(apertures)
        else:
            c.flux = np.repeat(self.flux, len(apertures)).reshape(len(c.flux), len(apertures))

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
