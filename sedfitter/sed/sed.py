from __future__ import print_function, division

import os

import numpy as np
from astropy import log
from astropy.io import fits
from scipy.interpolate import interp1d

from ..utils.validator import validate_array

C = 299792458
KPC = 3.086e21


class SED(object):

    def __init__(self):
        self.wav = None
        self.nu = None
        self.apertures = None
        self.flux = None
        self.error = None

    def copy(self):
        from copy import deepcopy
        return deepcopy(self)

    def scale_to_distance(self, distance):
        """
        Returns the SED scaled to distance `distance`
        
        Parameters
        ----------
        distance : float
            The distance in cm

        Returns
        -------
        sed : SED
            The SED, scaled to the new distance
        """
        sed = self.copy()
        sed.distance = distance
        sed.flux *= (self.distance / distance) ** 2
        sed.error *= (self.distance / distance) ** 2
        return sed

    def scale_to_av(self, av, law):
        sed = self.copy()
        sed.flux *= 10. ** (av * law(sed.wav))
        sed.error *= 10. ** (av * law(sed.wav))

    @property
    def wav(self):
        """
        The wavelengths at which the SED is defined
        """
        return self._wav

    @wav.setter
    def wav(self, value):
        if value is None:
            self._wav = None
        else:
            self._wav = validate_array('wav', value, domain='positive', ndim=1, shape=None if self.nu is None else (len(self.nu),))

    @property
    def nu(self):
        """
        The frequencies at which the SED is defined
        """
        return self._nu

    @nu.setter
    def nu(self, value):
        if value is None:
            self._nu = None
        else:
            self._nu = validate_array('nu', value, domain='positive', ndim=1, shape=None if self.wav is None else (len(self.wav),))

    @property
    def apertures(self):
        """
        The apertures at which the SED is defined
        """
        return self._apertures

    @apertures.setter
    def apertures(self, value):
        if value is None:
            self._apertures = None
        else:
            self._apertures = validate_array('apertures', value, domain='positive', ndim=1)

    @property
    def flux(self):
        """
        The SED fluxes
        """
        return self._flux

    @flux.setter
    def flux(self, value):
        if value is None:
            self._flux = value
        else:

            if self.n_ap is not None:
                self._flux = validate_array('flux', value, ndim=2, shape=(self.n_ap, self.n_wav))
            else:
                self._flux = validate_array('flux', value, ndim=1, shape=(self.n_wav, ))

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

            if self.n_ap is not None:
                self._error = validate_array('error', value, ndim=2, shape=(self.n_ap, self.n_wav))
            else:
                self._error = validate_array('error', value, ndim=1, shape=(self.n_wav, ))

    @property
    def n_ap(self):
        if self.apertures is None:
            return None
        else:
            return len(self.apertures)

    @property
    def n_wav(self):
        if self.wav is None:
            return None
        else:
            return len(self.wav)

    @classmethod
    def read(cls, filename, unit_wav='microns', unit_freq='Hz',
             unit_flux='ergs/cm^2/s', order='freq'):
        '''
        Read an SED from a FITS file.

        Parameters
        ----------
        filename: str
            The name of the file to read the SED from.
        unit_wav: str, optional
            The units to convert the wavelengths to.
        unit_freq: str, optional
            The units to convert the frequency to.
        unit_flux: str, optional
            The units to convert the flux to.
        order: str, optional
            Whether to sort the SED by increasing wavelength (`wav`) or
            frequency ('freq').
        '''

        # Instantiate SED class
        sed = cls()

        # Assume that the filename may be missing the .gz extension
        if not os.path.exists(filename) and os.path.exists(filename + '.gz'):
            filename += ".gz"

        # Open FILE file
        hdulist = fits.open(filename, memmap=False)

        # Extract model name
        sed.name = hdulist[0].header['MODEL']

        # Check if distance is specified in header, otherwise assume 1kpc
        if 'DISTANCE' in hdulist[0].header:
            sed.distance = hdulist[0].header['DISTANCE']
        else:
            log.debug("No distance found in SED file, assuming 1kpc")
            sed.distance = 3.086e21  # cm (=1kpc)

        # Extract SED values
        wav = hdulist[1].data.field('WAVELENGTH')
        nu = hdulist[1].data.field('FREQUENCY')
        ap = hdulist[2].data.field('APERTURE')
        flux = hdulist[3].data.field('TOTAL_FLUX')
        err = hdulist[3].data.field('TOTAL_FLUX_ERR')

        # Save apertures
        sed.apertures = ap

        # Extract units
        curr_unit_wav = hdulist[1].columns[0].unit.lower()
        curr_unit_freq = hdulist[1].columns[1].unit.lower()
        curr_unit_flux = hdulist[3].columns[0].unit.lower()

        # Convert requested units to lowercase for comparison
        unit_wav = unit_wav.lower()
        unit_freq = unit_freq.lower()
        unit_flux = unit_flux.lower()

        # Convert wavelength to microns
        if curr_unit_wav == 'm':
            wav_microns = wav * 1.e6
        elif curr_unit_wav == 'nm':
            wav_microns = wav / 1000.
        elif curr_unit_wav == 'microns':
            wav_microns = wav
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_wav, unit_wav))

        # Convert wavelength from microns to requested units
        if unit_wav == 'microns':
            sed.wav = wav_microns
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_wav, unit_wav))

        # Convert frequency to Hz
        if curr_unit_freq == 'hz':
            nu_hz = nu
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_freq, unit_freq))

        # Convert frequency from Hz to requested units
        if unit_freq == 'hz':
            sed.nu = nu_hz
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_freq, unit_freq))

        # Convert flux to ergs/cm^2/s
        if curr_unit_flux == 'ergs/s':
            flux_cgs = flux / sed.distance ** 2.
            error_cgs = err / sed.distance ** 2.
        elif curr_unit_flux == 'mjy':
            flux_cgs = flux * C / (wav_microns * 1.e-6) * 1.e-26
            error_cgs = err * C / (wav_microns * 1.e-6) * 1.e-26
        elif curr_unit_flux == 'ergs/cm^2/s':
            flux_cgs = flux
            error_cgs = err
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_flux, unit_flux))

        # Convert flux from ergs/cm^2/s to requested units
        if unit_flux == 'ergs/s':
            sed.flux = flux_cgs * sed.distance ** 2.
            sed.error = error_cgs * sed.distance ** 2.
        elif unit_flux == 'mjy':
            sed.flux = flux_cgs / C * (wav_microns * 1.e-6) / 1.e-26
            sed.error = error_cgs / C * (wav_microns * 1.e-6) / 1.e-26
        elif unit_flux == 'ergs/cm^2/s':
            sed.flux = flux_cgs
            sed.error = error_cgs
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_flux, unit_flux))

        # Sort SED
        if (order == 'freq' and sed.nu[0] > sed.nu[-1]) or \
           (order == 'wav' and sed.wav[0] > sed.wav[-1]):
            sed.wav = sed.wav[::-1]
            sed.nu = sed.nu[::-1]
            sed.flux = sed.flux[..., ::-1]
            sed.error = sed.error[..., ::-1]

        return sed

    def interpolate(self, apertures):
        '''
        Interpolate the SED to different apertures
        '''

        # If there is only one aperture, we can't interpolate, we can only repeat
        if self.n_ap == 1:
            return np.repeat(self.flux[0, :], len(apertures)).reshape(self.n_wav, len(apertures))

        # Create interpolating function
        flux_interp = interp1d(self.ap, self.flux.swapaxes(0, 1))

        # If any apertures are larger than the defined max, reset to max
        apertures[apertures > self.ap.max()] = self.ap.max()

        # If any apertures are smaller than the defined min, raise Exception
        if np.any(apertures < self.ap.min()):
            raise Exception("Aperture(s) requested too small")

        return flux_interp(apertures)

    def interpolate_variable(self, wavelengths, apertures):
        '''
        Interpolate the SED to a variable aperture as a function of
        wavelength. This method should be called with an interpolating
        function for aperture as a function of wavelength, in log10 space.
        '''

        # If any apertures are larger than the defined max, reset to max
        apertures[apertures > self.ap.max()] = self.ap.max() * 0.999

        # If any apertures are smaller than the defined min, raise Exception
        if np.any(apertures < self.ap.min()):
            raise Exception("Aperture(s) requested too small")

        if self.n_ap == 1:
            return self.flux[0, :]

        # Find wavelength order
        order = np.argsort(wavelengths)

        # Interpolate apertures vs wavelength
        log10_ap_interp = interp1d(np.log10(wavelengths[order]), np.log10(apertures[order]), bounds_error=False, fill_value=np.nan)

        # Create interpolating function
        flux_interp = interp1d(self.ap, self.flux.swapaxes(0, 1))

        # Interpolate the apertures
        apertures = 10. ** log10_ap_interp(np.log10(self.wav))

        # Extrapolate on either side
        apertures[np.log10(self.wav) < log10_ap_interp.x[0]] = 10. ** log10_ap_interp.y[0]
        apertures[np.log10(self.wav) > log10_ap_interp.x[-1]] = 10. ** log10_ap_interp.y[-1]

        # Interpolate and return only diagonal elements
        return flux_interp(apertures).diagonal()
