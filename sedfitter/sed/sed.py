import os

from scipy.interpolate import interp1d

# import atpy
import pyfits

from copy import copy
import numpy as np

c = 299792458


class SED(object):

    def __init__(self):
        pass

    def scale_to_distance(self, distance):
        self.distance = distance
        self._update_sed()

    def scale_to_av(self, av, law):
        self.av = av
        self.law = law
        self._update_sed()

    def _update_sed(self):
        self.flux = copy(self._flux)
        self.err = copy(self._err)
        if 'distance' in self.__dict__:
            self.flux /= self.distance ** 2
            self.err /= self.distance ** 2
        if 'law' in self.__dict__:
            self.flux *= 10. ** (self.av * self.law(self.wav))
            self.err *= 10. ** (self.av * self.law(self.wav))

    def read(self, filename, unit_wav='microns', unit_freq='Hz', unit_flux='ergs/cm^2/s', order='freq'):
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

        if not os.path.exists(filename) and os.path.exists(filename + '.gz'):
            filename += ".gz"

        hdulist = pyfits.open(filename, memmap=False)

        self.name = hdulist[0].header['MODEL']

        wav = hdulist[1].data.field('WAVELENGTH')
        nu = hdulist[1].data.field('FREQUENCY')
        ap = hdulist[2].data.field('APERTURE')
        flux = hdulist[3].data.field('TOTAL_FLUX')
        err = hdulist[3].data.field('TOTAL_FLUX_ERR')

        self.n_wav = len(wav)
        self.n_ap = len(ap)

        unit_wav = unit_wav.lower()
        unit_freq = unit_freq.lower()
        unit_flux = unit_flux.lower()

        curr_unit_wav = hdulist[1].columns[0].unit.lower()
        curr_unit_freq = hdulist[1].columns[1].unit.lower()
        curr_unit_flux = hdulist[3].columns[0].unit.lower()

        # ts = atpy.TableSet(filename, verbose=False)
        #
        # self.n_wav = len(ts[0])
        # self.n_ap = len(ts[1])
        #
        # self._wav = ts[0].WAVELENGTH
        # self.ap = ts[1].APERTURE
        # self._flux = ts[2].TOTAL_FLUX
        #
        # curr_unit_wav = ts[0].columns['WAVELENGTH'].unit
        # curr_unit_flux = ts[2].columns['TOTAL_FLUX'].unit

        # Convert wavelength to microns
        if curr_unit_wav == 'm':
            wav_microns = wav * 1.e6
        elif curr_unit_wav == 'nm':
            wav_microns = wav / 1000.
        elif curr_unit_wav == 'microns':
            wav_microns = wav
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_wav, unit_wav))

        # Convert wavelength to requested units
        if unit_wav == 'microns':
            self.wav = wav_microns
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_wav, unit_wav))

        # Convert frequency to Hz
        if curr_unit_freq == 'hz':
            nu_hz = nu
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_freq, unit_freq))

        # Convert frequency to requested units
        if unit_freq == 'hz':
            self.nu = nu_hz
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_freq, unit_freq))

        # Convert flux to ergs/cm^2/s
        if curr_unit_flux == 'ergs/s':
            flux_cgs = flux / 3.086e21 ** 2.
            err_cgs = err / 3.086e21 ** 2.
        elif curr_unit_flux == 'mjy':
            flux_cgs = flux * c / (wav_microns * 1.e-6) * 1.e-26
            err_cgs = err * c / (wav_microns * 1.e-6) * 1.e-26
        elif curr_unit_flux == 'ergs/cm^2/s':
            flux_cgs = flux
            err_cgs = err
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_flux, unit_flux))

        if unit_flux == 'ergs/s':
            self._flux = flux_cgs * 3.086e21 ** 2.
            self._err = err_cgs * 3.086e21 ** 2.
        elif unit_flux == 'mjy':
            self._flux = flux_cgs / c * (wav_microns * 1.e-6) / 1.e-26
            self._err = err_cgs / c * (wav_microns * 1.e-6) / 1.e-26
        elif unit_flux == 'ergs/cm^2/s':
            self._flux = flux_cgs
            self._err = err_cgs
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_flux, unit_flux))

        # Sort SED
        if order == 'freq' and self.nu[0] > self.nu[-1] or \
           order == 'wav' and self.wav[0] > self.wav[-1]:
            self.wav = self.wav[::-1]
            self.nu = self.nu[::-1]
            self._flux = self._flux[::-1]
            self._err = self._err[::-1]

        # Initialize distance and Av
        self.distance = 1.
        self.av = 0.
        self.ap = ap
        self.flux = copy(self._flux)
        self.err = copy(self._err)

    def interpolate(self, apertures):
        '''
        Interpolate the SED to different apertures
        '''

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
        apertures[apertures > self.ap.max()] = self.ap.max()

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
