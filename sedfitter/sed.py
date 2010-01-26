import os

from scipy.interpolate import interp1d

# import atpy
import pyfits

from copy import copy
import numpy as np

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
        if 'distance' in self.__dict__:
            self.flux /= self.distance**2
        if 'law' in self.__dict__:
            self.flux *= 10.**(self.av * self.law(self._wav))

    def read(self, filename, unit_wav='microns', unit_flux='ergs/cm^2/s'):

        if not os.path.exists(filename) and os.path.exists(filename + '.gz'):
            filename += ".gz"

        hdulist = pyfits.open(filename)
        self._wav = hdulist[1].data.field('WAVELENGTH')
        self.ap = hdulist[2].data.field('APERTURE')
        self._flux = hdulist[3].data.field('TOTAL_FLUX')
        
        self.n_wav = len(self._wav)
        self.n_ap = len(self.ap)

        curr_unit_wav = hdulist[1].columns[0].unit
        curr_unit_flux = hdulist[3].columns[0].unit

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

        # Convert wavelength units

        if unit_wav == 'microns':
            if curr_unit_wav.lower() == 'm':
                self._wav *= 1.e6
            elif curr_unit_wav.lower() == 'nm':
                self._wav /= 1000.
            elif curr_unit_wav.lower() == 'microns':
                pass
            else:
                raise Exception("Don't know how to convert %s to %s" % (curr_unit_wav, unit_wav))
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_wav, unit_wav))

        # Convert flux units

        if unit_flux == 'ergs/cm^2/s':
            if curr_unit_flux.lower() == 'ergs/s':
                self._flux /= 3.086e21**2.
            elif curr_unit_flux.lower() == 'mjy':
                self._flux *= 3.e8 / (self._wav * 1.e-6) * 1.e-26
            elif curr_unit_flux.lower() == 'ergs/cm^2/s':
                pass
            else:
                raise Exception("Don't know how to convert %s to %s" % (curr_unit_flux, unit_flux))
        else:
            raise Exception("Don't know how to convert %s to %s" % (curr_unit_flux, unit_flux))

        # Initialize distance and Av

        self.distance = 1.
        self.av = 0.
        self.wav = self._wav
        self.flux = self._flux

    def interpolate(self, apertures):
        '''
        Interpolate the SED to different apertures
        '''

        if self.n_ap==1:
            return np.repeat(self.flux[0,:],len(apertures)).reshape(self.n_wav, len(apertures))

        # Create interpolating function
        flux_interp = interp1d(self.ap, self.flux.swapaxes(0,1))

        return flux_interp(apertures)

    def interpolate_variable(self, wavelengths, apertures):
        '''
        Interpolate the SED to a variable aperture as a function of
        wavelength. This method should be called with an interpolating
        function for aperture as a function of wavelength, in log10 space.
        '''

        if np.any(apertures < self.ap.min()):
            raise Exception("Aperture(s) requested too small")

        if self.n_ap == 1:
            return self.flux[0,:]

        # Find wavelength order
        order = np.argsort(wavelengths)

        # Interpolate apertures vs wavelength
        log10_ap_interp = interp1d(np.log10(wavelengths[order]), np.log10(apertures[order]), bounds_error=False, fill_value=np.nan)

        # Create interpolating function
        flux_interp = interp1d(self.ap, self.flux.swapaxes(0,1))

        # Interpolate the apertures
        apertures = 10.**log10_ap_interp(np.log10(self.wav))

        # Extrapolate on either side
        apertures[np.log10(self.wav) < log10_ap_interp.x[0]] = self.ap[0]
        apertures[np.log10(self.wav) > log10_ap_interp.x[-1]] = self.ap[-1]

        # Interpolate and return only diagonal elements
        return flux_interp(apertures).diagonal()
