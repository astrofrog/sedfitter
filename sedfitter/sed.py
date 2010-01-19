import atpy
from copy import copy


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

        # Read file

        ts = atpy.TableSet(filename, verbose=False)

        self.n_wav = len(ts[0])
        self.n_ap = len(ts[1])

        self._wav = ts[0].WAVELENGTH
        self.ap = ts[1].APERTURE
        self._flux = ts[2].TOTAL_FLUX

        curr_unit_wav = ts[0].columns['WAVELENGTH'].unit
        curr_unit_flux = ts[2].columns['TOTAL_FLUX'].unit

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
