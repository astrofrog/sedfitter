import abc
import numpy as np

from astropy import units as u
from astropy.io import fits
from astropy.extern import six
from astropy.table import Table

from ..utils.validator import validate_scalar, validate_array


__all__ = ['SEDCube', 'PolarizationCube']


UNIT_MAPPING = {}
UNIT_MAPPING['MICRONS'] = u.micron
UNIT_MAPPING['HZ'] = u.Hz
UNIT_MAPPING['MJY'] = u.mJy
UNIT_MAPPING['ergs/cm^2/s'] = u.erg / u.cm ** 2 / u.s


def parse_unit_safe(unit_string):
    if unit_string in UNIT_MAPPING:
        return UNIT_MAPPING[unit_string]
    else:
        return u.Unit(unit_string, parse_strict=False)


def table_to_hdu(table):
    hdu = fits.BinTableHDU(np.array(table))
    for i in range(len(table.columns)):
        if table.columns[i].unit is not None:
            hdu.columns[i].unit = table.columns[i].unit.to_string(format='fits')
    return hdu


def assert_allclose_quantity(q1, q2):
    if q1 is None and q2 is None:
        return True
    if q1 is None or q2 is None:
        raise AssertionError("One of q1 or q2 is None")
    else:
        np.testing.assert_allclose(q1.value, q2.to(q1.unit).value)


@six.add_metaclass(abc.ABCMeta)
class BaseCube(object):
    """
    A cube to represent a cube of models.

    This consists of values and uncertainties as a function of wavelength,
    aperture, and models.

    Parameters
    ----------
    names : 1-d iterable, optional
        The names of all the models in the cube
    distance : `~astropy.units.Quantity`, optional
        The distance assumed for the values
    wav : 1-d `~astropy.units.Quantity`, optional
        The wavelengths at which the SEDs are defined (cannot be used with ``nu``)
    nu : 1-d `~astropy.units.Quantity`, optional
        The frequencies at which the SEDs are defined (cannot be used with ``wav``)
    apertures : 1-d `~astropy.units.Quantity`, optional
        The ap for which the SEDs are defined
    val : 3-d `~astropy.units.Quantity`, optional
         The values of the fluxes or polarization
    unc : 3-d `~astropy.units.Quantity`, optional
         The uncertainties in the fluxes or polarization
    """

    _physical_type = None

    def __init__(self, names=None, distance=None, wav=None, nu=None,
                 apertures=None, val=None, unc=None):

        # The names of all the models
        self.names = names

        # The distance at which the fluxes are defined
        self.distance = distance

        # The wavelengths and ap
        self.wav = wav
        self.nu = nu
        self.apertures = apertures

        # The value and uncertainties
        self.val = val
        self.unc = unc

    def __eq__(self, other):

        try:

            assert np.all(self.names == other.names)

            assert_allclose_quantity(self.distance, other.distance)

            assert_allclose_quantity(self.wav, other.wav)
            assert_allclose_quantity(self.nu, other.nu)

            assert_allclose_quantity(self.apertures, other.apertures)

            assert_allclose_quantity(self.val, other.val)
            assert_allclose_quantity(self.unc, other.unc)

        except AssertionError:
            raise
            return False
        else:
            return True

    @property
    def names(self):
        """
        The names of the models
        """
        return self._names

    @names.setter
    def names(self, value):
        if value is None:
            self._names = None
        else:
            if isinstance(value, np.ndarray):
                if value.dtype.kind == 'U':
                    value = np.char.encode(value, 'ascii')
            else:
                value = np.array(value, 'S')
            self._names = value

    @property
    def wav(self):
        """
        The wavelengths at which the SEDs are defined.
        """
        if self._wav is None and self._nu is not None:
            return self._nu.to(u.micron, equivalencies=u.spectral())
        else:
            return self._wav

    @wav.setter
    def wav(self, value):
        if value is None:
            self._wav = None
        else:
            self._nu = None
            self._wav = validate_array('wav', value, domain='positive', ndim=1,
                                       shape=None if self.nu is None else (len(self.nu),),
                                       physical_type='length')

    @property
    def nu(self):
        """
        The frequencies at which the SEDs are defined.
        """
        if self._nu is None and self._wav is not None:
            return self._wav.to(u.Hz, equivalencies=u.spectral())
        else:
            return self._nu

    @nu.setter
    def nu(self, value):
        if value is None:
            self._nu = None
        else:
            self._wav = None
            self._nu = validate_array('nu', value, domain='positive', ndim=1,
                                      shape=None if self.wav is None else (len(self.wav),),
                                      physical_type='frequency')

    @property
    def apertures(self):
        """
        The ap at which the SEDs are defined.
        """
        return self._apertures

    @apertures.setter
    def apertures(self, value):
        if value is None:
            self._apertures = None
        else:
            self._apertures = validate_array('apertures', value, domain='positive',
                                             ndim=1, physical_type='length')

    @property
    def distance(self):
        """
        The distance at which the SEDs are defined.
        """
        return self._distance

    @distance.setter
    def distance(self, value):
        if value is None:
            self._distance = None
        else:
            self._distance = validate_scalar('distance', value, domain='positive',
                                             physical_type='length')

    @property
    def val(self):
        """
        The fluxes or polarization values.
        """
        return self._val

    @val.setter
    def val(self, value):
        if value is None:
            self._val = value
        else:
            self._val = validate_array('val', value, ndim=3,
                                       shape=(self.n_ap, self.n_wav, self.n_models),
                                       physical_type=self._physical_type)

    @property
    def unc(self):
        """
        The uncertainties in the fluxes or polarization.
        """
        return self._unc

    @unc.setter
    def unc(self, value):
        if value is None:
            self._unc = value
        else:
            self._unc = validate_array('unc', value, ndim=3,
                                       shape=(self.n_ap, self.n_wav, self.n_models),
                                       physical_type=self._physical_type)

    @property
    def n_ap(self):
        if self.apertures is None:
            return 1
        else:
            return len(self.apertures)

    @property
    def n_wav(self):
        if self.wav is None:
            return None
        else:
            return len(self.wav)

    @property
    def n_models(self):
        if self.names is None:
            return None
        else:
            return len(self.names)

    @classmethod
    def read(cls, filename):
        """
        Read models from a FITS file.

        Parameters
        ----------
        filename: str
            The name of the file to read the cube from.
        """

        # Create class instance
        cube = cls()

        # Open FILE file
        hdulist = fits.open(filename, memmap=True)

        # Extract distance
        cube.distance = hdulist[0].header['DISTANCE'] * u.cm

        # Extract model names
        cube.names = hdulist['MODEL_NAMES'].data['MODEL_NAME']

        # Extract wavelengths
        cube.wav = hdulist['SPECTRAL_INFO'].data['WAVELENGTH'] * parse_unit_safe(hdulist[2].columns[0].unit)
        cube.nu = hdulist['SPECTRAL_INFO'].data['FREQUENCY'] * parse_unit_safe(hdulist[2].columns[1].unit)

        # Extract ap
        cube.apertures = hdulist['APERTURES'].data['APERTURE'] * parse_unit_safe(hdulist[3].columns[0].unit)

        # Extract value
        cube.val = hdulist['VALUES'].data * parse_unit_safe(hdulist[4].header['BUNIT'])

        # Extract uncertainty
        cube.unc = hdulist['UNCERTAINTIES'].data * parse_unit_safe(hdulist[5].header['BUNIT'])

        return cube

    def _check_all_set(self):
        if self.wav is None:
            raise ValueError("Wavelengths 'wav' are not set")
        if self.nu is None:
            raise ValueError("Frequencies 'nu' are not set")
        if self.val is None:
            raise ValueError("Values 'val' are not set")
        if self.distance is None:
            raise ValueError("Value 'distance' is not set")

    def write(self, filename, overwrite=False):
        """
        Write the models to a FITS file.

        Parameters
        ----------
        filename: str
            The name of the file to write the cube to.
        """

        self._check_all_set()

        hdulist = fits.HDUList()

        # Create empty first HDU and add distance
        hdu0 = fits.PrimaryHDU()
        hdu0.header['distance'] = (self.distance.to(u.cm).value, 'Distance assumed for the values, in cm')
        hdu0.header['NWAV'] = (self.n_wav, "Number of wavelengths")
        if self.apertures is not None:
            hdu0.header['NAP'] = (self.n_ap, "Number of apertures")
        hdulist.append(hdu0)

        # Create names table
        t1 = Table()
        t1['MODEL_NAME'] = self.names
        hdu1 = table_to_hdu(t1)
        hdu1.name = "MODEL_NAMES"
        hdulist.append(hdu1)

        # Create wavelength table
        t2 = Table()
        t2['WAVELENGTH'] = self.wav
        t2['FREQUENCY'] = self.nu
        hdu2 = table_to_hdu(t2)
        hdu2.name = "SPECTRAL_INFO"
        hdulist.append(hdu2)

        # Create aperture table
        if self.apertures is not None:
            t3 = Table()
            t3['APERTURE'] = self.apertures
            hdu3 = table_to_hdu(t3)
            hdu3.name = "APERTURES"
            hdulist.append(hdu3)

        # Create value HDU
        hdu4 = fits.ImageHDU(self.val.value)
        hdu4.header['BUNIT'] = self.val.unit.to_string()
        hdu4.name = 'VALUES'
        hdulist.append(hdu4)

        # Create uncertainty HDU
        if self.unc is not None:
            hdu5 = fits.ImageHDU(self.unc.value)
            hdu5.header['BUNIT'] = self.unc.unit.to_string()
            hdu5.name = 'UNCERTAINTIES'
            hdulist.append(hdu5)

        # Write out HDUList
        hdulist.writeto(filename, clobber=overwrite)


class SEDCube(BaseCube):
    _physical_type = ('power', 'flux', 'spectral flux density')


class PolarizationCube(BaseCube):
    _physical_type = ('dimensionless')
