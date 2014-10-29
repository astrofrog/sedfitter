import numpy as np
from astropy import units as u
from astropy.io import fits


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


def assert_allclose_quantity(q1, q2):
    if q1 is None and q2 is None:
        return True
    if q1 is None or q2 is None:
        raise AssertionError()
    else:
        np.testing.assert_allclose(q1.value, q2.to(q1.unit).value)


def convert_flux(nu, flux, target_unit, distance=None):
    """
    Convert flux to a target unit
    """

    curr_unit = flux.unit

    if curr_unit.is_equivalent(u.erg / u.s):
        flux = flux / distance ** 2
    elif curr_unit.is_equivalent(u.Jy):
        flux = flux * nu
    elif not curr_unit.is_equivalent(u.erg / u.cm ** 2 / u.s):
        raise Exception("Don't know how to convert {0} to ergs/cm^2/s" % (flux.unit))

    # Convert to requested unit

    if target_unit.is_equivalent(u.erg / u.s):
        flux = flux * distance ** 2
    elif target_unit.is_equivalent(u.Jy):
        flux = flux / nu
    elif not target_unit.is_equivalent(u.erg / u.cm ** 2 / u.s):
        raise Exception("Don't know how to convert %s to %s" % (curr_unit, target_unit))

    return flux.to(target_unit)
    
def table_to_hdu(table):
    hdu = fits.BinTableHDU(np.array(table))
    for i in range(len(table.columns)):
        if table.columns[i].unit is not None:
            hdu.columns[i].unit = table.columns[i].unit.to_string(format='fits')
    return hdu