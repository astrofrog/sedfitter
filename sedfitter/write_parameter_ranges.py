from __future__ import print_function, division

import numpy as np

from .fit_info import FitInfoFile
from .models import load_parameter_table, load_model_apertures
from .parameter_columns import (aperture_dependent_columns, aperture_to_au,
                                parameter_names, parameter_values,
                                validate_aperture)

__all__ = ['write_parameter_ranges']

NODATA = '-'.center(10)


def write_parameter_ranges(input_fits, output_file, select_format=("N", 1), additional={},
                           aperture=None):
    """
    Write out an ASCII file with ranges of parameters for each source.

    Parameters
    ----------
    input_fits : str or :class:`sedfitter.fit_info.FitInfo` or iterable
        This should be either a file containing the fit information, a
        :class:`sedfitter.fit_info.FitInfo` instance, or an iterable containing
        :class:`sedfitter.fit_info.FitInfo` instances.
    output_file : str, optional
        The output ASCII file containing the parameter ranges
    select_format : tuple, optional
        Tuple specifying which fits should be output. See the documentation
        for a description of the tuple syntax.
    additional : dict, optional
        A dictionary giving additional parameters for each model. This should
        be a dictionary where each key is a parameter, and each value is a
        dictionary mapping the model names to the parameter values.
    aperture : `~astropy.units.Quantity`, optional
        For models with aperture-dependent parameters (such as the Richardson
        et al. 2024 YSO models), by default the ranges of these parameters are
        written out separately for each model aperture. If an aperture radius
        is specified, the parameters are instead interpolated to that
        aperture, and the range of apertures used is given in an
        ``aperture_au`` column. The aperture can be given as a length or as an
        angle, in which case it is converted to a length using the distance of
        each fit.
    """

    aperture = validate_aperture(aperture)

    # Open input and output file
    fin = FitInfoFile(input_fits, 'r')
    fout = open(output_file, 'w')

    # Read in table of parameters for model grid
    t = load_parameter_table(fin.meta.model_dir)

    t['MODEL_NAME'] = np.char.strip(t['MODEL_NAME'])
    t.sort('MODEL_NAME')

    apertures = load_model_apertures(fin.meta.model_dir)
    names = parameter_names(t, apertures, interpolate=aperture is not None)
    names += list(additional.keys())

    # Only write out the aperture column if any parameters depend on it
    write_aperture = aperture is not None and len(aperture_dependent_columns(t, apertures)) > 0
    if write_aperture:
        names.insert(0, 'aperture_au')

    # First header line
    fout.write("%30s " % "")

    fout.write("%10s " % "")
    fout.write("%10s " % "")

    fout.write('chi2'.center(32) + ' ')
    fout.write('av'.center(32) + ' ')
    fout.write('scale'.center(32) + ' ')

    for par in names:
        fout.write(par.lower().center(32) + ' ')

    fout.write('\n')

    # Second header line
    fout.write("source_name".center(30) + ' ')
    fout.write("n_data".center(10) + ' ')
    fout.write("n_fits".center(10) + ' ')

    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")
    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")
    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")

    for par in names:
        fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")

    fout.write('\n')

    # Third header line
    fout.write('-' * 30 + ' ')
    fout.write('-' * 10 + ' ')
    fout.write('-' * 10 + ' ')

    fout.write('-' * 32 + ' ')
    fout.write('-' * 32 + ' ')
    fout.write('-' * 32 + ' ')

    for par in names:
        fout.write('-' * 32 + ' ')

    fout.write('\n')

    for info in fin:

        # Filter fits
        info.keep(select_format)

        # Get filtered and sorted table of parameters
        tsorted = info.filter_table(t, additional=additional)

        if aperture is None:
            aperture_au = None
        else:
            aperture_au = aperture_to_au(aperture, info.sc)

        values = parameter_values(tsorted, apertures, aperture_au=aperture_au)

        if write_aperture:
            values.insert(0, aperture_au.value)

        fout.write("%30s " % info.source.name)
        fout.write("%10i " % info.source.n_data)
        fout.write("%10i " % info.n_fits)

        if len(info.chi2) == 0:
            fout.write('%10s %10s %10s ' % (NODATA, NODATA, NODATA))
            fout.write('%10s %10s %10s ' % (NODATA, NODATA, NODATA))
            fout.write('%10s %10s %10s ' % (NODATA, NODATA, NODATA))
        else:
            fout.write('%10.3e %10.3e %10.3e ' % (np.nanmin(info.chi2), info.chi2[0], np.nanmax(info.chi2)))
            fout.write('%10.3e %10.3e %10.3e ' % (np.nanmin(info.av), info.av[0], np.nanmax(info.av)))
            fout.write('%10.3e %10.3e %10.3e ' % (np.nanmin(info.sc), info.sc[0], np.nanmax(info.sc)))

        for column in values:
            if len(info.chi2) == 0:
                fout.write('%10s %10s %10s ' % (NODATA, NODATA, NODATA))
            else:
                fout.write('%10.3e %10.3e %10.3e ' % (np.nanmin(column), column[0], np.nanmax(column)))

        fout.write('\n')

    # Close input and output files
    fin.close()
    fout.close()
    fout.close()
