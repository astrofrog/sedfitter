from __future__ import print_function, division

import numpy as np

from .fit_info import FitInfoFile
from .models import load_parameter_table, load_model_apertures
from .parameter_columns import (aperture_dependent_columns, aperture_to_au,
                                parameter_names, parameter_values,
                                validate_aperture)

__all__ = ['write_parameters']


def write_parameters(input_fits, output_file, select_format=("N", 1), additional={},
                     aperture=None):
    """
    Write out an ASCII file with the parameters for each source.

    Parameters
    ----------
    input_fits : str or :class:`sedfitter.fit_info.FitInfo` or iterable
        This should be either a file containing the fit information, a
        :class:`sedfitter.fit_info.FitInfo` instance, or an iterable containing
        :class:`sedfitter.fit_info.FitInfo` instances.
    output_file : str, optional
        The output ASCII file containing the parameters
    select_format : tuple, optional
        Tuple specifying which fits should be output. See the documentation
        for a description of the tuple syntax.
    additional : dict, optional
        A dictionary giving additional parameters for each model. This should
        be a dictionary where each key is a parameter, and each value is a
        dictionary mapping the model names to the parameter values.
    aperture : `~astropy.units.Quantity`, optional
        For models with aperture-dependent parameters (such as the Richardson
        et al. 2024 YSO models), by default these parameters are written out
        with one column per model aperture. If an aperture radius is
        specified, the parameters are instead interpolated to that aperture,
        and an ``aperture_au`` column gives the aperture used for each fit.
        The aperture can be given as a length or as an angle, in which case
        it is converted to a length using the distance of each fit.
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

    # First header line
    fout.write("source_name".center(30) + ' ')
    fout.write("n_data".center(10) + ' ')
    fout.write("n_fits".center(10) + ' ')

    fout.write('\n')

    # Second header line

    fout.write('fit_id'.center(10) + ' ')
    fout.write('model_name'.center(30) + ' ')
    fout.write('chi2'.center(10) + ' ')
    fout.write('av'.center(10) + ' ')
    fout.write('scale'.center(10) + ' ')

    if write_aperture:
        fout.write('aperture_au'.center(10) + ' ')

    for par in names:
        fout.write(par.lower().center(10) + ' ')

    fout.write('\n')

    fout.write('-' * (75 + 11 * (len(names) + write_aperture + 1)))
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

        fout.write("%30s " % info.source.name)
        fout.write("%10i " % info.source.n_data)
        fout.write("%10i " % info.n_fits)
        fout.write("\n")

        for fit_id in range(len(info.chi2)):
            fout.write('%10i ' % (fit_id + 1))
            fout.write('%30s ' % info.model_name[fit_id])
            fout.write('%10.3f ' % info.chi2[fit_id])
            fout.write('%10.3f ' % info.av[fit_id])
            fout.write('%10.3f ' % info.sc[fit_id])

            if write_aperture:
                fout.write('%10.3e ' % aperture_au[fit_id].value)

            for column in values:
                fout.write('%10.3e ' % column[fit_id])

            fout.write('\n')

    # Close input and output files
    fin.close()
    fout.close()
    fout.close()
