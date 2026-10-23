from __future__ import print_function, division

import numpy as np

from .fit_info import FitInfoFile
from .models import load_parameter_table, load_model_apertures
from .parameter_columns import (aperture_dependent_columns, aperture_to_au,
                                parameter_names, parameter_values,
                                validate_aperture)

__all__ = ['extract_parameters']


def extract_parameters(input=None, output_prefix=None, output_suffix=None,
                       parameters='all', select_format=("N", 1), header=True,
                       aperture=None):
    """
    Write out one ASCII file per source with the parameters of each fit.

    For models with aperture-dependent parameters (such as the Richardson et
    al. 2024 YSO models), these parameters are written out with one column
    per model aperture by default. If ``aperture`` is specified (as a length,
    or as an angle which is converted to a length using the distance of each
    fit), they are instead interpolated to that aperture, and an
    ``APERTURE_AU`` column gives the aperture used for each fit.
    """

    aperture = validate_aperture(aperture)

    fin = FitInfoFile(input, 'r')

    # Read in table of parameters for model grid
    t = load_parameter_table(fin.meta.model_dir)

    if parameters == 'all':
        parameters = t.colnames

    write_name = 'MODEL_NAME' in parameters
    t = t[['MODEL_NAME'] + [par for par in parameters if par != 'MODEL_NAME']]

    apertures = load_model_apertures(fin.meta.model_dir)
    names = parameter_names(t, apertures, interpolate=aperture is not None)

    # Only write out the aperture column if any parameters depend on it
    write_aperture = aperture is not None and len(aperture_dependent_columns(t, apertures)) > 0
    if write_aperture:
        names.insert(0, 'APERTURE_AU')

    for info in fin:

        # Filter fits
        info.keep(select_format)

        output = ""
        if output_prefix:
            output += output_prefix
        output += info.source.name
        if output_suffix:
            output += output_suffix

        fout = open(output, 'w')

        if header:
            fout.write("%11s %11s %11s " % ("CHI2", "AV", "SC"))
            fout.write(" ".join([("%11s" % par) for par in ['MODEL_NAME'] * write_name + names]))
            fout.write("\n")

        # Get filtered and sorted table of parameters
        tsorted = info.filter_table(t)

        if aperture is None:
            aperture_au = None
        else:
            aperture_au = aperture_to_au(aperture, info.sc)

        values = parameter_values(tsorted, apertures, aperture_au=aperture_au)

        if write_aperture:
            values.insert(0, aperture_au.value)

        for i in range(info.n_fits):

            basic = "%11.3e %11.3e %11.3e " % (info.chi2[i], info.av[i], info.sc[i])
            pars = ["%11s" % tsorted['MODEL_NAME'][i]] * write_name
            pars += [("%11.3e" % column[i]) for column in values]

            fout.write(basic + " ".join(pars) + "\n")

        fout.close()

    # Close input and output files
    fin.close()
