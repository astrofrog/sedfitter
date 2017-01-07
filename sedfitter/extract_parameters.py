from __future__ import print_function, division

import numpy as np

from .fit_info import FitInfoFile
from .models import load_parameter_table

__all__ = ['extract_parameters']


def extract_parameters(input=None, output_prefix=None, output_suffix=None,
                       parameters='all', select_format=("N", 1), header=True):

    fin = FitInfoFile(input, 'r')

    # Read in table of parameters for model grid
    t = load_parameter_table(fin.meta.model_dir)

    format = {}
    for par in t.dtype.names:
        if par == 'MODEL_NAME':
            format[par] = "%11s"
        else:
            format[par] = "%11.3e"

    if parameters == 'all':
        parameters = t.dtype.names

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
            fout.write(" ".join([("%11s" % par) for par in parameters]))
            fout.write("\n")

        # Get filtered and sorted table of parameters
        tsorted = info.filter_table(t)

        for i in range(info.n_fits):

            row = tsorted[i]

            basic = "%11.3e %11.3e %11.3e " % (info.chi2[i], info.av[i], info.sc[i])
            pars = " ".join([(format[par] % row[par]) for par in parameters])

            fout.write(basic + pars + "\n")

        fout.close()

    # Close input and output files
    fin.close()
