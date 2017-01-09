from __future__ import print_function, division

import numpy as np

from .fit_info import FitInfoFile
from .models import load_parameter_table

__all__ = ['write_parameter_ranges']

NODATA = '-'.center(10)


def write_parameter_ranges(input_fits, output_file, select_format=("N", 1), additional={}):
    """
    Write out an ASCII file with ranges of paramters for each source.

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
    """

    # Open input and output file
    fin = FitInfoFile(input_fits, 'r')
    fout = open(output_file, 'w')

    # Read in table of parameters for model grid
    t = load_parameter_table(fin.meta.model_dir)

    t['MODEL_NAME'] = np.char.strip(t['MODEL_NAME'])
    t.sort('MODEL_NAME')

    # First header line
    fout.write("%30s " % "")

    fout.write("%10s " % "")
    fout.write("%10s " % "")

    fout.write('chi2'.center(32) + ' ')
    fout.write('av'.center(32) + ' ')
    fout.write('scale'.center(32) + ' ')

    for par in list(t.columns.keys()) + list(additional.keys()):
        if par == 'MODEL_NAME':
            continue
        fout.write(par.lower().center(32) + ' ')

    fout.write('\n')

    # Second header line
    fout.write("source_name".center(30) + ' ')
    fout.write("n_data".center(10) + ' ')
    fout.write("n_fits".center(10) + ' ')

    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")
    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")
    fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")

    for par in list(t.columns.keys()) + list(additional.keys()):
        if par == 'MODEL_NAME':
            continue
        fout.write("min".center(10) + " " + "best".center(10) + " " + "max".center(10) + " ")

    fout.write('\n')

    # Third header line
    fout.write('-' * 30 + ' ')
    fout.write('-' * 10 + ' ')
    fout.write('-' * 10 + ' ')

    fout.write('-' * 32 + ' ')
    fout.write('-' * 32 + ' ')
    fout.write('-' * 32 + ' ')

    for par in list(t.columns.keys()) + list(additional.keys()):
        if par == 'MODEL_NAME':
            continue
        fout.write('-' * 32 + ' ')

    fout.write('\n')

    for info in fin:

        # Filter fits
        info.keep(select_format)

        # Get filtered and sorted table of parameters
        tsorted = info.filter_table(t, additional=additional)

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

        for par in tsorted.columns:
            if par == 'MODEL_NAME':
                continue
            if len(info.chi2) == 0:
                fout.write('%10s %10s %10s ' % (NODATA, NODATA, NODATA))
            else:
                fout.write('%10.3e %10.3e %10.3e ' % (np.nanmin(tsorted[par]), tsorted[par][0], np.nanmax(tsorted[par])))

        fout.write('\n')

    # Close input and output files
    fin.close()
    fout.close()
    fout.close()
