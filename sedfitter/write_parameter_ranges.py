from __future__ import print_function, division

import numpy as np

from copy import deepcopy

from .fit_info import FitInfoFile
from .models import load_parameter_table

__all__ = ['write_parameter_ranges']

NODATA = '-'.center(10)

#this exists for compatibility with the Richardson+ (2024) YSO models;
#it identifies aperture-dependent parameters with array values and
#pulls out a single value for a particular aperture
def standard_col(table,parameter,aperture):
    column = deepcopy(table[parameter])
    ndim = len(column.data.shape)
    if ndim == 2:
        if column.data.shape[-1] > 1:
            column = column[:,aperture]
        else:
            column = column[:,0]
    elif ndim == 3:
        column = column[:,0,aperture]
    return column

def write_parameter_ranges(input_fits, output_file, select_format=("N", 1), additional={}, aperture=None):
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
    aperture : int, optional
        The index of values to return for table columns with array values.
        Defaults to 5, corresponding to an aperture of radius ~1000 AU.
        Intended for use with the 'Richardson+ (2024) YSO SED models:
        <https://zenodo.org/records/10522816>'
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

    #set aperture to ~1000 AU if none is picked
    if aperture is None:
        aperture = 5

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
                col = standard_col(tsorted,par,aperture)
                
                fout.write('%10.3e %10.3e %10.3e ' % (np.nanmin(col), col[0], np.nanmax(col)))

        fout.write('\n')

    # Close input and output files
    fin.close()
    fout.close()
    fout.close()
