from __future__ import print_function, division

import numpy as np

from .fit_info import FitInfoFile
from .models import load_parameter_table

__all__ = ['write_parameters']


def write_parameters(input_fits, output_file, select_format=("N", 1), additional={}):
    """
    Write out an ASCII file with the paramters for each source.

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
    """

    # Open input and output file
    fin = FitInfoFile(input_fits, 'r')
    fout = open(output_file, 'w')

    # Read in table of parameters for model grid
    t = load_parameter_table(fin.meta.model_dir)

    t['MODEL_NAME'] = np.char.strip(t['MODEL_NAME'])
    t.sort('MODEL_NAME')

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

    for par in list(t.columns.keys()) + list(additional.keys()):
        if par == 'MODEL_NAME':
            continue
        fout.write(par.lower().center(10) + ' ')

    fout.write('\n')

    fout.write('-' * (75 + 11 * (len(list(t.columns.keys()) + list(additional.keys())))))
    fout.write('\n')

    for info in fin:

        # Filter fits
        info.keep(select_format)

        # Get filtered and sorted table of parameters
        tsorted = info.filter_table(t, additional=additional)

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

            for par in tsorted.columns:
                if par == 'MODEL_NAME':
                    continue
                fout.write('%10.3e ' % (tsorted[par][fit_id]))

            fout.write('\n')

    # Close input and output files
    fin.close()
    fout.close()
    fout.close()
