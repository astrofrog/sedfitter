from __future__ import print_function, division

from .fit_info import FitInfoFile
from .extinction import Extinction


def filter_output(input_file=None, output_good='auto', output_bad='auto', chi=None,
                  cpd=None):
    """
    Filter an output file into well and badly fit sources.

    The sources in the output file are split into ones where the best fit is
    good, and those where the best fit is bad, where the distinction between
    good and bad is made based on the chi^2 value or the chi^2 value per
    datapoint of the best fit.

    Parameters
    ----------
    input_file : str
        File containing the fit information
    output_good : str, optional
        The name of the file containing information about well-fit sources. If
        set to 'auto', then the output file will be set to the same as the
        input file but with a ``_good`` prefix.
    output_bad : str, optional
        The name of the file containing information about badly-fit sources. If
        set to 'auto', then the output file will be set to the same as the
        input file but with a ``_bad`` prefix.
    chi : float, optional
        If set, the separation between well and badly fit sources is done based
        on the total chi^2 of the best fit, using this value as a threshold.
    cpd : float, optional
        If set, the separation between well and badly fit sources is done based
        on the chi^2 per datapoint of the best fit, using this value as a
        threshold.
    """

    fin = FitInfoFile.open(input_file, 'r')

    if output_good == 'auto':
        fout_good = FitInfoFile.open(input_file + '_good', 'w')
    else:
        fout_good = FitInfoFile.open(output_good, 'w')

    if output_bad == 'auto':
        fout_bad = FitInfoFile.open(input_file + '_bad', 'w')
    else:
        fout_bad = FitInfoFile.open(output_bad, 'w')

    for info in fin:

        bestchi = info.chi2[0]
        bestcpd = info.chi2[0] / float(info.source.n_data)

        if (chi and bestchi < chi) or (cpd and bestcpd < cpd):
            fout_good.write(info)
        else:
            fout_bad.write(info)

    # Close input and output files
    fin.close()
    fout_good.close()
    fout_bad.close()
