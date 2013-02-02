from __future__ import print_function, division

try:
    import cPickle as pickle
except:
    import pickle
import string

import numpy as np
import atpy

from .fit_info import FitInfo
from .extinction import Extinction


def extract_parameters(input=None, output_prefix=None, output_suffix=None,
                       parameters='all', select_format=("N", 1), header=True):

    if input[-8:] != '.fitinfo':
        raise Exception("Extension of input file should be .fitinfo")

    fin = file(input, 'rb')

    # Read in header of output file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction_law = Extinction()
    extinction_law.read_binary(fin)

    # Read in parameters file
    tpar = atpy.Table(model_dir + '/parameters.fits.gz')
    try:
        par_model_names = np.char.strip(tpar.MODEL_NAME)
    except:
        par_model_names = np.array([x.strip() for x in tpar.MODEL_NAME],
                                   dtype=tpar.MODEL_NAME.dtype)
    tpar.MODEL_NAME = np.char.strip(tpar.MODEL_NAME)

    format = {}
    for par in tpar.names:
        if par == 'MODEL_NAME':
            format[par] = "%11s"
        else:
            format[par] = "%11.3e"

    if parameters == 'all':
        parameters = tpar.names

    info = FitInfo()

    while True:

        # Read in next fit
        try:
            info.read_binary(fin)
        except:
            break

        # Filter fits
        info.keep(select_format[0], select_format[1])

        output = ""
        if output_prefix:
            output += output_prefix
        output += info.source.name
        if output_suffix:
            output += output_suffix

        fout = file(output, 'wb')

        if header:
            fout.write("%11s %11s %11s " % ("CHI2", "AV", "SC"))
            fout.write(string.join([("%11s" % par) for par in parameters], " "))
            fout.write("\n")

        for i in range(info.n_fits):

            row = tpar.where(info.model_name[i] == par_model_names).row(0)

            basic = "%11.3e %11.3e %11.3e " % (info.chi2[i], info.av[i], info.sc[i])
            pars = string.join([(format[par] % row[par]) for par in parameters], " ")

            fout.write(basic + pars + "\n")

        fout.close()

    # Close input and output files
    fin.close()
