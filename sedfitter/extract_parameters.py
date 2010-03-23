import cPickle as pickle
import string

import numpy as np
import atpy

from sedfitter.fit_info import FitInfo
from sedfitter.extinction import Extinction


def extract_parameters(input=None, output_prefix=None, output_suffix=None, parameters=[], select_format=("N", 1)):

    if input[-8:] <> '.fitinfo':
        raise Exception("Extension of input file should be .fitinfo")

    fin = file(input, 'rb')

    # Read in header of output file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction_law = Extinction()
    extinction_law.read_binary(fin)

    # Read in parameters file
    tpar = atpy.Table(model_dir + '/parameters.fits.gz')
    par_model_names = np.char.strip(tpar.MODEL_NAME)

    format = {}
    for par in tpar.names:
        if par=='MODEL_NAME':
            format[par] = "%s"
        else:
            format[par] = "%11.3e"

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

        for i in range(info.n_fits):

                row = tpar.where(info.model_name[i] == par_model_names).row(0)
                
                basic = "%11.3e %11.3e %11.3e " % (info.chi2[i], info.av[i], info.sc[i]) 
                pars = string.join([(format[par] % row[par]) for par in parameters]," ")
 
                fout.write(basic + pars + "\n")

        fout.close()

    # Close input and output files
    fin.close()
