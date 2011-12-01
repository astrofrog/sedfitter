#!/usr/bin/env python

import os
import cPickle as pickle

import numpy as np

import atpy

from fit_info import FitInfo
from extinction import Extinction


def write_parameters(input_file, output_file, select_format=("N", 1), additional={}):

    # Open input and output file
    fin = file(input_file, 'rb')
    fout = file(output_file, 'wb')

    # Read in header of output file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction = Extinction()
    extinction.read_binary(fin)

    if os.path.exists(model_dir + '/parameters.fits'):
        t = atpy.Table(model_dir + '/parameters.fits')
    elif os.path.exists(model_dir + '/parameters.fits.gz'):
        t = atpy.Table(model_dir + '/parameters.fits.gz')
    else:
        raise Exception("Parameter file not found in %s" % model_dir)

    t['MODEL_NAME'] = np.char.strip(t['MODEL_NAME'])
    t.sort('MODEL_NAME')

    info = FitInfo()

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

    for par in t.columns.keys + additional.keys():
        if par == 'MODEL_NAME':
            continue
        fout.write(par.lower().center(10) + ' ')

    fout.write('\n')

    fout.write('-' * (75 + 11 * (len(t.columns.keys + additional.keys()))))
    fout.write('\n')

    while True:

        # Read in next fit
        try:
            info.read_binary(fin)
        except:
            break

        # Filter fits
        info.keep(select_format[0], select_format[1])

        subset = np.in1d(t['MODEL_NAME'], info.model_name)
        tsub = t.where(subset)
        index = np.argsort(np.argsort(info.model_name))
        tsorted = tsub.rows(index)
        if not np.all(info.model_name == tsorted['MODEL_NAME']):
            raise Exception("Parameter file sorting failed")

        # Add additional parameter columns if necessary
        for parameter in additional:
            if parameter in tsorted.columns:
                raise Exception("Parameter {} already exists in file".format(parameter))
            tsorted.add_empty_column(parameter, dtype=float)
            for i, name in enumerate(tsorted['MODEL_NAME']):
                tsorted[parameter][i] = additional[parameter][name.strip()]

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
