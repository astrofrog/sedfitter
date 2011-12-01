# Still to implement:
# - Everything

import cPickle as pickle

from sedfitter.fit_info import FitInfo
from sedfitter.extinction import Extinction


def filter_output(input=None, output_good='auto', output_bad='auto', chi=None,
                  cpd=None):

    if input[-8:] != '.fitinfo':
        raise Exception("Extension of input file should be .fitinfo")

    fin = file(input, 'rb')

    if output_good == 'auto':
        fout_good = file(input.replace('.fitinfo', '_good.fitinfo'), 'wb')
    else:
        fout_good = file(output_good, 'wb')

    if output_bad == 'auto':
        fout_bad = file(input.replace('.fitinfo', '_bad.fitinfo'), 'wb')
    else:
        fout_bad = file(output_bad, 'wb')

    # Read in header of output file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction_law = Extinction()
    extinction_law.read_binary(fin)

    # Output header to good fits file
    pickle.dump(model_dir, fout_good, 2)
    pickle.dump(filters, fout_good, 2)
    extinction_law.write_binary(fout_good)

    # Output header to bad fits file
    pickle.dump(model_dir, fout_bad, 2)
    pickle.dump(filters, fout_bad, 2)
    extinction_law.write_binary(fout_bad)

    info = FitInfo()

    while True:

        # Read in next fit
        try:
            info.read_binary(fin)
        except:
            break

        bestchi = info.chi2[0]
        bestcpd = info.chi2[0] / float(info.source.n_data)

        if (chi and bestchi < chi) or (cpd and bestcpd < cpd):
            info.write_binary(fout_good)
        else:
            info.write_binary(fout_bad)

    # Close input and output files
    fin.close()
    fout_good.close()
    fout_bad.close()
