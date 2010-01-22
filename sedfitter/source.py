import numpy as np


def log_fluxes(valid, flux, error):

    # Initialize arrays
    logflux = np.zeros(flux.shape, dtype=np.float32)
    logerror = np.zeros(error.shape, dtype=np.float32)
    weight = np.zeros(valid.shape, dtype=np.float32)

    # Fluxes
    r = valid == 1
    logflux[r] = np.log10(flux[r]) - 0.5 * (error[r]/flux[r]) ** 2. / np.log(10.)
    logerror[r] = np.abs(error[r]/flux[r]) / np.log(10.)
    weight[r] = 1./logerror[r]**2.

    # Lower and upper limits
    r = (valid == 2) | (valid == 3)
    logflux[r] = np.log10(flux[r])
    logerror[r] = error[r]

    # Log10[Fluxes]
    logflux[r] = flux[r]
    logerror[r] = error[r]
    weight[r] = 1./logerror[r]**2.

    return weight, logflux, logerror


class Source(object):

    def __init__(self, name, x, y, valid, flux, error):

        self.name = name
        self.x = x
        self.y = y
        self.valid = valid
        self.flux = flux
        self.error = error

        self.weight, self.logflux, self.logerror = log_fluxes(valid, flux, error)

        self.n_wav = len(self.valid)
        self.n_data = len(np.where((self.valid==1) | (self.valid==4))[0])

        return

    def __str__(self):

        string = "Source name : %s\n" % self.name
        string += "RA   / l    : %9.5f\n" % self.x
        string += "Decl / b    : %9.5f\n" % self.y
        for j in range(self.n_wav):
            string += "F = %12.4e +/- %12.4e mJy (%1i)  Log[F] = %8.5f+/-%8.5f\n" % \
                      (self.flux[j], self.error[j], self.valid[j], self.logflux[j], self.logerror[j])

        return string


def read_sources(filename, n_min_valid=0):

    result = []

    data = np.loadtxt(filename, dtype=str)

    try:
        n_sources = np.shape(data)[0]
        n_wav = (np.shape(data)[1]-3)/3
    except:
        n_sources = 1
        n_wav = (len(data)-3)/3

    n = n_wav

    name = np.loadtxt(filename, usecols=[0], dtype=str)
    x = np.loadtxt(filename, usecols=[1], dtype=np.float32)
    y = np.loadtxt(filename, usecols=[2], dtype=np.float32)
    valid = np.loadtxt(filename, usecols=range(3, 3+n), dtype=np.int32)
    flux = np.loadtxt(filename, usecols=range(3+n, 3+3*n-1, 2), dtype=np.float32)
    error = np.loadtxt(filename, usecols=range(3+n+1, 3+3*n, 2), dtype=np.float32)

    sources = []
    for i in range(n_sources):
        if np.sum((valid[i,:]==1) | (valid[i,:]==4)) > n_min_valid:
            sources.append(Source(name[i], x[i], y[i], valid[i,:], flux[i,:], error[i,:]))

    return sources

#def write(self, filename, mask=None):
#
#    f = file(filename, 'w')
#
#    if mask==None:
#        mask = np.ones((self.n_sources), bool)
#
#    for i in range(self.n_sources):
#        if mask[i]:
#            f.write("%30s " % self.name[i])
#            f.write("%9.5f %9.5f " % (self.x[i], self.y[i]))
#            f.write("%1i "*self.n_wav % tuple(self.valid[i, :].tolist()))
#            for j in range(self.n_wav):
#                f.write("%11.3e %11.3e "% (self.flux[i, j], self.error[i, j]))
#            f.write("\n")
#
#    f.close()
#
#    return
