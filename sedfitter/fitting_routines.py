import numpy as np

def linear_regression(data, weights, pattern1, pattern2):

    c1 = np.sum(data*pattern1*weights, axis=1)
    c2 = np.sum(data*pattern2*weights, axis=1)
    m11 = np.sum(pattern1*pattern1*weights, axis=1)
    m12 = np.sum(pattern1*pattern2*weights, axis=1)
    m22 = np.sum(pattern2*pattern2*weights, axis=1)

    inv_det = 1./(m11*m22-m12*m12)

    p1 = (m22*c1-m12*c2)*inv_det
    p2 = (m11*c2-m12*c1)*inv_det

    return p1, p2


def optimal_scaling(data, weights, pattern1):

    return np.sum(data*pattern1*weights) / np.sum(pattern1*pattern1*weights)


def chi_squared(valid, data, error, weight, model):

    chi2_array = np.zeros(data.shape, dtype=np.float32)

    elem = (valid==1) | (valid==4)
    chi2_array[elem] = (data[elem] - model[elem])**2 * weight[elem]

    elem = ((valid==2) & (data > model)) | ((valid==3) & (data > model))
    hard = error==1.
    chi2_array[elem & hard] = 1.e30
    chi2_array[elem & ~hard] = -2. * np.log10(1.-error[elem & ~hard])

    return np.sum(chi2_array, axis=1)