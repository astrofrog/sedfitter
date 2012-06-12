from __future__ import print_function, division

from scipy import weave
from scipy.weave import converters


def chi_squared_weave(valid, data, error, weight, model):

    nx = valid.shape[0]
    ny = valid.shape[1]

    chi2_array = np.zeros(data.shape, dtype=np.float32)

    code = """
    for (int i=0; i < nx; i++)
    {
        for (int j=0; j < ny; j++)
        {
            if (valid(i, j)==1 | valid(i, j)==4) {
                chi2_array(i, j) = ( data(i, j) - model(i, j) ) * ( data(i, j) - model(i, j) )* weight(i, j) ;
            } else if ((valid(i, j)==2 & data(i, j) > model(i, j)) | (valid(i, j)==3 & data(i, j) < model(i, j))) {
                if (error(i, j)==1.) {
                    chi2_array(i, j) = 1.e30;
                } else {
                    chi2_array(i, j) = -2. * log10(1-error(i, j));
                }
            }
        }
    }
    """

    weave.inline(code, ['valid', 'data', 'error', 'weight', 'model', 'nx', 'ny', 'chi2_array'], type_converters=converters.blitz)

    return np.sum(chi2_array, axis=1)

    # nx = weights.shape[0]
    # ny = weights.shape[1]
    # c1 = np.zeros(nx)
    # c2 = np.zeros(nx)
    # m11 = np.zeros(nx)
    # m12 = np.zeros(nx)
    # m22 = np.zeros(nx)
    #
    # code = """
    # for (int i=0; i < nx; i++)
    # {
    #     for (int j=0; j < ny; j++)
    #     {
    #         c1(i) = c1(i) + data(i, j) * pattern1(j) * weights(i, j);
    #         c2(i) = c2(i) + data(i, j) * pattern2(j) * weights(i, j);
    #         m11(i) = m11(i) + pattern1(j) * pattern1(j) * weights(i, j);
    #         m12(i) = m12(i) + pattern1(j) * pattern2(j) * weights(i, j);
    #         m22(i) = m22(i) + pattern2(j) * pattern2(j) * weights(i, j);
    #     }
    # }
    # """
    #
    # weave.inline(code, ['data', 'weights', 'pattern1', 'pattern2', 'nx', 'ny', 'c1', 'c2', 'm11', 'm12', 'm22'], type_converters=converters.blitz)
