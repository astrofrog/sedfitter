import convolved_fluxes as c
import numpy as np
import os,sys
from scipy import weave
from scipy.weave import converters

def fit(model_flux,valid,flux,flux_error,weight,av_law,sc_law):
        
    flux = np.tile(flux,model_flux.shape[0]).reshape(model_flux.shape)
    flux_error = np.tile(flux_error,model_flux.shape[0]).reshape(model_flux.shape)
    valid = np.tile(valid,model_flux.shape[0]).reshape(model_flux.shape)
    weight = np.tile(weight,model_flux.shape[0]).reshape(model_flux.shape)
                
    residual = flux - model_flux
    
    av_best,sc_best = linear_regression(residual,weight,av_law,sc_law)
    
    model = av_best[:,np.newaxis] * av_law[np.newaxis,:] + sc_best[:,np.newaxis] * sc_law[np.newaxis,:]
    
    ch_best = chi_squared(valid,residual,flux_error,weight,model)
    
    return av_best,sc_best,ch_best
    
def linear_regression(data,weights,pattern1,pattern2):
    
    nx = weights.shape[0]
    ny = weights.shape[1]

    c1 = np.zeros(nx)
    c2 = np.zeros(nx)
    m11 = np.zeros(nx)
    m12 = np.zeros(nx)
    m22 = np.zeros(nx)
    
#    print nx
    
    code = """
    for (int i=0; i < nx; i++)
    {
        for (int j=0; j < ny; j++)
        {
            c1(i) = c1(i) + data(i,j) * pattern1(j) * weights(i,j);
            c2(i) = c2(i) + data(i,j) * pattern2(j) * weights(i,j);
            m11(i) = m11(i) + pattern1(j) * pattern1(j) * weights(i,j);
            m12(i) = m12(i) + pattern1(j) * pattern2(j) * weights(i,j);
            m22(i) = m22(i) + pattern2(j) * pattern2(j) * weights(i,j);
        }
    }
    """
    
    weave.inline(code, ['data','weights','pattern1','pattern2','nx','ny','c1','c2','m11','m12','m22'], type_converters=converters.blitz)
    
#    c1 = np.sum(data*pattern1*weights,axis=1)
#    c2 = np.sum(data*pattern2*weights,axis=1)
#    m11 = np.sum(pattern1*pattern1*weights,axis=1)
#    m12 = np.sum(pattern1*pattern2*weights,axis=1)
#    m22 = np.sum(pattern2*pattern2*weights,axis=1)
    

 
    
    inv_det = 1./(m11*m22-m12*m12)
    
    p1 = (m22*c1-m12*c2)*inv_det
    p2 = (m11*c2-m12*c1)*inv_det
    
    return p1,p2
    
def optimal_scaling(data,weights,pattern1):
    
    return np.sum(data*pattern1*weights) / np.sum(pattern1*pattern1*weights)
    
def chi_squared(valid,data,error,weight,model):

    nx = valid.shape[0]
    ny = valid.shape[1]

    chi2_array = np.zeros(data.shape,dtype=np.float32)

    code = """
    for (int i=0; i < nx; i++)
    {
        for (int j=0; j < ny; j++)
        {
            if (valid(i,j)==1 | valid(i,j)==4) {
                chi2_array(i,j) = ( data(i,j) - model(i,j) ) * ( data(i,j) - model(i,j) )* weight(i,j) ;
            } else if ((valid(i,j)==2 & data(i,j) > model(i,j)) | (valid(i,j)==3 & data(i,j) < model(i,j))) {
                if (error(i,j)==1.) {
                    chi2_array(i,j) = 1.e30;
                } else {
                    chi2_array(i,j) = -2. * log10(1-error(i,j));
                }
            }
        }
    }
    """

    weave.inline(code, ['valid','data','error','weight','model','nx','ny','chi2_array'], type_converters=converters.blitz)

    return np.sum(chi2_array,axis=1)
    
def read(directory,filters):
        
    model_fluxes = []
    wavelengths = []
    
    for filter in filters:
                
        filename = directory + '/convolved/' + filter + '.fits'
        if not os.path.exists(filename):
            if os.path.exists(filename+'.gz'):
                filename += '.gz'
            else:
                raise Exception("File not found: "+filename)
        
        print "Reading "+filename
        
        results = c.read_convolved_fluxes(filename)
        
        wavelengths.append(results[0])
        model_fluxes.append(results[3])
        
    model_fluxes = np.column_stack(model_fluxes)
    model_names = results[2]
    
    valid = model_fluxes<>0
    
    model_fluxes[~valid] = -np.inf
    model_fluxes[valid] = np.log10(model_fluxes[valid])
        
    return wavelengths,model_fluxes,model_names