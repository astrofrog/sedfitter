import numpy as np
from scipy.interpolate import interp1d
import sys

def interpolate(filename,wavelengths):
    
     f = np.loadtxt(filename,dtype={'names':('wav','kl'),'formats':('f4','f4')},usecols=[0,3])
     
     exlaw = interp1d(f['wav'],f['kl'],bounds_error=False,fill_value=0)
          
     return -0.4 * exlaw(wavelengths) / exlaw(0.550000)
