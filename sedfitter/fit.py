# Still to implement:
# - Min/max A_v
# - Minimum number of datapoints
# - Performance monitoring
# - YSO version

import numpy as np

import parfile
import models
import source
import extinction
import output
import timer
        
def fit(parfile):
        
    # Read in fitting parameters
    par = parfile.read(parfile)

    # Read in data
    sources = source.read_sources(par['dfile'])

    # Read in data format
    f = file(par['dform'],'rb')
    f.readline()
    nwav = int(f.readline().split('=')[0])
    filter_names = f.readline().split('=')[0].strip().split()
    apertures = np.array(f.readline().split('=')[0].strip().split(),float)

    # Read in models
    wavelengths,model_fluxes,model_names = models.read(par['modir'],filter_names)

    # Construct filters dictionary
    filters = []
    for i in range(nwav):
        filters.append({'wav':wavelengths[i],'ap':apertures[i],'name':filter_names[i]})

    # Set extinction model
    av_law = extinction.interpolate(par['exlaw'],wavelengths)

    # Set scale model
    sc_law = -2. * np.ones(av_law.shape)

    # Cycle through sources

    av_all = []
    sc_all = []
    chi_all = []
    order_all = []
    fit_id = []

    t = timer.Timer()

    for s in sources:
    
        if s.n_data > int(par['drequ']):

            av,sc,chi = models.fit(model_fluxes,s.valid,s.logflux,s.logerror,s.weight,av_law,sc_law)
    
            av,sc,chi,order = output.order_fits(av,sc,chi,form=par['oform'],number=par['onumb'],nwav=nwav)
    
            av_all.append(av)
            sc_all.append(sc)
            chi_all.append(chi)
            order_all.append(order+1)
        
            t.display()
                
    t.display(force=True)

    sys.exit(0)

    output.write(par,filters,sources,order_all,model_names,av_all,sc_all,chi_all)

