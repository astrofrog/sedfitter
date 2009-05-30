import io
import numpy as np

import pyfits

def order_fits(av,sc,chi,form=None,number=None,nwav=None):
    
    order = np.argsort(chi)
    av = av[order]
    sc = sc[order]
    chi = chi[order]
    
    if not form:
        pass
    elif form=='A':
        n_fits = len(chi)
    elif form=='N':
        n_fits = int(number)
    elif form=='C':
        n_fits = len(chi <= number)
    elif form=='D':
        n_fits = len(chi-chi[0] <= number)
    elif form=='E':
        n_fits = len(chi/float(nwav) <= number)
    elif form=='F':
        n_fits = len((chi-chi[0])/float(nwav) <= number)
        
    return av[:n_fits],sc[:n_fits],chi[:n_fits],order[:n_fits]

def write(par,filters,sources,model_id,model_name,av,sc,chi):

    fit_id = []
    source_id = []
    
    for i in range(len(chi)):
        size = len(chi[i])
        source_id.append(np.ones(size)*(i+1))
        fit_id.append(range(1,size+1))
        
    first = [1]
    for i in range(len(chi)-1):
        first.append(first[i] + len(chi[i+1]))
    first = np.array(first)
    
    number = np.array([len(x) for x in chi])
    
    av = np.hstack(av)
    sc = np.hstack(sc)
    chi = np.hstack(chi)
    model_id = np.hstack(model_id)
    fit_id = np.hstack(fit_id)
    source_id = np.hstack(source_id)

    model_name = model_name[model_id-1]
    
    io.delete_file(par['ofile'])

    hdu0 = pyfits.PrimaryHDU()
    hdu0.header.update('MODELDIR',par['modir'])
    hdu0.header.update('EXLAW',par['exlaw'])
    hdu0.header.update('OUT_FORM',par['oform'])
    hdu0.header.update('OUT_NUMB',par['onumb'])    
    
    nwav = len(sources[0].flux)
    
    cols = [ pyfits.Column(name='SOURCE_NAME',format='30A',array=np.array([sources[i].name for i in range(len(sources))])),
             pyfits.Column(name='X',format='1D',unit='deg',array=np.array([sources[i].x for i in range(len(sources))])),
             pyfits.Column(name='Y',format='1D',unit='deg',array=np.array([sources[i].y for i in range(len(sources))])),
             pyfits.Column(name='SOURCE_ID',format='1J',array=np.array(range(1,len(sources)+1))),
             pyfits.Column(name='FIRST_ROW',format='1J',array=first),
             pyfits.Column(name='NUMBER_FITS',format='1J',array=number),
             pyfits.Column(name='VALID',format=str(nwav)+"J",array=np.vstack([sources[i].valid for i in range(len(sources))])),
             pyfits.Column(name='FLUX',format=str(nwav)+"E",unit='mJy',array=np.vstack([sources[i].flux for i in range(len(sources))])),
             pyfits.Column(name='FLUX_ERROR',format=str(nwav)+"E",unit='mJy',array=np.vstack([sources[i].error for i in range(len(sources))])),
             pyfits.Column(name='LOG10FLUX',format=str(nwav)+"E",unit='mJy',array=np.vstack([sources[i].logflux for i in range(len(sources))])),
             pyfits.Column(name='LOG10FLUX_ERROR',format=str(nwav)+"E",unit='mJy',array=np.vstack([sources[i].logerror for i in range(len(sources))]))
             ]

    hdu1 = pyfits.new_table(cols)
    hdu1.name = "SOURCES"
    hdu1.header.update('NWAV',nwav)
    
    for j in range(nwav):
       hdu1.header.update("FILT"+str(j+1),filters[j]['name'],comment="Filter code")
       hdu1.header.update("WAV"+str(j+1),filters[j]['wav'],comment="Wavelength (microns)")
       hdu1.header.update("AP"+str(j+1),filters[j]['ap'],comment="Aperture (arcsec)")
    
    print chi.shape
    
    cols = [ pyfits.Column(name='SOURCE_ID',format='1J',array=source_id),
             pyfits.Column(name='FIT_ID',format='1J',array=fit_id),
             pyfits.Column(name='MODEL_ID',format='1J',array=model_id),
             pyfits.Column(name='MODEL_NAME',format='30A',array=model_name),
             pyfits.Column(name='CHI2',format='1E',array=chi),
             pyfits.Column(name='AV',format='1E',array=av,unit='mag'),
             pyfits.Column(name='SCALE',format='1E',array=sc) ]

    hdu2 = pyfits.new_table(cols)
    hdu2.name = "FITS"
    oconv = par['oconv'] in ['Y','y']
    hdu2.header.update('MODELFLX',oconv)    
    
    
    hdulist = pyfits.HDUList([hdu0,hdu1,hdu2])
    hdulist.writeto(par['ofile'])