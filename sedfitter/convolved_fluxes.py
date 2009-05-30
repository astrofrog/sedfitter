import pyfits
import io

def read_convolved_fluxes(filename):
    
    hdu0,hdu1,hdu2 = pyfits.open(filename)

    if hdu0.header.has_key('FILTWAV'):
        wavelength = hdu0.header['FILTWAV']
    else:
        wavelength = None
        
    n_models = hdu0.header['NMODELS']
    n_wav = hdu0.header['NAP']
    
    model_names = hdu1.data.field('MODEL_NAME')
    fluxes = hdu1.data.field('TOTAL_FLUX')
    flux_errors = hdu1.data.field('TOTAL_FLUX_ERR')
    radius_sigma_50 = hdu1.data.field('RADIUS_SIGMA_50')
    radius_cumul_99 = hdu1.data.field('RADIUS_CUMUL_99')
    
    apertures = hdu2.data.field('APERTURE')
    
    return (wavelength,apertures,model_names,fluxes,flux_errors,radius_sigma_50,radius_cumul_99)
    
#def interpolate_convolved_fluxes(filename,apertures):
    
#    models = read_convolved_fluxes(filename)
    
#    resolved = np.zeros(fluxes.shape,bool)
#    extended = np.zeros(fluxes.shape,bool)
    
    
    
#    for aperture in apertures:
        
    
def write(filename,wavelength,model_names,fluxes,fluxes_err,apertures):

	io.delete_file(filename)

	hdu0 = pyfits.PrimaryHDU()
	hdu0.header.update('FILTWAV',wavelength)
	hdu0.header.update('NMODELS',len(fluxes))
	hdu0.header.update('NAP',len(apertures))
	
	c1 = pyfits.Column(name='MODEL_NAME',format='30A',array=model_names)
	c2 = pyfits.Column(name='TOTAL_FLUX',format='1E',array=fluxes)
	c3 = pyfits.Column(name='TOTAL_FLUX_ERR',format='1E',array=fluxes_err)
#	c4 = pyfits.Column(name='RADIUS_SIGMA_50',format='E',array=)
#	c5 = pyfits.Column(name='RADIUS_CUMUL_99',format='E',array=)
	hdu1 = pyfits.new_table([c1,c2,c3])
	
	c1 = pyfits.Column(name='APERTURE',format='1E',array=apertures)
	hdu2 = pyfits.new_table([c1])
	
	hdulist = pyfits.HDUList([hdu0,hdu1,hdu2])
	hdulist.writeto(filename)
	