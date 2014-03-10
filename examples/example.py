import os
import glob

from astropy import units as u

from sedfitter import (fit, plot, plot_params_1d, plot_params_2d,
                       write_parameters, write_parameter_ranges)
from sedfitter.extinction import Extinction

# Define path to models
model_dir = '/Volumes/Data/models/models_r06'

# Read in extinction law. We read in columns 1 (the wavelength in microns) and
# 4 (the opacity in cm^2/g)
extinction = Extinction.from_file('kmh94.par', columns=[0, 3],
                                  wav_unit=u.micron, chi_unit=u.cm**2 / u.g)

# Define filters and apertures
filters = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']
apertures = [3., 3., 3., 3., 3., 3., 3.] * u.arcsec

# Run the fitting
fit('data_glimpse', filters, apertures, model_dir,
    'output.fitinfo',
    extinction_law=extinction,
    distance_range=[1., 2.] * u.kpc,
    av_range=[0., 40.])

# For the remaining commands, we always select the models with chi^2-chi_best^2
# per datapoint less than 3.
select_format = ('F', 3)

# Make SED plots
plot('output.fitinfo', 'plots_seds', plot_max=100, select_format=select_format)

# Make histograms of the disk mass
plot_params_1d('output.fitinfo', 'MDISK', 'plots_mdisk',
               log_x=True, select_format=select_format)

# Make 2-d plots of the envelope infall rate vs disk mass
plot_params_2d('output.fitinfo', 'MDISK', 'MDOT', 'plots_mdot_mdisk',
               log_x=True, log_y=True, select_format=select_format)

# Write out all models with a delta chi^2-chi_best^2 per datapoint < 3
write_parameters('output.fitinfo', 'parameters.txt',
                 select_format=select_format)

# Write out the min/max ranges corresponding to the above file
write_parameter_ranges('output.fitinfo', 'parameter_ranges.txt',
                       select_format=select_format)