import os
import glob

import matplotlib
matplotlib.use('Agg')

import sedfitter
from sedfitter.extinction import Extinction

model_dir = '/Users/tom/Models/models_r06'

# Read in extinction law

extinction = Extinction.from_file('kmh94.par')

# Define filters and apertures

filters = ['2J', '2H', '2K', 'I1', 'I2', 'I3', 'I4']

apertures = [3., 3., 3., 3., 3., 3., 3.]

# Run the fitting

sedfitter.fit('data_glimpse', filters, apertures, model_dir,
              'output_glimpse.fitinfo',
              extinction_law=extinction,
              distance_range=[1., 2.],
              av_range=[0., 40.],
              output_format=('F', 3.),
              output_convolved=False,
              remove_resolved=True)

sedfitter.plot('output_glimpse.fitinfo',
               plot_mode='A',
               output_dir='plots_glimpse',
               select_format=('F', 2.),
               format='png',
               show_convolved=False, show_sed=True,
               x_mode='M', x_range=(0.1, 2000),
               y_mode='M', y_range=(1.e-14, 2e-8),
               plot_max=100)

sedfitter.plot_params_1d('output_glimpse.fitinfo', 'TSTAR',
                         log_x=True,
                         output_dir='plots_glimpse_1d',
                         select_format=('F', 2.), format='png')

sedfitter.plot_params_2d('output_glimpse.fitinfo', 'RSTAR', 'TSTAR',
                         log_x=True, log_y=True,
                         output_dir='plots_glimpse_2d',
                         select_format=('F', 2.), format='png')
