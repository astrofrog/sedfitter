# Error bars
# Interpolation of SEDs
# Color of lines
# Text on plots
# Grayscale SEDs

import parfile
import cPickle as pickle
from sed import SED

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mpl

from extinction import Extinction

import numpy as np


def plot_source_data(ax, source, filters):

    wav = np.array([f['wav'] for f in filters])
    plot_wav = np.log10(wav)
    plot_flux = source.logflux - 26. + np.log10(3.e8/(wav * 1.e-6))
    plot_error = source.logerror

    for j in range(source.n_wav):

        if source.valid[j] in [1, 4]:
            ax.scatter(plot_wav[j], plot_flux[j], marker='o', edgecolor='black', facecolor='black', zorder=100)
            ax.errorbar(plot_wav[j], plot_flux[j], yerr=plot_error[j], color='black', zorder=100)
        elif source.valid[j] == 2:
            ax.scatter(plot_wav[j], plot_flux[j], marker='^', edgecolor='black', facecolor='black', zorder=100)
        elif source.valid[j] == 3:
            ax.scatter(plot_wav[j], plot_flux[j], marker='v', edgecolor='black', facecolor='black', zorder=100)
        elif source.valid[j] == 9:
            ax.scatter(plot_wav[j], plot_flux[j], marker='o', edgecolor='black', facecolor='none', zorder=100)

    return ax


def set_view_limits(ax, par, wav, source):

    if par['xmode'] == 'A':
        xmin = np.log10(np.min(wav[source.valid<>0])) - par['xmina']
        xmax = np.log10(np.max(wav[source.valid<>0])) + par['xmaxa']
    else:
        xmin = np.log10(par['xminm'])
        xmax = np.log10(par['xmaxm'])

    if par['ymode'] == 'A':
        plot_flux = source.logflux - 26. + np.log10(3.e8/(wav * 1.e-6))
        ymin = np.min(plot_flux[source.valid<>0]) - par['ymina']
        ymax = np.max(plot_flux[source.valid<>0]) + par['ymaxa']
    else:
        ymin = np.log10(par['yminm'])
        ymax = np.log10(par['ymaxm'])

    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))

    return ax


def get_axes(fig):
    vxmin, vxmax = 1.5, 5.5
    vymin, vymax = 1.0, 4.0
    width, height = fig.get_figwidth(), fig.get_figheight()
    rect = [vxmin/width, vymin/width, (vxmax-vxmin)/width, (vymax-vymin)/height]
    return fig.add_axes(rect)

# call pgsci(1)
#
# do j=1, s%n_wav
#    select case(s%valid(j))
#    case(1, 4)
#       call pgslw(15)
#       call pgpt1(plot_wav(j), plot_flux(j), -1)
#       call pgslw(1)
#       call pgerr1(6, plot_wav(j), plot_flux(j), plot_error(j), 0.)
#    case(2)
#       call pgsch(1.20) ; call pgslw(3)
#       call pgpt1(plot_wav(j), plot_flux(j), 852)
#       call pgsch(0.75) ; call pgslw(1)
#    case(3)
#       call pgsch(1.20) ; call pgslw(3)
#       call pgpt1(plot_wav(j), plot_flux(j), 854)
#       call pgsch(0.75) ; call pgslw(1)
#    case(9)
#       call pgpt1(plot_wav(j), plot_flux(j), 4)
#    end select
# end do

def plot(parameter_file, input_file, output_dir):

    # Read in plotting parameters
    par = parfile.read(parameter_file)

    fin = file(input_file, 'rb')
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    model_names = pickle.load(fin)
    extinction_file = pickle.load(fin)
    extinction = Extinction(extinction_file)

    while True:

        try:
            info = pickle.load(fin)
        except:
            break

        for i in range(info.n_fits):

            s = SED()
            s.read(model_dir + '/seds/' + model_names[info.model_id[i]].strip() + '_sed.fits.gz')
            s.scale_to_distance(10.**info.sc[i])
            s.scale_to_av(info.av[i], extinction.av)
            wav = np.array([f['wav'] for f in filters])

            fig = mpl.figure()

            ax = get_axes(fig)

            ax.plot(np.log10(s.wav), np.log10(s.flux[0,:]))

            ax = plot_source_data(ax, info.source, filters)
            ax = set_view_limits(ax, par, wav, info.source)
            fig.savefig(output_dir + '/' + info.source.name)
