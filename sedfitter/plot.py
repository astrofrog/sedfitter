# Axis numbers
# Interpolation of SEDs
# Color of lines
# Grayscale SEDs

import parfile
import cPickle as pickle
from sed import SED

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mpl
from matplotlib.font_manager import FontProperties
from extinction import Extinction

import numpy as np


mpl.rc('axes', titlesize='small')
mpl.rc('axes', labelsize='small')
mpl.rc('xtick', labelsize='x-small')
mpl.rc('ytick', labelsize='x-small')
mpl.rc('font', family='serif')
mpl.rc('axes', linewidth=0.5)
mpl.rc('patch', linewidth=0.5)

fp = FontProperties(size='x-small')


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

            model_name = model_names[info.model_id[i]].strip()

            s = SED()
            s.read(model_dir + '/seds/' + model_name + '_sed.fits.gz')
            s.scale_to_distance(10.**info.sc[i])
            s.scale_to_av(info.av[i], extinction.av)
            wav = np.array([f['wav'] for f in filters])

            fig = mpl.figure()

            ax = get_axes(fig)

            ax.plot(np.log10(s.wav), np.log10(s.flux[0,:]))

            ax = plot_source_data(ax, info.source, filters)
            ax = set_view_limits(ax, par, wav, info.source)

            labels = []

            if par['pname'].lower() == 'y':
                labels.append(info.source.name)

            if par['pinfo'].lower() == 'y':
                labels.append("Model: %s" % model_name)
                if i==0:
                    labels.append("Best fit")
                else:
                    labels.append("Fit: %i" % fit_id)
                labels.append("$\chi^2$ = %10.3f    Av = %5.1f   Scale = %5.2f" % (info.chi2[i], info.av[i], info.sc[i]))

            pos = 0.95
            for label in labels:
                ax.text(0.5, pos, label, horizontalalignment='center',
                                          verticalalignment='center',
                                          transform = ax.transAxes,
                                          fontproperties=fp)
                pos -= 0.06

            fig.savefig(output_dir + '/' + info.source.name)
