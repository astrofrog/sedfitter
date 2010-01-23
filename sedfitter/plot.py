# Still to implement:
# - Axis numbers
# - Interpolation of SEDs (stype) - needs testing
# - Use number/format to know how many to plot (can use info.keep)

# Low priority
# - Grayscale SEDs
# - Color of lines
# - Option to plot/not plot SEDs (pseds)
# - Plot convolved fluxes (pconv)
# - Overplot stellar photosphere (patmo)
# - Color caption for multi-aperture plot
# - Ignore device (or translate to matplotlib)
# - Esthetic parameters (warn that they are being ignored for now)

# Optimization:
# - Compute interpolating functions as little as possible
# - Create LineCollection and PatchCollection and do the plotting at the very end

import parfile
import cPickle as pickle
from sed import SED

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mpl
from matplotlib.collections import LineCollection
from matplotlib.font_manager import FontProperties
from extinction import Extinction
from fit_info import FitInfo

import numpy as np

mpl.rc('axes', titlesize='small')
mpl.rc('axes', labelsize='small')
mpl.rc('xtick', labelsize='x-small')
mpl.rc('ytick', labelsize='x-small')
mpl.rc('font', family='serif')
mpl.rc('axes', linewidth=0.5)
mpl.rc('patch', linewidth=0.5)

fp = FontProperties(size='x-small')


def plot_source_info(ax, i, par, info, model_name):

    labels = []

    if par['pname']:
        labels.append(info.source.name)

    if par['pinfo']:
        labels.append("Model: %s" % model_name)
        if i==0:
            labels.append("Best fit")
        else:
            labels.append("Fit: %i" % (i+1))
        labels.append("$\chi^2$ = %10.3f    Av = %5.1f   Scale = %5.2f" % (info.chi2[i], info.av[i], info.sc[i]))

    pos = 0.95
    for label in labels:
        ax.text(0.5, pos, label, horizontalalignment='center',
                                  verticalalignment='center',
                                  transform = ax.transAxes,
                                  fontproperties=fp)
        pos -= 0.06

    return ax


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


def plot(parameter_file, input_file, output_dir):

    fin = file(input_file, 'rb')
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    model_names = pickle.load(fin)
    extinction_file = pickle.load(fin)
    extinction = Extinction(extinction_file)

    # Read in plotting parameters
    par = parfile.read(parameter_file, 'par')

    # Read in model parameters
    modpar = parfile.read("%s/models.conf" % model_dir, 'conf')

    info = FitInfo()

    while True:

        try:
            info.read(fin)
        except:
            break

        lines = []

        for i in range(info.n_fits):

            if (par['pmode'] == 'A' and i == 0) or par['pmode'] == 'I':
                fig = mpl.figure()
                ax = get_axes(fig)

            model_name = model_names[info.model_id[i]].strip()

            s = SED()
            if modpar['length_subdir'] == 0:
                s.read(model_dir + '/seds/' + model_name + '_sed.fits.gz')
            else:
                s.read(model_dir + '/seds/%s/%s_sed.fits.gz' % (model_name[:modpar['length_subdir']], model_name))

            s.scale_to_distance(10.**info.sc[i])
            s.scale_to_av(info.av[i], extinction.av)
            wav = np.array([f['wav'] for f in filters])
            ap = np.array([f['ap'] for f in filters])

            if par['stype'] == 'interp':
                apertures = ap * 10.**info.sc[i] * 1000.
                flux = s.interpolate_variable(wav, apertures)
            elif par['stype'] == 'largest':
                apertures = np.array([ap.max()]) * 10.**info.sc[i] * 1000.
                flux = s.interpolate(apertures)
            elif par['stype'] == 'largest+smallest':
                apertures = np.array([ap.min(), ap.max()]) * 10.**info.sc[i] * 1000.
                flux = s.interpolate(apertures)
            elif par['stype'] == 'all':
                raise Exception("stype=all not implemented")

            if flux.ndim > 1:
                for j in range(flux.shape[1]):
                    lines.append(np.column_stack([np.log10(s.wav),np.log10(flux[:,j])]))
            else:
                lines.append(np.column_stack([np.log10(s.wav),np.log10(flux)]))

            if (par['pmode'] == 'A' and i == info.n_fits-1) or par['pmode'] == 'I':

                ax.add_collection(LineCollection(lines, colors='0.75'))

                ax = plot_source_data(ax, info.source, filters)
                ax = set_view_limits(ax, par, wav, info.source)
                ax = plot_source_info(ax, i, par, info, model_name)

                ax.set_xlabel('$\lambda$ ($\mu$m)')
                ax.set_ylabel('$\lambda$F$_\lambda$ (ergs/cm$^2$/s)')

                if par['pmode'] == 'A':
                    fig.savefig("%s/%s.png" % (output_dir, info.source.name))
                else:
                    fig.savefig("%s/%s_%05i.png" % (output_dir, info.source.name, i+1))
