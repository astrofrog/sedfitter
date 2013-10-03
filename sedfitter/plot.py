from __future__ import print_function, division

# Still to implement:
# - Axis numbers - done
# - Interpolation of SEDs (stype) - done - needs testing
# - Use number/format to know how many to plot (can use info.keep) - done - needs testing

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

try:
    import cPickle as pickle
except ImportError:
    import pickle

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker

from .extinction import Extinction
from .fit_info import FitInfo
from .sed import SED
from .utils import io
from .utils import parfile
from .utils.formatter import LogFormatterMathtextAuto

KPC = 3.086e21

plt.rc('text', usetex=True)
plt.rc('axes', titlesize='small')
plt.rc('axes', labelsize='small')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('font', family='serif')
plt.rc('axes', linewidth=0.5)
plt.rc('patch', linewidth=0.5)

fp = FontProperties(size='x-small')

color = {}

color['gray'] = '0.75'
color['black'] = '0.00'

color['full'] = []
color['full'].append((0.65, 0.00, 0.00))
color['full'].append((0.20, 0.30, 0.80))
color['full'].append((1.00, 0.30, 0.30))
color['full'].append((1.00, 0.30, 0.60))
color['full'].append((0.30, 0.80, 0.30))
color['full'].append((0.50, 0.10, 0.80))
color['full'].append((0.20, 0.60, 0.80))
color['full'].append((1.00, 0.00, 0.00))
color['full'].append((0.50, 0.25, 0.00))
color['full'].append((0.90, 0.90, 0.00))
color['full'].append((0.00, 0.50, 0.00))

color['faded'] = []
color['faded'].append((1.00, 0.70, 0.70))
color['faded'].append((0.70, 0.70, 0.80))
color['faded'].append((1.00, 0.80, 0.70))
color['faded'].append((1.00, 0.75, 1.00))
color['faded'].append((0.70, 0.80, 0.70))
color['faded'].append((0.75, 0.60, 0.80))
color['faded'].append((0.70, 0.75, 0.80))
color['faded'].append((1.00, 0.70, 0.80))
color['faded'].append((0.90, 0.80, 0.70))
color['faded'].append((0.90, 0.90, 0.70))
color['faded'].append((0.50, 0.90, 0.50))


def tex_friendly(string):
    return string
    # return string.replace('_', '\_').replace('%', '\%')


def plot_source_info(ax, i, info, plot_name, plot_info):

    labels = []

    if plot_name:
        labels.append(tex_friendly(info.source.name))

    if plot_info:
        labels.append("Model: %s" % tex_friendly(info.model_name[i]))
        if i == 0:
            labels.append("Best fit")
        else:
            labels.append("Fit: %i" % (i + 1))
        labels.append("$\chi^2$ = %10.3f    A$_{\\rm V}$ = %5.1f    Scale = %5.2f" % (info.chi2[i], info.av[i], info.sc[i]))

    pos = 0.95
    for label in labels:
        ax.text(0.5, pos, label, horizontalalignment='center',
                verticalalignment='center',
                                  transform=ax.transAxes,
                                  fontproperties=fp)
        pos -= 0.06

    return ax


def plot_source_data(ax, source, filters, size=20, capsize=3):

    weight, log_flux, log_error = source.get_log_fluxes()

    wav = np.array([f['wav'] for f in filters])
    plot_wav = wav
    plot_flux = log_flux - 26. + np.log10(3.e8 / (wav * 1.e-6))
    plot_flux_up = 10. ** (plot_flux + log_error)
    plot_flux_down = 10. ** (plot_flux - log_error)
    plot_flux = 10. ** plot_flux
    plot_error = np.vstack([plot_flux - plot_flux_down, plot_flux_up - plot_flux])

    for j in range(source.n_wav):

        if source.valid[j] in [1, 4]:
            ax.scatter(plot_wav[j], plot_flux[j], marker='o', edgecolor='black', facecolor='black', zorder=100, s=size)
            ax.errorbar(plot_wav[j], plot_flux[j], yerr=plot_error[:, j:j + 1], color='black', zorder=100, capsize=capsize)
        elif source.valid[j] == 2:
            ax.scatter(plot_wav[j], plot_flux[j], marker='^', edgecolor='black', facecolor='black', zorder=100, s=size)
        elif source.valid[j] == 3:
            ax.scatter(plot_wav[j], plot_flux[j], marker='v', edgecolor='black', facecolor='black', zorder=100, s=size)
        elif source.valid[j] == 9:
            ax.scatter(plot_wav[j], plot_flux[j], marker='o', edgecolor='black', facecolor='none', zorder=100, s=size)

    return ax


def set_view_limits(ax, wav, source, x_mode, y_mode, x_range, y_range):

    if x_mode == 'A':
        xmin = np.min(wav[source.valid != 0]) * 10. ** (-x_range[0])
        xmax = np.max(wav[source.valid != 0]) * 10. ** (+x_range[1])
    else:
        xmin = x_range[0]
        xmax = x_range[1]

    if y_mode == 'A':
        weight, log_flux, log_error = source.get_log_fluxes()
        plot_flux = 10. ** (log_flux - 26. + np.log10(3.e8 / (wav * 1.e-6)))
        ymin = np.min(plot_flux[source.valid != 0]) * 10. ** (-y_range[0])
        ymax = np.max(plot_flux[source.valid != 0]) * 10. ** (+y_range[1])
    else:
        ymin = y_range[0]
        ymax = y_range[1]

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))
    ax.xaxis.set_major_formatter(LogFormatterMathtextAuto())


    return ax


def get_axes(fig):
    vxmin, vxmax = 1.5, 5.5
    vymin, vymax = 1.0, 4.0
    width, height = fig.get_figwidth(), fig.get_figheight()
    rect = [vxmin / width, vymin / width, (vxmax - vxmin) / width, (vymax - vymin) / height]
    return fig.add_axes(rect)


def plot(input_file, output_dir=None, select_format=("N", 1), plot_max=None,
         plot_mode="A", sed_type="interp", show_sed=True, show_convolved=False,
         x_mode='A', y_mode='A', x_range=(1., 1.), y_range=(1., 2.),
         plot_name=True, plot_info=True, format='eps'):
    '''
    Make SED plots

    Parameters
    ----------
    input_file : str
        File containing the fit information
    output_dir : str, optional
        If specified, plots are written to that directory
    select_format : tuple, optional
        Tuple specifying which fits should be plotted. See the documentation
        for a description of the tuple syntax.
    plot_max : int, optional
        Maximum number of fits to plot
    plot_mode : str, optional
        Whether to plot all fits in a single plot ('A') or one fit per plot
        ('I')
    sed_type : str, optional
        Which SED should be shown:
            * `largest`: show the SED for the largest aperture specified.
            * `smallest`: show the SED for the smallest aperture specified.
            * `largest+smallest`: show the SEDs for the largest and smallest
              apertures specified.
            * `all`: show the SEDs for all apertures specified.
            * `interp`: interpolate the SEDs to the correct aperture at each
              wavelength (and interpolated apertures in between), so that a
              single composite SED is shown.
    show_sed : bool, optional
        Show the SEDs
    show_convolved : bool, optional
        Show convolved model fluxes
    x_mode : str, optional
        Whether to automatically select the wavelength range ('A'), or whether
        to use manually set values ('M').
    x_range : tuple, optional
        If x_mode is set to 'M', this is the range of wavelengths to show. If
        x_mode is set to 'A', this is the marging to add around the wavelength
        range (in dex).
    y_mode : str, optional
        Whether to automatically select the flux range ('A'), or whether to
        use manually set values ('M').
    y_range : tuple, optional
        If y_mode is set to 'M', this is the range of fluxes to show. If
        y_mode is set to 'A', this is the marging to add around the flux
        range (in dex).
    plot_name : bool, optional
        Whether to show the source name on the plot(s).
    plot_info : bool, optional
        Whether to show the fit information on the plot(s).
    format : str, optional
        The file format to use for the plot, if output_dir is specified.
    '''

    if output_dir:
        io.create_dir(output_dir)
    else:
        figures = {}

    fin = open(input_file, 'rb')

    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction = pickle.load(fin)

    wav = np.array([f['wav'] for f in filters])
    ap = np.array([f['aperture_arcsec'] for f in filters])

    unique_ap = np.unique(ap)
    unique_ap.sort()

    # Read in model parameters
    modpar = parfile.read("%s/models.conf" % model_dir, 'conf')

    while True:

        # Read in next fit
        try:
            info = pickle.load(fin)
        except EOFError:
            break

        # Filter fits
        info.keep(select_format[0], select_format[1])

        if plot_max:
            info.keep('N', plot_max)

        if show_convolved and info.model_fluxes is None:
            raise Exception("Cannot plot convolved fluxes as these are not included in the input file")

        for i in range(info.n_fits - 1, -1, -1):

            if (plot_mode == 'A' and i == info.n_fits - 1) or plot_mode == 'I':
                fig = plt.figure()
                ax = get_axes(fig)

            # Initalize lines and colors list
            if (plot_mode == 'A' and i == info.n_fits - 1) or plot_mode == 'I':
                lines = []
                colors = []
                if show_convolved:
                    conv = []

            if (plot_mode == 'A' and i == 0) or plot_mode == 'I':
                if sed_type in ['interp', 'largest']:
                    color_type = 'black'
                else:
                    color_type = 'full'
            else:
                if sed_type in ['interp', 'largest']:
                    color_type = 'gray'
                else:
                    color_type = 'faded'

            if modpar['length_subdir'] == 0:
                s = SED.read(model_dir + '/seds/' + info.model_name[i] + '_sed.fits')
            else:
                s = SED.read(model_dir + '/seds/%s/%s_sed.fits' % (info.model_name[i][:modpar['length_subdir']], info.model_name[i]))

            s = s.scale_to_distance(10. ** info.sc[i] * KPC)
            s = s.scale_to_av(info.av[i], extinction.get_av)

            if sed_type == 'interp':
                apertures = ap * 10. ** info.sc[i] * 1000.
                flux = s.interpolate_variable(wav, apertures)
            elif sed_type == 'largest':
                apertures = np.array([ap.max()]) * 10. ** info.sc[i] * 1000.
                flux = s.interpolate(apertures)
            elif sed_type == 'largest+smallest':
                apertures = np.array([ap.min(), ap.max()]) * 10. ** info.sc[i] * 1000.
                flux = s.interpolate(apertures)
            elif sed_type == 'all':
                apertures = unique_ap * 10. ** info.sc[i] * 1000.
                flux = s.interpolate(apertures)

            if flux.ndim > 1:
                for j in range(flux.shape[1]):
                    lines.append(np.column_stack([s.wav, flux[:, j]]))
                    colors.append(color[color_type][j])
            else:
                lines.append(np.column_stack([s.wav, flux]))
                colors.append(color[color_type])

            if show_convolved:
                conv.append(10. ** (info.model_fluxes[i, :] - 26. + np.log10(3.e8 / (wav * 1.e-6))))

            if (plot_mode == 'A' and i == 0) or plot_mode == 'I':

                if show_sed:
                    ax.add_collection(LineCollection(lines, colors=colors))

                if show_convolved:
                    for j in range(len(conv)):
                        ax.plot(wav, conv[j], color=colors[j], linestyle='solid', marker='o', markerfacecolor='none', markeredgecolor=colors[j])

                ax = plot_source_data(ax, info.source, filters)

                if plot_mode == 'A':
                    ax = plot_source_info(ax, 0, info, plot_name, plot_info)
                else:
                    ax = plot_source_info(ax, i, info, plot_name, plot_info)

                ax.set_xlabel('$\lambda$ ($\mu$m)')
                ax.set_ylabel('$\lambda$F$_\lambda$ (ergs/cm$^2$/s)')

                ax = set_view_limits(ax, wav, info.source, x_mode, y_mode, x_range, y_range)

                if output_dir:
                    if plot_mode == 'A':
                        filename = "%s/%s.%s" % (output_dir, info.source.name, format)
                    else:
                        filename = "%s/%s_%05i.%s" % (output_dir, info.source.name, i + 1, format)
                    fig.savefig(filename, bbox_inches='tight')
                    plt.close(fig)
                else:
                    figures[info.source.name] = {'source': info.source, 'filters': filters, 'lines': LineCollection(lines, colors=colors)}

    if not output_dir:
        return figures
