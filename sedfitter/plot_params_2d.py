from __future__ import print_function, division

import os
import cPickle as pickle
from copy import deepcopy

import atpy
import numpy as np
from scipy.ndimage import convolve

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import LogFormatterMathtext

from .fit_info import FitInfo
from .extinction import Extinction
from . import util


# KERNEL = np.array([[ 0.  ,  0.39,  0.87,  1.  ,  0.87,  0.39,  0.  ],
#                    [ 0.39,  0.98,  1.  ,  1.  ,  1.  ,  0.98,  0.39],
#                    [ 0.87,  1.  ,  1.  ,  1.  ,  1.  ,  1.  ,  0.87],
#                    [ 1.  ,  1.  ,  1.  ,  1.  ,  1.  ,  1.  ,  1.  ],
#                    [ 0.87,  1.  ,  1.  ,  1.  ,  1.  ,  1.  ,  0.87],
#                    [ 0.39,  0.98,  1.  ,  1.  ,  1.  ,  0.98,  0.38],
#                    [ 0.  ,  0.39,  0.87,  1.  ,  0.87,  0.38,  0.01]])

KERNEL = np.array([[0., 0., 0.15, 0.58, 0.91, 1., 0.91, 0.58, 0.15, 0., 0.],
                   [0., 0.26, 0.98, 1., 1., 1., 1., 1., 0.98, 0.26, 0.],
                   [0.15, 0.98, 1., 1., 1., 1., 1., 1., 1., 0.98, 0.15],
                   [0.58, 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.58],
                   [0.91, 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.91],
                   [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
                   [0.91, 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.91],
                   [0.58, 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.58],
                   [0.15, 0.98, 1., 1., 1., 1., 1., 1., 1., 0.98, 0.15],
                   [0., 0.26, 0.98, 1., 1., 1., 1., 1., 0.98, 0.26, 0.],
                   [0., 0., 0.15, 0.58, 0.91, 1., 0.91, 0.58, 0.15, 0., 0.]])


class LogFormatterMathtextAuto(LogFormatterMathtext):

    def __call__(self, x, pos=None):
        if x in [0.001, 0.01, 0.1]:
            return str(x)
        elif x in [1., 10., 100., 1000.]:
            return str(int(x))
        else:
            return LogFormatterMathtext.__call__(self, x, pos=pos)


plt.rc('text', usetex=False)
plt.rc('axes', titlesize='small')
plt.rc('axes', labelsize='small')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('font', family='serif')
plt.rc('axes', linewidth=0.5)
plt.rc('patch', linewidth=0.5)

fp = FontProperties(size='small')


def get_axes(fig):
    vxmin, vxmax = 1.5, 4.5
    vymin, vymax = 1.0, 4.0
    width, height = fig.get_figwidth(), fig.get_figheight()
    rect = [vxmin / width, vymin / width, (vxmax - vxmin) / width, (vymax - vymin) / height]
    return fig.add_axes(rect)


def plot_params_2d(input_file, parameter_x, parameter_y, output_dir=None,
                   select_format=("N", 1), log_x=False, log_y=True,
                   label_x=None, label_y=None, additional={}, plot_name=True,
                   format='eps'):
    '''
    Make histogram plots of parameters

    Parameters
    ----------
    input_file : str
        File containing the fit information
    parameter_x : str
        The parameter to plot on the x-axis
    parameter_y : str
        The parameter to plot on the y-axis
    output_dir : str, optional
        If specified, plots are written to that directory
    select_format : tuple, optional
        Tuple specifying which fits should be plotted. See the documentation
        for a description of the tuple syntax.
    log_x : bool, optional
        Whether to plot the x-axis values in log space
    log_y : bool, optional
        Whether to plot the y-axis values in log space
    label_x : str, optional
        The x-axis label (if not specified, the parameter name is used)
    label_y : str, optional
        The y-axis label (if not specified, the parameter name is used)
    additional : dict, optional
        A dictionary specifying additional parameters not listed in the
        parameter list for the models. Each item of the dictionary should
        itself be a dictionary giving the values for each model (where the key
        is the model name).
    plot_name: bool, optional
        Whether to show the source name on the plot(s).
    format: str, optional
        The file format to use for the plot, if output_dir is specified.

    '''

    npix = 1024

    if output_dir is None:
        raise ValueError("No output directory has been specified")

    # Create output directory
    util.create_dir(output_dir)

    # Open output file
    fin = file(input_file, 'rb')

    # Read in header of output file
    model_dir = pickle.load(fin)
    filters = pickle.load(fin)
    extinction = Extinction()
    extinction.read_binary(fin)

    # Read in table of parameters for model grid
    if os.path.exists(model_dir + '/parameters.fits'):
        t = atpy.Table(model_dir + '/parameters.fits')
    elif os.path.exists(model_dir + '/parameters.fits.gz'):
        t = atpy.Table(model_dir + '/parameters.fits.gz')
    else:
        raise Exception("Parameter file not found in %s" % model_dir)

    # Sort alphabetically
    t['MODEL_NAME'] = np.char.strip(t['MODEL_NAME'])
    t.sort('MODEL_NAME')
    tpos = deepcopy(t)
    if log_x:
        tpos = tpos.where(tpos[parameter_x] > 0.)
    if log_y:
        tpos = tpos.where(tpos[parameter_y] > 0.)

    info = FitInfo()

    # Initialize figure
    fig = plt.figure()
    ax = get_axes(fig)

    # Find range of values
    xmin, xmax = tpos[parameter_x].min(), tpos[parameter_x].max()
    ymin, ymax = tpos[parameter_y].min(), tpos[parameter_y].max()

    # Compute histogram
    if log_x and log_y:
        gray_all, ex, ey = np.histogram2d(np.log10(tpos[parameter_x]),
                                          np.log10(tpos[parameter_y]), bins=npix,
                                          range=[[np.log10(xmin), np.log10(xmax)],
                                                 [np.log10(ymin), np.log10(ymax)]])
        ex, ey = 10. ** ex, 10. ** ey
    elif log_x:
        gray_all, ex, ey = np.histogram2d(np.log10(tpos[parameter_x]),
                                          tpos[parameter_y], bins=npix,
                                          range=[[np.log10(xmin), np.log10(xmax)],
                                                 [ymin, ymax]])
        ex = 10. ** ex
    elif log_y:
        gray_all, ex, ey = np.histogram2d(tpos[parameter_x],
                                          np.log10(tpos[parameter_y]), bins=npix,
                                          range=[[xmin, xmax],
                                                 [np.log10(ymin), np.log10(ymax)]])
        ey = 10. ** ey
    else:
        gray_all, ex, ey = np.histogram2d(tpos[parameter_x],
                                          tpos[parameter_y], bins=npix,
                                          range=[[xmin, xmax],
                                                [ymin, ymax]])

    gray_all = convolve(gray_all, KERNEL)
    gray_all = np.clip(gray_all, 0., 13.)

    # Grayscale showing all models
    ax.pcolormesh(ex, ey, gray_all.transpose(), cmap='binary', vmin=0, vmax=40.)

    ax.set_xlabel(parameter_x if label_x is None else label_x)
    ax.set_ylabel(parameter_y if label_y is None else label_y)

    if log_x:
        ax.xaxis.set_major_formatter(LogFormatterMathtextAuto())
        ax.set_xscale('log')
    if log_y:
        ax.yaxis.set_major_formatter(LogFormatterMathtextAuto())
        ax.set_yscale('log')

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.set_autoscale_on(False)

    pfits = None
    source_label = None

    while True:  # Loop over the fits

        # Read in next fit
        try:
            info.read_binary(fin)
        except:
            break

        # Remove previous histogram
        if pfits is not None:
            pfits.remove()
        if source_label is not None:
            source_label.remove()

        # Filter fits
        info.keep(select_format[0], select_format[1])

        # Match good-fitting models to parameter list
        subset = np.in1d(t['MODEL_NAME'], info.model_name)
        tsub = t.where(subset)
        index = np.argsort(np.argsort(info.model_name))
        tsorted = tsub.rows(index)
        if not np.all(info.model_name == tsorted['MODEL_NAME']):
            raise Exception("Parameter file sorting failed")

        # Add additional parameter columns if necessary
        for par in additional:
            if par in tsorted.columns:
                raise Exception("Parameter {} already exists in file".format(par))
            tsorted.add_empty_column(par, dtype=float)
            for i, name in enumerate(tsorted['MODEL_NAME']):
                tsorted[par][i] = additional[par][name.strip()]

        pfits = ax.scatter(tsorted[parameter_x], tsorted[parameter_y], c='black', s=10)

        if plot_name:
            source_label = ax.text(0.5, 0.95, info.source.name,
                                   horizontalalignment='center',
                                   verticalalignment='center',
                                   transform=ax.transAxes,
                                   fontproperties=fp, zorder=200)
        # Save to file
        filename = "%s/%s.%s" % (output_dir, info.source.name, format)
        fig.savefig(filename, bbox_inches='tight', facecolor='none', dpi=300)

    # Close input and output files
    fin.close()

    plt.close(fig)
