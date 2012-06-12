#!/usr/bin/env python

import os
import cPickle as pickle

import numpy as np

import atpy

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon

from fit_info import FitInfo
from extinction import Extinction

from matplotlib.ticker import LogFormatterMathtext

import util

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

fp = FontProperties(size='x-small')


def get_axes(fig):
    vxmin, vxmax = 1.5, 5.5
    vymin, vymax = 1.0, 4.0
    width, height = fig.get_figwidth(), fig.get_figheight()
    rect = [vxmin / width, vymin / width, (vxmax - vxmin) / width, (vymax - vymin) / height]
    return fig.add_axes(rect)


def plot_params_1d(input_file, parameter, output_dir=None,
                   select_format=("N", 1), log_x=False, log_y=True,
                   label=None, bins=30, additional={}, plot_name=True,
                   format='eps'):
    '''
    Make histogram plots of parameters

    Parameters
    ----------
    input_file : str
        File containing the fit information
    parameter : str
        The parameter to plot a histogram of
    output_dir : str, optional
        If specified, plots are written to that directory
    select_format : tuple, optional
        Tuple specifying which fits should be plotted. See the documentation
        for a description of the tuple syntax.
    log_x : bool, optional
        Whether to plot the x-axis values in log space
    log_y : bool, optional
        Whether to plot the y-axis values in log space
    label : str, optional
        The x-axis label (if not specified, the parameter name is used)
    bins : int, optional
        The number of bins for the histogram
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

    info = FitInfo()

    # Initialize figure
    fig = plt.figure()
    ax = get_axes(fig)

    # Find range of values
    pmin, pmax = t[parameter].min(), t[parameter].max()

    # Compute histogram
    if log_x:
        hist_all, edges = np.histogram(np.log10(t[parameter]), bins=bins, range=[np.log10(pmin), np.log10(pmax)])
        center = (edges[1:] + edges[:-1]) / 2.
        edges, center = 10. ** edges, 10. ** center
    else:
        hist_all, edges = np.histogram(t[parameter], bins=bins, range=[pmin, pmax])
        center = (edges[1:] + edges[:-1]) / 2.

    # Grayscale showing all models
    p = []
    for i in range(len(hist_all)):
        p.append((edges[i], max(hist_all[i], 0.01)))
        p.append((edges[i+1], max(hist_all[i], 0.01)))
    p.append((edges[-1], 0.01))
    p.append((edges[0], 0.01))

    p = Polygon(p, facecolor='0.8', edgecolor='none')
    ax.add_patch(p)

    ax.xaxis.set_major_formatter(LogFormatterMathtextAuto())
    ax.yaxis.set_major_formatter(LogFormatterMathtextAuto())
    ax.set_xlabel(parameter if label is None else label)

    if log_x:
        ax.set_xscale('log')
    if log_y:
        ax.set_yscale('log')

    ax.set_xlim(pmin, pmax)
    ax.set_ylim(0.1, hist_all.max() * 10.)

    ax.set_autoscale_on(False)

    pfits = None

    while True:  # Loop over the fits

        # Read in next fit
        try:
            info.read_binary(fin)
        except:
            break

        # Remove previous histogram
        if pfits is not None:
            pfits.remove()

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

        # Compute histogram
        if log_x:
            hist, _ = np.histogram(np.log10(tsorted[parameter]), bins=bins, range=[np.log10(pmin), np.log10(pmax)])
        else:
            hist, _ = np.histogram(tsorted[parameter], bins=bins, range=[pmin, pmax])

        # Histogram showing values for good-fitting models
        p = []
        for i in range(len(hist)):
            p.append((edges[i], max(hist[i], 0.01)))
            p.append((edges[i+1], max(hist[i], 0.01)))
        p.append((edges[-1], 0.01))
        p.append((edges[0], 0.01))
        pfits = Polygon(p, hatch='/', edgecolor='black', facecolor='none')
        ax.add_patch(pfits)

        # Save to file
        filename = "%s/%s.%s" % (output_dir, info.source.name, format)
        fig.savefig(filename, bbox_inches='tight')
        plt.close(fig)

    # Close input and output files
    fin.close()
