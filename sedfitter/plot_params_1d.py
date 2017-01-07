from __future__ import print_function, division

from copy import deepcopy

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon

from .fit_info import FitInfoFile
from .models import load_parameter_table
from .utils import io
from .utils.formatter import LogFormatterMathtextAuto
from .plot_helpers import tex_friendly

__all__ = ['plot_params_1d']

fp = FontProperties(size='small')


def get_axes(fig):
    vxmin, vxmax = 1.5, 5.5
    vymin, vymax = 1.0, 4.0
    width, height = fig.get_figwidth(), fig.get_figheight()
    rect = [vxmin / width, vymin / width, (vxmax - vxmin) / width, (vymax - vymin) / height]
    return fig.add_axes(rect)


def plot_params_1d(input_fits, parameter, output_dir=None,
                   select_format=("N", 1), log_x=False, log_y=True,
                   label=None, bins=30, hist_range=None, additional={}, plot_name=True,
                   format='pdf', dpi=None):
    """
    Make histogram plots of parameters

    Parameters
    ----------
    input_fits : str or :class:`sedfitter.fit_info.FitInfo` or iterable
        This should be either a file containing the fit information, a
        :class:`sedfitter.fit_info.FitInfo` instance, or an iterable containing
        :class:`sedfitter.fit_info.FitInfo` instances.
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
    hist_range : tuple, optional
        The range of values to show
    additional : dict, optional
        A dictionary specifying additional parameters not listed in the
        parameter list for the models. Each item of the dictionary should
        itself be a dictionary giving the values for each model (where the key
        is the model name).
    plot_name: bool, optional
        Whether to show the source name on the plot(s).
    format: str, optional
        The file format to use for the plot, if output_dir is specified.
    dpi : int, optional
        The resolution of the figure to save
    """

    if output_dir is None:
        raise ValueError("No output directory has been specified")
    # Create output directory
    io.create_dir(output_dir)

    # Open input file
    fin = FitInfoFile(input_fits, 'r')

    # Read in table of parameters for model grid
    t = load_parameter_table(fin.meta.model_dir)

    # Sort alphabetically
    t['MODEL_NAME'] = np.char.strip(t['MODEL_NAME'])
    t.sort('MODEL_NAME')
    tpos = deepcopy(t)
    if log_x:
        tpos = tpos[tpos[parameter] > 0.]

    # Initialize figure
    fig = plt.figure()
    ax = get_axes(fig)

    # Find range of values
    if hist_range is None:
        pmin, pmax = tpos[parameter].min(), tpos[parameter].max()
    else:
        pmin, pmax = hist_range

    # Compute histogram
    if log_x:
        hist_all, edges = np.histogram(np.log10(tpos[parameter]), bins=bins, range=[np.log10(pmin), np.log10(pmax)])
        center = (edges[1:] + edges[:-1]) / 2.
        edges, center = 10. ** edges, 10. ** center
    else:
        hist_all, edges = np.histogram(t[parameter], bins=bins, range=[pmin, pmax])
        center = (edges[1:] + edges[:-1]) / 2.

    # Grayscale showing all models
    p = []
    for i in range(len(hist_all)):
        p.append((edges[i], max(hist_all[i], 0.01)))
        p.append((edges[i + 1], max(hist_all[i], 0.01)))
    p.append((edges[-1], 0.01))
    p.append((edges[0], 0.01))

    p = Polygon(p, facecolor='0.8', edgecolor='none')
    ax.add_patch(p)

    ax.set_xlabel(parameter if label is None else label)

    if log_x:
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(LogFormatterMathtextAuto())
    if log_y:
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(LogFormatterMathtextAuto())

    ax.set_xlim(pmin, pmax)
    ax.set_ylim(0.1, hist_all.max() * 10.)

    ax.set_autoscale_on(False)

    pfits = None
    source_label = None

    for info in fin:

        # Remove previous histogram
        if pfits is not None:
            pfits.remove()
        if source_label is not None:
            source_label.remove()

        # Filter fits
        info.keep(select_format)

        # Get filtered and sorted table of parameters
        tsorted = info.filter_table(t, additional=additional)

        # Compute histogram
        if log_x:
            hist, _ = np.histogram(np.log10(tsorted[parameter]), bins=bins, range=[np.log10(pmin), np.log10(pmax)])
        else:
            hist, _ = np.histogram(tsorted[parameter], bins=bins, range=[pmin, pmax])

        # Histogram showing values for good-fitting models
        p = []
        for i in range(len(hist)):
            p.append((edges[i], max(hist[i], 0.01)))
            p.append((edges[i + 1], max(hist[i], 0.01)))
        p.append((edges[-1], 0.01))
        p.append((edges[0], 0.01))
        pfits = Polygon(p, hatch='/', edgecolor='black', facecolor='none')
        ax.add_patch(pfits)

        if plot_name:
            source_label = ax.text(0.5, 0.95, tex_friendly(info.source.name),
                                   horizontalalignment='center',
                                   verticalalignment='center',
                                   transform=ax.transAxes,
                                   fontproperties=fp, zorder=200)

        # Save to file
        filename = "%s/%s.%s" % (output_dir, info.source.name, format)
        fig.savefig(filename, bbox_inches='tight', dpi=dpi)
        plt.close(fig)

    # Close input and output files
    fin.close()
