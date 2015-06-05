import matplotlib.pyplot as plt


def set_rc_params():

    plt.rc('axes', titlesize='small')
    plt.rc('axes', labelsize='small')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.rc('font', family='serif')
    plt.rc('axes', linewidth=0.5)
    plt.rc('patch', linewidth=0.5)


def tex_friendly(string):
    if plt.rcParams['text.usetex']:
        return string.replace('_', '\_').replace('%', '\%')
    else:
        return string