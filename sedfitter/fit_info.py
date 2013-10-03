from __future__ import print_function, division

try:
    import cPickle as pickle
except ImportError:
    import pickle

import numpy as np

from .source import Source


class FitInfo(object):

    def __init__(self, source=None):

        self.source = source
        self.av = None
        self.sc = None
        self.chi2 = None
        self.model_id = None
        self.model_name = None
        self.model_fluxes = None

    def sort(self):

        order = np.argsort(self.chi2)
        self.av = self.av[order]
        self.sc = self.sc[order]
        self.chi2 = self.chi2[order]
        self.model_name = self.model_name[order]
        if self.model_fluxes is not None:
            self.model_fluxes = self.model_fluxes[order, :]
        self.model_id = order

    def keep(self, form, number):

        if form == 'A':
            n_fits = len(self.chi2)
        elif form == 'N' or not form:  # Form N is parsed as boolean
            n_fits = int(number)
        elif form == 'C':
            n_fits = np.sum(self.chi2 <= number)
        elif form == 'D':
            n_fits = np.sum(self.chi2 - self.chi2[0] <= number)
        elif form == 'E':
            n_fits = np.sum((self.chi2 / self.source.n_wav) <= number)
        elif form == 'F':
            n_fits = np.sum((self.chi2 - self.chi2[0]) / self.source.n_wav <= number)
        else:
            raise Exception("Unknown format: %s" % form)

        self.av = self.av[:n_fits]
        self.sc = self.sc[:n_fits]
        self.chi2 = self.chi2[:n_fits]
        self.model_name = self.model_name[:n_fits]
        if self.model_fluxes is not None:
            self.model_fluxes = self.model_fluxes[:n_fits]
        self.model_id = self.model_id[:n_fits]

    def __getattr__(self, attribute):
        if attribute == 'n_fits':
            return len(self.chi2)
        else:
            raise AttributeError(attribute)

    def __getstate__(self):
        return {
            'source': self.source,
            'av': self.av,
            'sc': self.sc,
            'chi2': self.chi2,
            'model_id': self.model_id,
            'model_name': self.model_name,
            'model_fluxes': self.model_fluxes,
        }

    def __setstate__(self, d):
        self.__init__()
        self.source = d['source']
        self.av = d['av']
        self.sc = d['sc']
        self.chi2 = d['chi2']
        self.model_id = d['model_id']
        self.model_name = d['model_name']
        self.model_fluxes = d['model_fluxes']

    def filter_table(self, input_table):
        """
        Given an input table, return only the rows matching the FitInfo object, and in the same order.
        """

        if not "MODEL_NAME" in input_table.dtype.names:
            raise ValueError("Input table should contain a MODEL_NAME column")

        subset = np.in1d(input_table['MODEL_NAME'], self.model_name)
        table_subset = input_table[subset]
        index = np.argsort(np.argsort(self.model_name))
        table_sorted = table_subset[index]

        # Double check that the sorting worked
        if not np.all(self.model_name == table_sorted['MODEL_NAME']):
            raise Exception("Parameter file sorting failed")

        return table_sorted
