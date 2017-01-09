from __future__ import print_function, division

try:
    import cPickle as pickle
except ImportError:
    import pickle

from . import six

import numpy as np


class FitInfoFile(object):

    def __init__(self, fits, mode=None):

        if isinstance(fits, six.string_types):

            if mode not in 'wr':
                raise ValueError('mode should be r or w')

            self._handle = open(fits, mode + 'b')
            self._mode = mode

            if mode == 'r':
                self._first_meta = FitInfoMeta()
                self._first_meta.model_dir = pickle.load(self._handle)
                self._first_meta.filters = pickle.load(self._handle)
                self._first_meta.extinction_law = pickle.load(self._handle)
            else:
                self._first_meta = None

            self._fits = None

        elif isinstance(fits, FitInfo):

            self._fits = [fits]

        elif isinstance(fits, (list, tuple)):

            for info in self._fits[1:]:
                if info.meta != self._fits[0].meta:
                    raise ValueError("The meta property of all FitInfo instances should match")

            self._fits = fits

        else:

            raise TypeError('fits should be a string, FitInfo instance, or iterable of FitInfo instances')

    @property
    def meta(self):
        if self._fits is None:
            if self._mode != 'r':
                raise ValueError("meta property is only available in read mode")
            return self._first_meta
        else:
            return self._fits[0].meta

    def write(self, info):

        if self._mode != 'w':
            raise ValueError("File not open for writing")

        # We only write the metadata for the first source, and we then check
        # the metadata of other sources against the first one to make sure it
        # matches.

        if self._first_meta is None:
            pickle.dump(info.meta.model_dir, self._handle, 2)
            pickle.dump(info.meta.filters, self._handle, 2)
            pickle.dump(info.meta.extinction_law, self._handle, 2)
            self._first_meta = info.meta
        else:
            if not info.meta == self._first_meta:
                raise ValueError("meta does not match previously written value")

        pickle.dump(info, self._handle, 2)

    def close(self):
        if self._fits is None:
            self._handle.close()

    def __iter__(self):
        if self._fits is None:
            if self._mode != 'r':
                raise ValueError("File not open for reading")
            while True:
                try:
                    info = pickle.load(self._handle)
                except EOFError:
                    return
                else:
                    info.meta = self._first_meta
                    yield info
        else:
            for info in self._fits:
                yield info


class FitInfoMeta(object):
    def __eq__(self, other):
        return (self.model_dir == other.model_dir and
                self.filters == other.filters and
                self.extinction_law == other.extinction_law)


class FitInfo(object):
    """
    Results from a fit of a set of models to a source.
    """

    def __init__(self, source=None):

        self.source = source
        self.av = None
        self.sc = None
        self.chi2 = None
        self.model_id = None
        self.model_name = None
        self.model_fluxes = None
        self.meta = FitInfoMeta()

    def sort(self):
        """
        Sort the fit results from best to worst based on the chi^2 value.
        """

        order = np.argsort(self.chi2)
        self.av = self.av[order]
        self.sc = self.sc[order]
        self.chi2 = self.chi2[order]
        self.model_name = self.model_name[order]
        if self.model_fluxes is not None:
            self.model_fluxes = self.model_fluxes[order, :]
        self.model_id = order

    def keep(self, select_format):
        """
        Keep only a fraction of fits.

        Parameters
        ----------
        select_format : tuple, optional
            Tuple specifying which fits should be output. See the documentation
            for a description of the tuple syntax.
        """

        form, number = select_format

        if len(self.chi2) == 0:
            n_fits = 0
        elif form == 'A':
            n_fits = len(self.chi2)
        elif form == 'N' or not form:  # Form N is parsed as boolean
            n_fits = int(number)
        elif form == 'C':
            n_fits = np.sum(self.chi2 <= number)
        elif form == 'D':
            n_fits = np.sum(self.chi2 - self.chi2[0] <= number)
        elif form == 'E':
            n_fits = np.sum((self.chi2 / self.source.n_data) <= number)
        elif form == 'F':
            n_fits = np.sum((self.chi2 - self.chi2[0]) / self.source.n_data <= number)
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

    def filter_table(self, input_table, additional={}):
        """
        Given an input table, return only the rows matching the FitInfo object,
        and in the same order.

        Parameters
        ----------
        input_table : `~astropy.table.Table`
            The input table to filter and sort
        additional : dict, optional
            Additional columns to include in the table. This can be specified
            as a dictionary of Numpy arrays, where the key of the dictionary
            is the name of the column to add to the table.
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

        # Add additional parameter columns if necessary
        for par in additional:
            if par in table_sorted.columns:
                raise Exception("Parameter {0} already exists in table".format(par))
            table_sorted[par] = np.zeros(len(table_sorted), dtype=float)
            for i, name in enumerate(table_sorted['MODEL_NAME']):
                table_sorted[par][i] = additional[par][name.strip()]

        return table_sorted
