import numpy as np
import cPickle as pickle
from source import Source

class FitInfo(object):

    def __init__(self, source=None):

        self.source = source
        self.av = None
        self.sc = None
        self.chi2 = None
        self.model_id = None

        return

    def sort(self):

        order = np.argsort(self.chi2)
        self.av = self.av[order]
        self.sc = self.sc[order]
        self.chi2 = self.chi2[order]
        self.model_id = order

        return

    def keep(self, form, number):

        if form=='A':
            n_fits = len(self.chi2)
        elif form=='N' or not form: # Form N is parsed as boolean
            n_fits = int(number)
        elif form=='C':
            n_fits = len(self.chi2 <= number)
        elif form=='D':
            n_fits = len(self.chi2 - self.chi2[0] <= number)
        elif form=='E':
            n_fits = len((self.chi2 / self.source.n_wav) <= number)
        elif form=='F':
            n_fits = len((self.chi2 - self.chi2[0]) / self.source.n_wav <= number)
        else:
            raise Exception("Unknown format: %s" % form)

        self.av = self.av[:n_fits]
        self.sc = self.sc[:n_fits]
        self.chi2 = self.chi2[:n_fits]
        self.model_id = self.model_id[:n_fits]

        return

    def __getattr__(self, attribute):
        if attribute == 'n_fits':
            return len(self.model_id)
        else:
            raise AttributeError(attribute)

    def write(self, file_handle):
        self.source.write_binary(file_handle)
        pickle.dump(self.n_fits, file_handle)
        file_handle.write(self.av.astype(np.float32).tostring())
        file_handle.write(self.sc.astype(np.float32).tostring())
        file_handle.write(self.chi2.astype(np.float32).tostring())
        file_handle.write(self.model_id.astype(np.int32).tostring())

    def read(self, file_handle):
        self.source = Source()
        self.source.read_binary(file_handle)
        n_fits = pickle.load(file_handle)
        self.av = np.fromstring(file_handle.read(n_fits*4), dtype=np.float32)
        self.sc = np.fromstring(file_handle.read(n_fits*4), dtype=np.float32)
        self.chi2 = np.fromstring(file_handle.read(n_fits*4), dtype=np.float32)
        self.model_id = np.fromstring(file_handle.read(n_fits*4), dtype=np.int32)
