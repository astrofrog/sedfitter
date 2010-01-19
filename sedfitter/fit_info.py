import numpy as np


class FitInfo(object):

    def __init__(self, source):

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
        elif form=='N':
            n_fits = int(number)
        elif form=='C':
            n_fits = len(self.chi2 <= number)
        elif form=='D':
            n_fits = len(self.chi2 - self.chi2[0] <= number)
        elif form=='E':
            n_fits = len((self.chi2 / self.source.n_wav) <= number)
        elif form=='F':
            n_fits = len((self.chi2 - self.chi2[0]) / self.source.n_wav <= number)

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
