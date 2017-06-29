from __future__ import print_function, division

import time

import numpy as np


class Timer(object):

    def __init__(self):
        self.time1 = time.time()
        self.n = 0
        self.step = 1
        print("   # Sources    CPU time (sec)    Sources/sec  ")
        print(" ----------------------------------------------")

    def display(self, force=False):
        self.n += 1
        if np.mod(self.n, self.step) == 0:
            self.time2 = time.time()
            if self.time2 - self.time1 < 1.:
                self.step *= 10
            else:
                print("    %7i       %10.1f        %7.2f" % (self.n, self.time2 - self.time1, self.n / (self.time2 - self.time1)))
        elif force:
            self.time2 = time.time()
            if self.time2 == self.time1:
                print("    %7i       %10.1f         -------" % (self.n, self.time2 - self.time1))
            else:
                print("    %7i       %10.1f        %7.2f" % (self.n, self.time2 - self.time1, self.n / (self.time2 - self.time1)))
