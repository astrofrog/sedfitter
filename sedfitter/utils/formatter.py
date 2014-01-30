from matplotlib.ticker import LogFormatterMathtext


class LogFormatterMathtextAuto(LogFormatterMathtext):

    def __call__(self, x, pos=None):
        if x in [0.001, 0.01, 0.1]:
            return str(x)
        elif x in [1., 10., 100., 1000.]:
            return str(int(x))
        else:
            return LogFormatterMathtext.__call__(self, x, pos=pos)
