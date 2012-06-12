import atpy
import numpy as np


class OutputFile(object):

    def __init__(self):
        self.fitinfo = []

    def append(self, fitinfo):
        self.fitinfo.append(fitinfo)

    def write(self, par, filters, model_name):

        # Compute useful quantities
        n_wav = self.fitinfo[0].source.n_wav

        ts = atpy.TableSet()

        # Overall header

        ts.add_keyword('MODELDIR', par['modir'])
        ts.add_keyword('EXLAW', par['exlaw'])
        ts.add_keyword('OUT_FORM', par['oform'])
        ts.add_keyword('OUT_NUMB', par['onumb'])

        # Source table

        first = [1]
        for i in range(len(self.fitinfo) - 1):
            first.append(first[i] + self.fitinfo[i + 1].n_fits)

        number = [f.n_fits for f in self.fitinfo]

        ts.append(atpy.Table(name='SOURCES'))

        ts[0].add_keyword('NWAV', n_wav)

        for j in range(n_wav):
            ts[0].add_keyword("FILT%i" % (j + 1), filters[j]['name'])
            ts[0].add_keyword("WAV%i" % (j + 1), filters[j]['wav'])
            ts[0].add_keyword("AP%i" % (j + 1), filters[j]['ap'])

        ts[0].add_column('SOURCE_NAME', [f.source.name for f in self.fitinfo])
        ts[0].add_column('X', [f.source.x for f in self.fitinfo], unit='deg')
        ts[0].add_column('Y', [f.source.y for f in self.fitinfo], unit='deg')
        ts[0].add_column('SOURCE_ID', range(1, len(self.fitinfo) + 1))
        ts[0].add_column('FIRST_ROW', first)
        ts[0].add_column('NUMBER_FITS', number)
        ts[0].add_column('VALID', np.vstack([f.source.valid for f in self.fitinfo]))
        ts[0].add_column('FLUX', np.vstack([f.source.flux for f in self.fitinfo]), unit='mJy')
        ts[0].add_column('FLUX_ERROR', np.vstack([f.source.error for f in self.fitinfo]), unit='mJy')
        ts[0].add_column('LOG10FLUX', np.vstack([f.source.logflux for f in self.fitinfo]), unit='mJy')
        ts[0].add_column('LOG10FLUX_ERROR', np.vstack([f.source.logerror for f in self.fitinfo]), unit='mJy')

        # Fits table

        fit_id = []
        source_id = []

        for i, f in enumerate(self.fitinfo):
            source_id.append(np.ones(f.n_fits) * (i + 1))
            fit_id.append(range(1, f.n_fits + 1))
        source_id = np.hstack(source_id)
        fit_id = np.hstack(fit_id)

        model_id = np.hstack([f.order for f in self.fitinfo])

        model_name = model_name[model_id - 1]

        ts.append(atpy.Table(name='FITS'))

        ts[1].add_column('SOURCE_ID', source_id)
        ts[1].add_column('FIT_ID', fit_id)
        ts[1].add_column('MODEL_ID', model_id)
        ts[1].add_column('MODEL_NAME', model_name)
        ts[1].add_column('CHI2', np.hstack([f.chi2 for f in self.fitinfo]))
        ts[1].add_column('AV', np.hstack([f.av for f in self.fitinfo]), unit='mag')
        ts[1].add_column('SCALE', np.hstack([f.sc for f in self.fitinfo]))

        output_convolved = par['oconv'] in ['Y', 'y']
        if output_convolved:
            raise Exception("Not implemented")
            ts[1].add_keyword('MODELFLX', output_convolved)

        ts.write(par['ofile'])
