import atpy


class ConvolvedFluxes(object):

    def __init__(self, *args):
        if args:
            self.read(*args)

    def read(self, filename):
        '''
        Read convolved fluxes from a FITS file
        '''

        t = atpy.TableSet(filename, verbose=False)

        if 'FILTWAV' in t.keywords:
            self.wavelength = t.keywords['FILTWAV']
        else:
            self.wavelength = None

        self.n_models = t.keywords['NMODELS']
        self.n_ap = t.keywords['NAP']

        self.model_names = t[0].MODEL_NAME
        self.fluxes = t[0].TOTAL_FLUX
        self.flux_errors = t[0].TOTAL_FLUX_ERR
        self.radius_sigma_50 = t[0].RADIUS_SIGMA_50
        self.radius_cumul_99 = t[0].RADIUS_CUMUL_99

        self.apertures = t[1].APERTURE

        return

    def write(self, filename):
        '''
        Write convolved fluxes to a FITS file
        '''

        ts = atpy.TableSet()

        ts.add_keyword('FILTWAV', self.wavelength)
        ts.add_keyword('NMODELS', self.n_models)
        ts.add_keyword('NAP', self.n_ap)

        ts.append(atpy.Table(name='CONVOLVED FLUXES'))
        ts[0].add_column('MODEL_NAME', self.model_names)
        ts[0].add_column('TOTAL_FLUX', self.fluxes)
        ts[0].add_column('TOTAL_FLUX_ERR', self.fluxes_err)
        ts[0].add_column('RADIUS_SIGMA_50', self.radius_sigma_50)
        ts[0].add_column('RADIUS_CUMUL_99', self.radius_cumul_50)

        ts.append(atpy.Table(name='APERTURES'))
        ts[1].add_column("APERTURE", self.apertures)

        ts.write(filename, verbose=False)
