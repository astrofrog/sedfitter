#!/usr/bin/env python

from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, Command

from distutils.command.build_py import build_py

from sedfitter import __version__

class SEDFitterTest(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):

        import os
        import shutil
        import tempfile

        # First ensure that we build the package so that 2to3 gets executed
        self.reinitialize_command('build', inplace=False)
        self.run_command('build')
        build_cmd = self.get_finalized_command('build')
        new_path = os.path.abspath(build_cmd.build_lib)

        # Copy the build to a temporary directory for the purposes of testing
        # - this avoids creating pyc and __pycache__ directories inside the
        # build directory
        tmp_dir = tempfile.mkdtemp(prefix='sedfitter-test-')
        testing_path = os.path.join(tmp_dir, os.path.basename(new_path))
        shutil.copytree(new_path, testing_path)

        import sys
        import subprocess

        errno = subprocess.call([sys.executable, os.path.abspath('runtests.py')], cwd=testing_path)
        raise SystemExit(errno)

setup(name='sedfitter',
      version=__version__,
      description='SED Fitter in Python',
      author='Thomas Robitaille',
      author_email='thomas.robitaille@gmail.com',
      packages=['sedfitter',
                'sedfitter.convolve',
                'sedfitter.convolved_fluxes',
                'sedfitter.convolved_fluxes.tests',
                'sedfitter.extinction',
                'sedfitter.extinction.tests',
                'sedfitter.filter',
                'sedfitter.filter.tests',
                'sedfitter.sed',
                'sedfitter.sed.tests',
                'sedfitter.source',
                'sedfitter.source.tests',
                'sedfitter.tests',
                'sedfitter.utils',
                'sedfitter.utils.tests'],
      package_data={'sedfitter.sed.tests':['data/*.fits.gz'],
                    'sedfitter.filter.tests':['data/*.txt'],
                    'sedfitter.utils.tests':['data/*.conf', 'data/*.par']},
      provides=['sedfitter'],
      requires=['numpy', 'scipy', 'matplotlib', 'astropy'],
      cmdclass={'build_py': build_py, 'test':SEDFitterTest},
      keywords=['Scientific/Engineering'],
     )
