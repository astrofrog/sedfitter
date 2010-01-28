#!/usr/bin/env python

from distutils.core import setup

scripts = ['sed_fit', 'sed_plot', 'sed_filter_output', 'sed_fitinfo2data', 'sed_fitinfo2ascii']

setup(name='sedfitter',
      version='0.1.0',
      description='SED Fitter in python',
      author='Thomas Robitaille',
      author_email='trobitaille@cfa.harvard.edu',
      packages=['sedfitter'],
      scripts=['scripts/' + x for x in scripts]
     )
