About
-----

This is a Python port of the SED fitter from [Robitaille et al., 2007, ApJS 169
328](http://adsabs.harvard.edu/abs/2007ApJS..169..328R). It replaces the
original Fortran code which can be found
[here](https://github.com/astrofrog/sedfitter-legacy).

At this stage, it should be used cautiously as it is still beta-quality
software and may contain bugs. Please make sure that you verify your results
using an independent method if at all possible (for example using the fitting
tool from http://caravan.astro.wisc.edu/protostars/ which uses an independent
Fortran code).

The documentation is available at http://sedfitter.readthedocs.org

Status
------

[![Build Status](https://travis-ci.org/astrofrog/sedfitter.png?branch=refactor-non-affiliated)](https://travis-ci.org/astrofrog/sedfitter)

[![Coverage Status](https://coveralls.io/repos/astrofrog/sedfitter/badge.png)](https://coveralls.io/r/astrofrog/sedfitter)
