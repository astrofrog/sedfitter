## v1.4 - 2023-02-06

### What's Changed

- Small astropy 4.0 compatibility refactor by @keflavich in https://github.com/astrofrog/sedfitter/pull/67
- Update package infrastructure and fix compatibility with latest astropy versions by @astrofrog in https://github.com/astrofrog/sedfitter/pull/79
- ndistances can't be a float in recent versions of numpy by @keflavich in https://github.com/astrofrog/sedfitter/pull/73
- Use memmap when reading in model grids by @keflavich in https://github.com/astrofrog/sedfitter/pull/78
- More infrastructure improvements by @astrofrog in https://github.com/astrofrog/sedfitter/pull/80

**Full Changelog**: https://github.com/astrofrog/sedfitter/compare/v1.3...v1.4

## v1.3 - 2019-04-25

- Fixed an issue that occurred during the convolution if not all SEDs
- had the same number of wavelengths. [#66]

## v1.2 - 2018-06-16

- Use install_requires instead of requires.

## v1.1 - 2018-06-16

- Fix color of lines in SED plots. [#52]
- Fix Windows compatibility. [#53]
- Fix compatibility with Astropy 2.x and above. [#60]

## v1.0 - 2016-01-09

- Fixed a bug in the calculation of the reduced chi^2: this should
- have used the number of valid data points (valid flag 1 or 4) but
- instead was using the total number of valid and invalid points.
- [#49]
- Fixed compatibility with Python 3.6. [#49]

## v0.9.6 - 2016-07-31

- A valid flag of '9' can now be used to denote a flux that should
- be plotted but not fit. [#32]
- Fix issues that occurred if no valid fits were present for a given
- source. [#38]
- Switch the order of the dimensions in SED flux cubes to optimize
- performance. This means that this version of the SED fitter will be
- incompatible with flux cubes generated with earlier versions, but
- will result in significant performance improvements when plotting
- SEDs. [#44]
- Add a `memmap` option to `convolve_model_dir` that controls whether
- memory mapping is used to read the flux cubes (in the case where
- SEDs are stored in cubes rather than individual files). If the cubes
- can fit in memory, the convolution is much faster if the memory
- mapping is explicitly turned off.

## v0.9.5 - 2015-06-07

- Fixed calculation of indices for monochromatic 'convolution'.

## v0.9.4 - 2015-06-05

- Added missing ez_setup.py file.

## v0.9.3 - 2015-06-05

- Fixed a bug with filter normalization if filter was given in
- increasing wavelength. [#26]
- Fixed plotting of source names when using Latex. [#25]

## v0.9.2 - 2014-12-01

- Added support for a new model directory format that stores all SEDs
- in a single file. [#16]
- Added a new Fitter class that provides an OO interface to the
- fitter. This makes it easier to fit sources one by one without
- reloading models. [#12]
- Various bug fixes.

## v0.9.1 - 2014-03-10

- Added documentation on convolution.
- Added page about accessing model packages.
- Renamed `wavelength` attribute on `Filter` and `ConvolvedFluxes` to
- `central_wavelength`, renamed `r` attribute on `Filter` to
- `response`, and removed `wav` attribute for `Filter` since it was
- reduntant with `nu`.

## v0.9.0 - 2014-01-30

- Initial release
