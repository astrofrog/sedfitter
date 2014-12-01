0.9.2 (2014-12-01)
------------------

- Added support for a new model directory format that stores all SEDs in a
  single file. [#16]

- Added a new Fitter class that provides an OO interface to the fitter. This
  makes it easier to fit sources one by one without reloading models. [#12]

- Various bug fixes.

0.9.1 (2014-03-10)
------------------

- Added documentation on convolution.

- Added page about accessing model packages.

- Renamed ``wavelength`` attribute on ``Filter`` and ``ConvolvedFluxes`` to
  ``central_wavelength``, renamed ``r`` attribute on ``Filter`` to
  ``response``, and removed ``wav`` attribute for ``Filter`` since it was
  reduntant with ``nu``.

0.9.0 (2014-01-30)
------------------

- Initial release
