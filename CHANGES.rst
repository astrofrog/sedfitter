0.9.2 (unreleased)
------------------

- Added support for a new model directory format that stores all SEDs in a
  single file. [#16]

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
