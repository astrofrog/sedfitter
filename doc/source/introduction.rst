============
Introduction
============

The SED fitter was originally written to help analyze GLIMPSE data, and can now be used with any set of photometric data. The concept is very simple: consider a set of sources to be studied. For each source, the SED fitter can fit models, including model stellar photospheres, YSO model SEDs, as well as galaxy and AGB templates, to the multi-wavelength photometry measurements of this particular source using linear regression.

The scale factor :math:`S` (which is related to the luminosity and distance of the sources) and the extinction :math:`A_{\rm V}` are used as free parameters in the fitting process. The result for any given source is a value for the goodness of fit and best-fit values :math:`S` and :math:`A_{\rm V}` for every single model. These fits can then be analyzed to derive properties of the source.

To quantify the goodness/badness of each fit to each source, we calculate the :math:`\chi^2` value:

.. math::
    \chi^2 = \sum_{i=1}^{n}\left(\frac{F_{i}-M_{i}}{\sigma_{i}^2}\right)^2

where :math:`F_{i}` are the flux values at a given wavelength :math:`\lambda_{i}`, :math:`\sigma_{i}` are the flux errors, and :math:`M_{i}` are the extincted and scaled model values. In this manual we sometimes refer to the :math:`\chi^2` per datapoint, :math:`n_{\rm data}`, where the number of datapoints does *not* include upper and lower limits.

There are two kinds of models that can be used with the SED fitter:

* Models for which there is no aperture dependence on the flux and for which the absolute distance cannot be determined (e.g. stellar photosphere models)
* Models for which the flux depends on the aperture chosen and/or for which the distance can be found from the scalefactor of the fit. (e.g. YSO models). In this case, the fitting procedure is more complex as it involves computing the aperture-dependent SEDs for a fine grid of distances, and optionally removing models that would clearly be extended relative to the aperture chosen. This is described in more detail in Paper II.

Two binaries are available for fitting, optimized for each case: ``fit_stellar`` should be used for models with no aperture dependence and no distance information, while ``fit`` should be used for models with aperture-dependent fluxes and/or absolutely scaled models. Unlike the previous versions of the fitter, only one set of models can be used at a time, however, it is easy to produce a new set of models without changing a single line of code in the fitter distribution. If you wish to produce your own package of models, see the ``model_packages.pdf`` file.

