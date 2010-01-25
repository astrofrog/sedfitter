==================
Preparing the data
==================

The data file you will need to feed to the fitter should contain one source per line. Each line should contain ``3*(n+1)`` columns, where ``n`` is the number of wavelengths/filters for which fluxes are given:

* Column ``1`` should contain the source name, with no spaces
* Columns ``2`` and ``3`` should contain the source coordinates
* Columns ``4`` to ``3 + n`` should contain a flag for each flux, where the flag can take the following values:

  * ``0`` for data points which are not valid or used
  * ``1`` for valid data points
  * ``2`` for lower limits
  * ``3`` for upper limits
  * ``4`` for valid data points (as for ``1``) but where the fluxes and errors are expressed as :math:`\log_{10}{\rm [flux]}` and :math:`\log_{10}{\rm [error]}`. This is useful for sources which have a magnitude with a large error, which can then be converted to :math:`\log_{10}{\rm [flux]}` and :math:`\log_{10}{\rm [error]}` without worrying about the conversion from magnitude and magnitude error to linear fluxes and errors.

* Columns ``4 + n`` to ``3(n + 1)`` should contain the fluxes and flux uncertainties, alternated, in mJy. For lower and upper limits, the error should be set to the confidence placed on the limit, with 0 corresponding to 0% and 1 corresponding to 100%. For example, in the case of an upper limit, setting the confidence to 1 effectively forbids any models to be larger than the upper limit, while setting the confidence between 0 and 1 allows such models, but introduces a :math:`\chi^2` penalty dependent on the confidence. Setting the confidence to 0 effectively disables the upper limit.

We refer to this format as the **fitter data format**. Note that in a single data file, all lines should have the same number of columns. If one or more fluxes are not available for a specific source, the flag for these fluxes can be set to ``0``.
