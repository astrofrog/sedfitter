Model selection syntax
======================

Many of the functions in the SED fitter accept a tuple containing a single character string and a value. This tuple is used to indicate how many models to plot, output, keep, etc.

The character should indicate what kind of selection is being done, and the number quantifies the selection. The options are:

* ``('A', value)``: select all models (``value`` is ignored).

* ``('N', value)``: select a fixed number of fits given by ``value``.

* ``('C', value)``: select all models with a :math:`\chi^2` below a threshold given by ``value``.

* ``('D', value)``: select all models with a :math:`\chi^2-\chi^2_{\rm best}` value below a threshold given by ``value``.

* ``('E', value)``: select all models with a :math:`\chi^2` per datapoint below a threshold given by ``value``.

* ``('F', value)``: select all models with a :math:`\chi^2-\chi^2_{\rm best}` value per datapoint below a threshold given by ``value``.


