"""
Helpers for model parameter tables that include vector-valued columns.

Some model grids (for example the Richardson et al. 2024 YSO models) include
parameters that are defined for each aperture, stored as columns whose last
dimension matches the apertures of the flux cube. These helpers turn such
tables into one scalar value per model, either by interpolating the
aperture-dependent columns to a given aperture or by splitting them into one
column per aperture.
"""

import numpy as np
from astropy import units as u

__all__ = ['aperture_dependent_columns', 'aperture_to_au', 'interpolate_to_aperture', 'parameter_names',
           'parameter_values', 'select_parameter', 'validate_aperture']

# Placeholder index meaning that a column should be interpolated to the
# requested aperture rather than indexed.
INTERPOLATE = 'interpolate'


def _is_aperture_dependent(shape, apertures):
    return apertures is not None and len(shape) == 1 and shape[0] == len(apertures)


def _layout(table, apertures, interpolate):
    """
    Yield ``(output_name, column_name, index)`` for each scalar output column.

    ``index`` is a tuple of indices into the trailing dimensions of the
    column, or `INTERPOLATE` for aperture-dependent columns that should be
    interpolated.
    """
    for column_name in table.colnames:
        if column_name == 'MODEL_NAME':
            continue
        shape = table[column_name].shape[1:]
        if shape == ():
            yield column_name, column_name, ()
        elif _is_aperture_dependent(shape, apertures):
            if interpolate:
                yield column_name, column_name, INTERPOLATE
            else:
                for i, ap in enumerate(apertures.to_value(u.au)):
                    yield '{0}_{1:.0f}au'.format(column_name, ap), column_name, (i,)
        else:
            for index in np.ndindex(shape):
                yield '{0}_{1}'.format(column_name, '_'.join(map(str, index))), column_name, index


def aperture_dependent_columns(table, apertures):
    """
    Return the names of the columns that are defined for each aperture.
    """
    return [name for name in table.colnames
            if _is_aperture_dependent(table[name].shape[1:], apertures)]


def parameter_names(table, apertures, interpolate=False):
    """
    Return the names of the scalar columns that `parameter_values` returns.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table of model parameters
    apertures : `~astropy.units.Quantity` or `None`
        The apertures the model fluxes are defined at, used to identify
        aperture-dependent columns
    interpolate : bool, optional
        Whether aperture-dependent columns will be interpolated to a single
        aperture (one output column each) rather than split into one column
        per aperture
    """
    return [name for name, _, _ in _layout(table, apertures, interpolate)]


def parameter_values(table, apertures, aperture_au=None):
    """
    Return one array of scalar values per model for each output column.

    Scalar columns are returned as-is. Aperture-dependent columns are
    interpolated to ``aperture_au`` if specified, and are otherwise split into
    one column per aperture. Any other vector columns are split into one column
    per element. The order matches `parameter_names`.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table of model parameters
    apertures : `~astropy.units.Quantity` or `None`
        The apertures the model fluxes are defined at
    aperture_au : `~astropy.units.Quantity`, optional
        The aperture to interpolate aperture-dependent columns to, either as a
        scalar or as one value per row of ``table``
    """
    values = []
    for _, column_name, index in _layout(table, apertures, aperture_au is not None):
        column = table[column_name]
        if index is INTERPOLATE:
            values.append(interpolate_to_aperture(column, apertures, aperture_au))
        else:
            values.append(column[(slice(None),) + index])
    return values


def select_parameter(table, name, apertures, aperture_au=None):
    """
    Return the values of a single parameter, with one scalar per model.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table of model parameters
    name : str
        The parameter to return
    apertures : `~astropy.units.Quantity` or `None`
        The apertures the model fluxes are defined at
    aperture_au : `~astropy.units.Quantity`, optional
        The aperture to interpolate the parameter to if it is aperture-dependent
    """
    column = table[name]
    shape = column.shape[1:]
    if shape == ():
        return column
    elif _is_aperture_dependent(shape, apertures):
        if aperture_au is None:
            raise ValueError("Parameter '{0}' depends on aperture, so an aperture "
                             "needs to be specified".format(name))
        return interpolate_to_aperture(column, apertures, aperture_au)
    else:
        raise ValueError("Parameter '{0}' has shape {1} per model and cannot "
                         "be reduced to a single value".format(name, shape))


def aperture_to_au(aperture, log_distance_kpc):
    """
    Convert an aperture radius to au for each fit.

    Parameters
    ----------
    aperture : `~astropy.units.Quantity`
        A scalar aperture radius, either as a length or as an angle. Angles
        are converted to lengths using the distance of each fit.
    log_distance_kpc : `~numpy.ndarray`
        The log10 of the distance of each fit in kpc (the ``sc`` attribute of
        :class:`~sedfitter.fit_info.FitInfo`).
    """
    log_distance_kpc = np.asarray(log_distance_kpc)
    if aperture.unit.physical_type == 'angle':
        return aperture.to_value(u.arcsec) * 10. ** log_distance_kpc * 1000. * u.au
    else:
        return np.broadcast_to(aperture.to(u.au), log_distance_kpc.shape, subok=True)


def validate_aperture(aperture, allow_angle=True):
    """
    Check that ``aperture`` is a scalar length (or angle if allowed).
    """
    if aperture is None:
        return None
    allowed = ('length', 'angle') if allow_angle else ('length',)
    if (not isinstance(aperture, u.Quantity) or not aperture.isscalar or
            aperture.unit.physical_type not in allowed):
        raise TypeError("aperture should be a scalar Quantity with units of "
                        "{0}".format(' or '.join(allowed)))
    return aperture


def interpolate_to_aperture(values, apertures, aperture_au):
    """
    Interpolate per-aperture values to a given aperture.

    The interpolation is linear in log10(aperture). Apertures larger than the
    largest model aperture use the values at the largest aperture, matching
    how model fluxes are treated, while apertures smaller than the smallest
    model aperture raise an error.

    Parameters
    ----------
    values : array-like
        The values, with shape ``(n_models, n_apertures)``
    apertures : `~astropy.units.Quantity`
        The model apertures, in increasing order
    aperture_au : `~astropy.units.Quantity`
        The aperture to interpolate to, either as a scalar or with one value
        per model
    """
    values = np.ma.filled(np.ma.asarray(values, dtype=float), np.nan)
    log_ap = np.log10(apertures.to_value(u.au))
    target = np.log10(np.broadcast_to(aperture_au.to_value(u.au), values.shape[:1]))
    if np.any(target < log_ap[0] - 1e-10):
        raise ValueError("Aperture(s) requested are smaller than the smallest "
                         "model aperture ({0:g})".format(apertures.min()))
    target = np.clip(target, log_ap[0], log_ap[-1])
    upper = np.clip(np.searchsorted(log_ap, target), 1, len(log_ap) - 1)
    weight = (target - log_ap[upper - 1]) / (log_ap[upper] - log_ap[upper - 1])
    rows = np.arange(len(values))
    lower_values, upper_values = values[rows, upper - 1], values[rows, upper]
    # Avoid propagating NaN values from a neighbouring aperture that has zero weight
    return np.where(weight == 0, lower_values,
                    np.where(weight == 1, upper_values,
                             lower_values * (1. - weight) + upper_values * weight))
