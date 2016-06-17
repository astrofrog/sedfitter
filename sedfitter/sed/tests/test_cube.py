import numpy as np
from astropy import units as u
from ..cube import SEDCube


def test_roundrip(tmpdir):

    n_models = 30
    n_ap = 3
    n_wav = 10

    s = SEDCube()

    s.names = ['name_{0:02d}'.format(i) for i in range(n_models)]

    s.apertures = np.linspace(10, 100, n_ap) * u.au

    s.wav = np.linspace(0.01, 5000, n_wav)[::-1] * u.micron

    s.distance = 1. * u.kpc

    s.val = np.random.random((n_models, n_ap, n_wav)) * u.mJy
    s.unc = np.random.random((n_models, n_ap, n_wav)) * u.mJy

    temp_file = tmpdir.join('test_roundtrip_sedcube').strpath

    s.write(temp_file)

    s2 = SEDCube.read(temp_file)

    assert s == s2


def test_roundrip_missing_optional(tmpdir):

    n_models = 30
    n_wav = 10

    s = SEDCube()

    s.names = ['name_{0:02d}'.format(i) for i in range(n_models)]

    s.wav = np.linspace(0.01, 5000, n_wav)[::-1] * u.micron

    s.distance = 1. * u.kpc

    s.val = np.random.random((n_models, 1, n_wav)) * u.mJy

    temp_file = tmpdir.join('test_roundtrip_sedcube').strpath

    s.write(temp_file)

    s2 = SEDCube.read(temp_file)

    assert s == s2
