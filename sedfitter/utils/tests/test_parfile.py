import os
import pytest

from ..parfile import read

ROOT = os.path.dirname(__file__)

@pytest.mark.parametrize('format', ('par', 'conf'))
def test_parfile(format):
    p = read(os.path.join(ROOT, 'data', 'test.' + format), format)
    assert p['a'] == 'hello, world!'
    assert p['b'] == 1
    assert p['c'] == 2.3
    assert p['d'] == '1 + 2'
