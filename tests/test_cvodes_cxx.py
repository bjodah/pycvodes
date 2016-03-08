import pyximport
pyximport.install()

from _cvodes_cxx import PyDecay

from math import exp

def test_PyDecay():
    pd = PyDecay(1.0)
    tout, yout = pd.adaptive(1.0, 1.0)
    for t, y in zip(tout, yout):
        assert abs(y - exp(-t)) < 1e-9
