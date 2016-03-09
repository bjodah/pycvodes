from math import exp


if __name__ == '__main__':
    import pyximport
    pyximport.install()
    from _cvodes_cxx import PyDecay

    pd = PyDecay(1.0)
    tout, yout = pd.adaptive(1.0, 1.0)
    for t, y in zip(tout, yout):
        assert abs(y - exp(-t)) < 1e-9
