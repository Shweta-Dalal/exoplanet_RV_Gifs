import numpy as np
from exoplanet_plottingRV import radvelocity

def test_radvelocity_zero_eccentricity():
    t = np.linspace(0, 10, 5)
    # For e=0 and omega=0, vr = K * cos(2π t/P)
    P, K = 5.0, 2.0
    vr = radvelocity(t, P, K, 0.0, 0.0, 0.0)
    expected = K * np.cos(2*np.pi*t/P)
    assert np.allclose(vr, expected)
