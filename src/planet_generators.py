"""A small collection of density profile generators with consistent syntax."""

from __future__ import division
import sys
import numpy as np

def piecewise_quadratic_planet(N, x):
    """Piecewise quadratic parameterization of planet's density profile.

    piecewise_quadratic_planet(N, x) returns an N-point density profile, i.e. a
    (zvec,dvec) tuple, where the density D is approximated by a
    piecewise-quadratic function in the normalized mean radius Z=r/R_m. Each
    quadratic segment is defined by its end points and a curvature coefficient.
    The segment breakpoints are at z1 (upper) and z2 (lower). The density at z=1
    is d10 (and is usually zero). The density at z1 is d11 (the right-limit) and
    d21 the left-limit). The density at z2 is d22 (the right-limit) and d32 (the
    left-limit). The density at z=0 is d33. The curvatures of the segments are a1,
    a2, and a3. The input parameters are ordered as follows:

        x = [a1, y10, y11, a2, y21, y22, a3, y32, y33, z1, z2]
    where
        a1: curvature of first segment (upper envelope)
        y10: d10 - 0.0 (usually zero, or the 1-bar density)
        y11: d11 - d10
        a2: curvature of second segment (lower envelope)
        y21: d21 - d11
        y22: d22 - d21
        a3: curvature of third segment (core)
        y32: d32 - d22
        y33: d33 - d32
        z1: normalized mean radius of first (upper) break point
        z2: normalized mean radius of second (lower) break point

    NOTE: The densities must be specified in real units and the resulting model
    should have approximately the correct mass when the normalized radii are
    multiplied by planet's outer mean radius.
    """

    # Some minimal input control
    assert np.isscalar(N) and N > 0, "Input 1 should be positive scalar (N)"
    x = np.array(x)
    assert x.ndim == 1 and len(x) == 11, "Input 2 should be 11-vector (x)"

    # Interpreting the parameters
    a1 = x[0]; y10 = x[1]; y11 = x[2]
    a2 = x[3]; y21 = x[4]; y22 = x[5]
    a3 = x[6]; y32 = x[7]; y33 = x[8]
    z1 = x[9]; z2 = x[10]

    d10 = y10
    d11 = d10 + y11
    d21 = d11 + y21
    d22 = d21 + y22
    d32 = d22 + y32
    d33 = d32 + y33

    # Upper envelope region
    b1 = (d11 - d10)/(z1 - 1.0) - a1*(z1 + 1.0)
    c1 = d10 - a1 - b1

    # Lower envelope region
    b2 = (d22 - d21)/(z2 - z1) - a2*(z2 + z1)
    c2 = d21 - a2*z1**2 - b2*z1

    # Core region
    b3 = (d32 - d33)/z2 - a3*z2
    c3 = d33

    # Write the profile
    zvec = np.linspace(1, 1/N, N)
    dvec = np.zeros(N)
    for k in range(N):
        if zvec[k] > z1:
            dvec[k] = a1*zvec[k]**2 + b1*zvec[k] + c1
        elif zvec[k] > z2:
            dvec[k] = a2*zvec[k]**2 + b2*zvec[k] + c2
        else:
            dvec[k] = a3*zvec[k]**2 + b3*zvec[k] + c3

    # Return
    return (zvec, dvec)

def mass(svec, dvec):
    dro = np.hstack((dvec[0], np.diff((dvec))))
    m = sum(dro*(svec)**3)
    return 4*np.pi/3*m

if __name__ == '__main__':
    print("alo world")
