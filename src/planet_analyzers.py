"""Some functions to calculate physical properties from density profiles."""

import sys
import numpy as np

def mass(svec, dvec):
    """Return total mass implied by density profile.

    Parameters
    ----------
    svec : vector
        Mean radii of const. density surfaces, in real units.
    dvec : vector
        Density on corresponding radius svec, in real units

    Returns
    -------
    M : scalar
        Total mass implied by the density profile.

    NOTE: svec, dvec are typically the output from one of the functions in module
    planet_generators.py.
    """

    dro = np.hstack((dvec[0], np.diff((dvec))))
    m = 4*np.pi/3*sum(dro*(svec)**3)
    return m

def mass_variable(svec, dvec):
    """Return cumulative mass implied by density profile.

    Parameters
    ----------
    svec : vector
        Mean radii of const. density surfaces, in real units.
    dvec : vector
        Density on corresponding radius svec, in real units

    Returns
    -------
    mvec : vector
        mvec[i] is the mass inside of mean radius svec[i]. Our density profiles
        are given with the outer radius in svec[0], so that mvec[0] is the total
        planetary mass.

    NOTE: svec, dvec are typically the output from one of the functions in module
    planet_generators.py.
    """

    mvec = np.zeros(svec.shape)
    N = mvec.shape[0]
    mvec[-1] = 4*np.pi/3*dvec[-1]*svec[-1]**3 # the last layer is a sphere
    for k in range(N-1,0,-1):
        mvec[k-1] = mvec[k] + 4*np.pi/3*dvec[k-1]*(svec[k-1]**3 - svec[k]**3)
    return mvec

def _test():
    print("alo world")
    N = 12
    svec = np.linspace(1, 1/N, N)
    dvec = np.ones(svec.shape)
    m = mass_variable(svec, dvec)
    print(m)

if __name__ == "__main__":
    _test()
