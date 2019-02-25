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

    m = 0 # ALICE, YOUR CODE GOES HERE
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

    mvec = 0 # EMMA, YOUR CODE GOES HERE
    return mvec

def _test():
    print("alo world")

if __name__ == "__main__":
    _test()
