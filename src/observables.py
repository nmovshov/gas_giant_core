"""Observed planetary values in consistent format.

The purpose of the observables module is to create a consistent and predictable
format for a structure containing observed values of a planet's vitals. The
typical usage is:

    from observables import <planet_name>

The returned struct has the following fields:

  obs.M, obs.dM        -  planet's mass in kg, with uncertainty
  obs.a0               -  planet's equatorial radius in km.
  obs.s0               -  planet's mean radius in km.
  obs.P0               -  reference surface pressure on obs.a0
  obs.q                -  planet's dimensionless rotation parameter, q=w^2*a0^2/(GM)
  obs.m                -  planet's dimensionless rotation parameter, m=w^2*s0^2/(GM)
  obs.J<n>, obs.dJ<n>  -  planet's n-th gravity harmonic with uncertainty
  obs.name
"""

class Jupiter:
    M  = 1.8986112e27   # Guillot et al. (1994) Table I
    dM = 4.7e-5*M       # https://physics.nist.gov/cuu/index.html
    a0 = 7.1492e7       # Seidelmann et al. (2007) table 4
    s0 = 6.9911e7       # Seidelmann et al. (2007) table 4
    P0 = 1e5            # The reference radius is the 1 bar level
    q  = 0.088822426    # Hubbard (2013) Table 1
    m  = q*0.935292     # My models' typical (s/a)^3 value
    J2  = 14696.572e-6  # Iess et al. 2018 Table 1
    J4  =  -586.609e-6  # Iess et al. 2018 Table1
    J6  =    34.198e-6  # Iess et al. 2018 Table1
    J8  =    -2.426e-6  # Iess et al. 2018 Table1
    J10 =     0.172e-6  # Iess et al. 2018 Table1
