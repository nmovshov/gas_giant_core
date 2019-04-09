"""A small collection of plotting methods specific to planet interior models."""

import sys
import numpy as np
import matplotlib.pyplot as plt

def rho_of_s(svec, dvec, *args, **kwargs):
    """Plot density vs. normalized mean radius."""

    # Some minimal input control
    assert(svec.shape == dvec.shape)

    # Prepare the canvas
    plt.figure()

    # Prepare the data
    x = svec/svec[0] # normalized radius
    y = dvec/1000    # density in 1000 kg/m^3

    # Plot a line
    plt.plot(x, y, *args, **kwargs)

    # Style and annotate
    plt.xlabel('Level surface mean radius [normalized]')
    plt.ylabel('density [1000 kg/m$^3$]')

    # Display and return
    plt.show()

if __name__ == '__main__':
    print("alo world")
