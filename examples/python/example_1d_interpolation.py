"""
Example of the Python binding of the Gpuspline library for the
calculation of multidimensional cubic splines.

Interpolates 1D data. The data is upsampled, cut, stretched and shifted.

Requires pyGpuspline, Numpy and Matplotlib
"""

import numpy as np
from matplotlib import pyplot as plt
import pygpuspline.gpuspline as gs

if __name__ == '__main__':
    # input data
    y = np.array([0, 0, 0.2, 1, 1.1, 1.3, 2, 2.5, 3, 4, 4.25, 4, 3, 2.5, 2, 1.3, 1.1, 1, 0.2, 0, 0], np.float32)
    x = np.arange(y.size)
    center = x[-1] / 2

    # interpolation parameters
    edge = 1.4
    width = 1.1
    shift = 1.2
    sampling_factor = 0.5

    # interpolation
    xq = np.arange(x[0], x[-1], sampling_factor, np.float32)
    xq = xq[np.logical_and(xq >= edge, xq <= np.amax(xq) - edge)]
    xq /= width
    xq += center * (1 - 1 / width) - shift
    yq = gs.spline_interpolate(y, xq)

    # show result
    fig, ax = plt.subplots()
    ax.plot(x, y, color='blue', label='original')
    ax.plot(xq + shift, yq, color='red', marker='x', label='interpolated')
    ax.grid()
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 1.1 * np.amax(y))
    ax.legend()
    plt.show()
