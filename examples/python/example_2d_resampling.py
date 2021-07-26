"""

"""

import numpy as np
from matplotlib import pyplot as plt
import pygpuspline.gpuspline as gs


def calculate_psf(x, y, p):
    """

    """
    sx = p[3] - 0.2
    sy = p[3] + 0.2

    psf = p[0] * np.exp(-0.5*(((x-p[1])/sx)**2+((y-p[2])/sy)**2))+p[4]

    return psf


if __name__ == '__main__':

    # PSF size
    size_x = 10
    size_y = 20

    # derived values
    x = np.arange(size_x, dtype=np.float32).reshape((size_x, 1))
    y = np.arange(size_y, dtype=np.float32).reshape((1, size_y))

    x_up = np.arange(size_x, step=0.1, dtype=np.float32)
    y_up = np.arange(size_y, step=0.1, dtype=np.float32)

    x_down = np.arange(size_x, step=2, dtype=np.float32)
    y_down = np.arange(size_y, step=2, dtype=np.float32)

    x_shift = x - 1.3
    y_shift = y + 2.7

    # PSF parameters
    psf_parameters = (100, (size_x-1)/2, (size_y-1)/2, 3, 10)

    # calculate PSF
    psf = calculate_psf(x, y, psf_parameters)

    # calculate spline coefficients
    coefficients = gs.spline_coefficients(psf)

    # generate upsampled PSF
    psf_up = gs.spline_values(coefficients, x_up, y_up)

    # generate downsampled PSF
    psf_down = gs.spline_values(coefficients, x_down, y_down)

    # generate shifted PSF
    psf_shift = gs.spline_values(coefficients, x_shift, y_shift)

    # display results
    fig, axs = plt.subplots(2,2)
    fig.tight_layout()
    axs = axs.flat
    axs[0].imshow(psf)
    axs[0].set_title('PSF')
    axs[1].imshow(psf_up)
    axs[1].set_title('Upsampled PSF')
    axs[2].imshow(psf_down)
    axs[2].set_title('Downsampled PSF')
    axs[3].imshow(psf_shift)
    axs[3].set_title('Shifted PSF')
    plt.show()
