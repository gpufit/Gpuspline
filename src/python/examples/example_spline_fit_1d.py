"""

"""

import numpy as np
from matplotlib import pyplot as plt
import pygpuspline.gpuspline as gs
import pygpufit.gpufit as gf

if __name__ == '__main__':
    # data size
    size_x = 25

    # fit parameter
    tolerance = 1e-30
    max_n_iterations = 100
