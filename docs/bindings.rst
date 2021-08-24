.. _external-bindings:

=================
External bindings
=================

This sections describes the Gpuspline bindings to other programming languages. The bindings to Python and Matlab aim to
emulate the :ref:`c-interface` as closely as possible.

Python
------

The Python binding for Gpuspline consists of a Python package pyGpuspline that provides various functions
that call the C interface of the Gpuspline library. In general the routines expect data as NumPy arrays.

Installation
++++++++++++

Wheel files for Python (x64) on Windows are included in the binary package. NumPy is required.

Install the wheel file with

.. code-block:: bash

    pip install --no-index --find-links=LocalPathToWheelFile pyGpuspline

Python Interface
++++++++++++++++

The Python interface is a thin wrapper around the C interface. Please see the API documentation of the
:ref:`c-interface` for more details on the interpretation of input and output parameters.

spline_coefficients
...................

The signature of the spline_coefficients method is

.. code-block:: python

    def spline_coefficients(data)

The data must be a 1-3D NumPy array of data type single. This method is equivalent to call the C interface
functions calculate_coefficients_Xd. The return value is a NumPy array of size 4^d (d=dimension of data) times product of
number of pixels in each direction minus one.

spline_interpolate
..................

The signature of the spline_interpolate method is

.. code-block:: python

    def spline_interpolate(data, x, y=None, z=None):

The data must be a 1-3D NumPy array of data type single. This method is equivalent to call the C interface
functions interpolate_Xd. In case of a 1D data array only x values should be specified, in case of a 2D array
x and y values and for a 3D data array x, y and z values.

The output is the interpolated data, a NumPy array with product of elements in x times elements in y  times
elements in z entries.

spline_values
.............

The signature of the spline_values method is

.. code-block:: python

    def spline_values(coefficients, x, y=None, z=None):

The coefficients are a (d+1) dimensional NumPy array of type single (d=dimension of original data) describing the spline
associated with the data. In the first dimension with 4^d entries are the coefficients of a single spline interval.
The spline intervals for each data pixel are stored from the second dimension on. For a 1D spline, only specify x, for a
2D spline only x and y and for a 3D spline only x, y and z. The output is a NumPy array holding the interpolated
values of the data represented by the splines at the positions specified by x, y and z. This method is equivalent
to call functions calculate_values_Xd in the C interface.

Python Examples
+++++++++++++++

1D interpolation example
........................

An example for interpolating data points calling a cubic spline interpolation routine implemented in C.1D data is
upsampled, cut, stretched and shifted. The example can be found at `example_1d_interpolation.py`_.

.. code-block:: python

    """
    Example of the Python binding of the Gpuspline library for the
    calculation of multidimensional cubic splines.
    https://github.com/gpufit/Gpuspline
    https://gpuspline.readthedocs.io/en/latest/bindings.html#python

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
        yq = gs.spline_interpolate(y, xq)  # call to the spline library

        # show result
        fig, ax = plt.subplots()
        ax.plot(x, y, color='blue', label='original')
        ax.plot(xq + shift, yq, color='red', marker='x', label='interpolated')
        ax.grid()
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 1.1 * np.amax(y))
        ax.legend()
        plt.show()

2D resampling and shifting example
..................................

The example can be found at `example_2d_resampling.py`_.

.. code-block:: python

    """
    Example of the Matlab binding of the Gpuspline library for the
    calculation of multidimensional cubic splines.
    https://github.com/gpufit/Gpuspline
    https://gpuspline.readthedocs.io/en/latest/bindings.html#python

    2D data is interpolated (up- and downsampled and shifted).

    Requires pyGpuspline, Numpy and Matplotlib
    """

    import numpy as np
    from matplotlib import pyplot as plt
    import pygpuspline.gpuspline as gs


    def calculate_psf(x, y, p):
        """
        Calculates an elliptic 2D Gaussian peak function.
        """
        sx = p[3] - 0.2
        sy = p[3] + 0.2

        psf = p[0] * np.exp(-0.5 * (((x - p[1]) / sx) ** 2 + ((y - p[2]) / sy) ** 2)) + p[4]

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
        psf_parameters = (100, (size_x - 1) / 2, (size_y - 1) / 2, 3, 10)

        # calculate PSF
        psf = calculate_psf(x, y, psf_parameters)

        # calculate spline coefficients
        coefficients = gs.spline_coefficients(psf)  # call to spline library

        # generate upsampled PSF
        psf_up = gs.spline_values(coefficients, x_up, y_up)  # call to spline library

        # generate downsampled PSF
        psf_down = gs.spline_values(coefficients, x_down, y_down)  # call to spline library

        # generate shifted PSF
        psf_shift = gs.spline_values(coefficients, x_shift, y_shift)  # call to spline library

        # display results
        fig, axs = plt.subplots(2, 2)
        fig.tight_layout()
        axs = axs.flat
        axs[0].imshow(psf, cmap='hot')
        axs[0].set_title('Original data')
        axs[1].imshow(psf_up, cmap='hot')
        axs[1].set_title('Upsampled')
        axs[2].imshow(psf_down, cmap='hot')
        axs[2].set_title('Downsampled')
        axs[3].imshow(psf_shift, cmap='hot')
        axs[3].set_title('Shifted')
        plt.show()


Matlab
------

The Matlab binding for Gpuspline consists of Matlab scripts (spline_coefficients.m, spline_values.m,
spline_interpolate.m). These scripts check the input data and call the :ref:`c-interface` of the Gpuspline library, via
compiled .mex files. Please note, that before using the Matlab binding, the path to the .m and .mex files must be added
to the Matlab path.

Matlab Interfaces
+++++++++++++++++

spline_coefficients
...................

The data dimensions are deduced from the dimensions of the input data.

The signature of the spline_coefficients function is

.. code-block:: matlab

    function [coefficients, time] = spline_coefficients(data)

*Input parameters*

:data: Data |br|
    1D, 2D or 3D matrix of data type single.

*Output parameters*

:coefficients: Cubic spline coefficients |br|
    2D, 3D or 4D matrix of data type single of size: |br|
    [N\ :sub:`c`, N\ :sub:`x`, N\ :sub:`y`, N\ :sub:`z`] |br|
    ,where N\ :sub:`c` represents the number of coefficients per spline interval and
    N\ :sub:`x`, N\ :sub:`y`, N\ :sub:`z`
    the number of spline intervals in x, y, z.
:time: Execution time of call to spline_coefficientsMex in seconds.

Errors are raised if checks on parameters fail or if the execution of the function fails.

spline_values
.............

The spline dimensions are deduced from the dimensions of the input data and the number of input arguments. This function
calculates function values based on cubic spline coefficients. Optionally it calculates values for multiple splines, if
the 5th dimension of the input spline coefficients is greater than 1.

The signature of the spline_values function is

.. code-block:: matlab

    function [coefficients, time] = spline_values(coefficients, x, y, z)

*Input parameters*

:coeficients: Cubic spline coefficients |br|
    2D, 3D, 4D or 5D matrix of data type single.

    :dimension1: Number of spline coefficients per spline interval depending on the number of dimensions of the spline (4, 16 or 64)

    :dimension2: Number of spline intervals in x

    :dimension3: Number of spline intervals in y

    :dimension4: Number of spline intervals in z

    :dimension5: Number of splines/channels

:x: Independent variable x values

:y: Independent variable y values (optional)

:z: Independent variable z values (optional)

y and z parameter are optional (for 2D/3D data)

*Output parameters*

:values: Output values |br|
    1D, 2D, 3D or 4D matrix of data type single of size: |br|
    [N\ :sub:`x`, N\ :sub:`y`, N\ :sub:`z`, N\ :sub:`ch`] |br|
    ,where N\ :sub:`x`, N\ :sub:`y`, N\ :sub:`z` represent the number of output data points in x, y, z and
    N\ :sub:`ch` the number of channels.
:time: Execution time of call to spline_valuesMex in seconds.

Errors are raised if checks on parameters fail or if the execution of the function fails.

spline_interpolate
..................

The data dimensions are deduced from the number of input arguments.

The signature of the spline_interpolate function is

.. code-block:: matlab

    function [interpolated_data, time] = spline_interpolate(data, x, y, z)

*Input parameters*

:data: Input data values |br|
    1D, 2D or 3D matrix of data type single.

:x: Independent variable x values |br|
    1D matrix of data type single.

:y: Independent variable y values (optional) |br|
    1D matrix of data type single.

:z: Independent variable z values (optional) |br|
    1D matrix of data type single.

y and z parameter are optional (for 2D/3D interpolation)

*Output parameters*

:values: Interpolated data values |br|
    1D, 2D or 3D matrix of data type single of size: |br|
    [N\ :sub:`x`, N\ :sub:`y`, N\ :sub:`z`] |br|
    ,where N\ :sub:`x`, N\ :sub:`y`, N\ :sub:`z` represent the number of output data points in x, y, z.
:time: Execution time of call to spline_interpolateMex in seconds.

Errors are raised if checks on parameters fail or if the execution of the function fails.


Matlab Examples
+++++++++++++++

1D interpolation example
........................

An example for interpolating data points calling a cubic spline interpolation routine implemented in C.1D data is
upsampled, cut, stretched and shifted. The example can be found at `example_1d_interpolation.m`_.

2D resampling example
.....................

Example can be found at `example_2d_resampling.m`_.

.. code-block:: matlab

    function example_2d_resampling()
    % Example of the Matlab binding of the Gpuspline library for the
    % calculation of multidimensional cubic splines.
    % https://github.com/gpufit/Gpuspline
    %
    % 2D data is interpolated (up- and downsampled).
    % https://gpuspline.readthedocs.io/en/latest/bindings.html#matlab

    %% psf size
    size_x = 15;
    size_y = 20;

    %% derived values
    x = single(0 : size_x - 1)';
    y = single(0 : size_y - 1);

    x_up = single(0 : 0.1 : size_x - 1)';
    y_up = single(0 : 0.1 : size_y - 1)';

    x_down = single(0 : 2 : size_x - 1)';
    y_down = single(0 : 2 : size_y - 1)';

    %% PSF parameters
    psf_parameters = single([100, (size_x-1)/2, (size_y-1)/2, 3, 10]);

    %% calculate PSF
    psf = calculate_psf(x, y, psf_parameters);

    %% calculate spline coefficients
    coefficients = spline_coefficients(psf);

    %% generate upsampled psf
    psf_up = spline_values(coefficients, x_up, y_up);

    %% generate downsampled psf
    psf_down = spline_values(coefficients, x_down, y_down);

    %% figure
    figure;
    subplot(131); imagesc(x, y, psf);
    axis image; title('Original data');
    subplot(132); imagesc(x_up, y_up, psf_up);
    axis image; title('Upsampled');
    subplot(133); imagesc(x_down, y_down, psf_down);
    axis image; title('Downsampled');
    colormap('hot');

    end

    function psf = calculate_psf(x, y, p)
    % PSF consists of an elliptic 2D Gaussian

    % p(1) - amplitude
    % p(2) - center x
    % p(3) - center y
    % p(4) - Standard deviation
    % p(5) - constant background
    assert(nargin == 3);

    sx = p(4) - 0.2;
    sy = p(4) + 0.2;

    arg_ex = exp(-1/2*((x-p(2))/sx).^2-1/2*((y-p(3))/sy).^2);

    psf = p(1) .* arg_ex + p(5); % scale with amplitude and background

    end


example_2d_shift()
..................

Example can be found at `example_2d_shift.m`_.

.. code-block:: matlab

    function example_2d_shift()
    % Example of the Matlab binding of the Gpuspline library for the
    % calculation of multidimensional cubic splines.
    % https://github.com/gpufit/Gpuspline
    %
    % 2D data is interpolated (shifted).
    % https://gpuspline.readthedocs.io/en/latest/bindings.html#matlab

    %% psf size
    size_x = 20;
    size_y = 30;

    %% derived values
    x = single(0 : size_x - 1)';
    y = single(0 : size_y - 1);

    x_shifted = x - 1.3;
    y_shifted = y + 2.7;

    %% PSF parameters
    psf_parameters = single([100, (size_x-1)/2, (size_y-1)/2, 2, 10]);

    %% calculate PSF
    psf = calculate_psf(x, y, psf_parameters);

    %% calculate spline coefficients
    coefficients = spline_coefficients(psf);

    %% generate upsampled psf
    psf_shifted = spline_values(coefficients, x_shifted, y_shifted);

    %% figure
    figure;
    subplot(121); imagesc(x, y, psf);
    axis image; title('Original');
    subplot(122); imagesc(x_shifted, y_shifted, psf_shifted);
    axis image; title('Shifted');
    colormap('hot');

    end

    function psf = calculate_psf(x, y, p)
    % PSF consists of an elliptic 2D Gaussian

    % p(1) - amplitude
    % p(2) - center x
    % p(3) - center y
    % p(4) - Standard deviation
    % p(5) - constant background
    assert(nargin == 3);

    sx = p(4) - 0.2;
    sy = p(4) + 0.2;

    arg_ex = exp(-1/2*((x-p(2))/sx).^2-1/2*((y-p(3))/sy).^2);

    psf = p(1) .* arg_ex + p(5); % scale with amplitude and background

    end