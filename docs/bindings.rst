.. _external-bindings:

=================
External bindings
=================

This sections describes the spline_analysis bindings to other programming languages. The bindings (to Python, Matlab or Java) aim to
emulate the :ref:`c-interface` as closely as possible.
	
Python
------

 ...

Installation
++++++++++++

 ...

Python Interface
++++++++++++++++

 ...

Python Examples
+++++++++++++++

 ...

Matlab
------

The Matlab binding for spline_analysis consists of Matlab scripts (spline_coefficients.m, spline_values.m,
spline_interpolate.m). This scripts check the input data and call the C interfaces of the spline_analysis library, via
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

Errors are raised if checks on parameters fail or if the execution of gpufit fails.

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

Errors are raised if checks on parameters fail or if the execution of gpufit fails.

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

Errors are raised if checks on parameters fail or if the execution of gpufit fails.


Matlab Examples
+++++++++++++++

1D interpolation example
........................

An example for interpolating data points calling a cubic spline interpolation routine implemented in C.1D data is
upsampled, cut, stretched and shifted. The example can be found at `example_1d_interpolation()`_.

2D resampling example
.....................

Example can be found at `example_2d_resampling()`_.

.. code-block:: matlab

    function example_2d_resampling()

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
    axis image; title('PSF');
    subplot(132); imagesc(x_up, y_up, psf_up);
    axis image; title('Upsampled PSF');
    subplot(133); imagesc(x_down, y_down, psf_down);
    axis image; title('Downsampled PSF');
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

Example can be found at `example_2d_shift()`_.

.. code-block:: matlab

    function example_2d_shift()

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
    axis image; title('PSF');
    subplot(122); imagesc(x_shifted, y_shifted, psf_shifted);
    axis image; title('shifted PSF');
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

Java
----

Installation
++++++++++++

 ...

Java Interface
++++++++++++++

spline_analysis.spline_c
........................

Java Example
++++++++++++

example_name
............

...