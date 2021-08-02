"""
Python binding for Gpuspline, a library for the calculation of multidimensional cubic splines
See https://github.com/gpufit/Gpuspline

The binding is based on ctypes.
See https://docs.python.org/3.5/library/ctypes.html, http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
"""

import os
from ctypes import cdll, POINTER, c_float, c_size_t
import numpy as np

# define library loader (actual loading is lazy)
package_dir = os.path.dirname(os.path.realpath(__file__))

if os.name == 'nt':
    lib_path = os.path.join(package_dir, 'splines.dll')  # library name on Windows
elif os.name == 'posix':
    lib_path = os.path.join(package_dir, 'libsplines.so')  # library name on Unix
else:
    raise RuntimeError('OS {} not supported by pyGpuspline.'.format(os.name))

lib = cdll.LoadLibrary(lib_path)


def convert_ndarray(x):
    """
    We may get lists instead of NumPy arrays from outside. Just silently convert them to (single-valued) NumPy arrays.
    :param x: List or NumPy array
    :return: If x was a NumPy array just returns itself, otherwise a NumPy array initialized with x.
    """
    if not isinstance(x, np.ndarray):
        x = np.array(x, dtype=np.float32)
    return x


def spline_coefficients(data):
    """
    Calculates spline coefficients from given 1-3D data
    """
    data = convert_ndarray(data)
    if data.dtype != np.float32:
        raise RuntimeError('Type of data is not single!')
    size = data.shape

    number_dimensions = len(size)
    number_coefficients_per_interval = pow(4, number_dimensions)

    # preallocate coefficients array
    intervals = [n - 1 for n in size]
    coefficients = np.zeros([number_coefficients_per_interval] + intervals, np.float32)

    # call to calculate_coefficients_xd in the library
    if number_dimensions == 1:
        status = lib.calculate_coefficients_1d(
            data.ctypes.data_as(POINTER(c_float)),
            c_size_t(size[0]),
            coefficients.ctypes.data_as(POINTER(c_float))
        )
    elif number_dimensions == 2:
        status = lib.calculate_coefficients_2d(
            data.ctypes.data_as(POINTER(c_float)),
            c_size_t(size[1]),
            c_size_t(size[0]),
            coefficients.ctypes.data_as(POINTER(c_float))
        )
    elif number_dimensions == 3:
        status = lib.calculate_coefficients_3d(
            data.ctypes.data_as(POINTER(c_float)),
            c_size_t(size[2]),
            c_size_t(size[1]),
            c_size_t(size[0]),
            coefficients.ctypes.data_as(POINTER(c_float))
        )

    # check status
    if status != 0:
        raise RuntimeError('Call to Gpuspline library failed, exit code {}'.format(status))

    # return coefficients
    return coefficients


def spline_interpolate(data, x, y=None, z=None):
    """
    Interpolates 1-3D data at positions x, y, z.

    Note: Internally calls spline_coefficients first to create a spline representation
    and then spline_value to obtain interpolated values. If repeated interpolation
    of the same data set is needed, it's better to call spline_coefficients only once
    and then use spline_value instead of this function.
    """
    data = convert_ndarray(data)
    if data.dtype != np.float32:
        raise RuntimeError('Type of data is not single!')
    size = data.shape

    x = convert_ndarray(x)
    if x.dtype != np.float32:
        raise RuntimeError('Type of x values is not single!')

    new_size_x = x.size
    if y is None:
        # 1D case
        number_dimensions = 1
        new_size_y = 1
        new_size_z = 1
    elif z is None:
        # 2D case
        number_dimensions = 2
        y = convert_ndarray(y)
        if y.dtype != np.float32:
            raise RuntimeError('Type of y values is not single!')
        new_size_y = y.size
        new_size_z = 1
    else:
        # 3D case
        number_dimensions = 3
        z = convert_ndarray(z)
        if z.dtype != np.float32:
            raise RuntimeError('Type of z values is not single!')
        new_size_y = y.size
        new_size_z = z.size
    new_number_points = new_size_x * new_size_y * new_size_z

    # preallocate interpolated_data array
    interpolated_size = [new_size_x, new_size_y, new_size_z][:number_dimensions]
    interpolated_data = np.zeros(interpolated_size, np.float32)

    # call to interpolate_xd in the library
    if number_dimensions == 1:
        status = lib.interpolate_1d(
            data.ctypes.data_as(POINTER(c_float)),
            c_size_t(size[0]),
            c_size_t(new_size_x),
            x.ctypes.data_as(POINTER(c_float)),
            interpolated_data.ctypes.data_as(POINTER(c_float))
        )
    elif number_dimensions == 2:
        status = lib.interpolate_2d(
            data.ctypes.data_as(POINTER(c_float)),
            c_size_t(size[1]),
            c_size_t(size[0]),
            c_size_t(new_size_y),
            c_size_t(new_size_x),
            y.ctypes.data_as(POINTER(c_float)),
            x.ctypes.data_as(POINTER(c_float)),
            interpolated_data.ctypes.data_as(POINTER(c_float))
        )
    elif number_dimensions == 3:
        status = lib.interpolate_3d(
            data.ctypes.data_as(POINTER(c_float)),
            c_size_t(size[2]),
            c_size_t(size[1]),
            c_size_t(size[0]),
            c_size_t(new_size_z),
            c_size_t(new_size_y),
            c_size_t(new_size_x),
            z.ctypes.data_as(POINTER(c_float)),
            y.ctypes.data_as(POINTER(c_float)),
            x.ctypes.data_as(POINTER(c_float)),
            interpolated_data.ctypes.data_as(POINTER(c_float))
        )

    # check status
    if status != 0:
        raise RuntimeError('Call to Gpuspline library failed, exit code {}'.format(status))

    # return interpolated data
    return interpolated_data


def spline_values(coefficients, x, y=None, z=None):
    """
    Calculates spline values given a spline representation and x, y and z values.

    y and z parameter is optional (only if you want 2D/3D values)

    Note: x, y, z are a range of values along their coordinate axis and the total
    number is the cartesian product of them all!
    """
    # if not numpy arrays, convert to them
    coefficients = convert_ndarray(coefficients)
    if coefficients.dtype != np.float32:
        raise RuntimeError('Type of coefficients is not single!')

    x = convert_ndarray(x)
    if x.dtype != np.float32:
        raise RuntimeError('Type of x values is not single!')

    values_size_x = x.size
    if y is None:
        # 1D case
        number_dimensions = 1
        values_size_y = 1
        values_size_z = 1
        values_shape = (values_size_x,)
    elif z is None:
        # 2D case
        number_dimensions = 2
        y = convert_ndarray(y)
        if y.dtype != np.float32:
            raise RuntimeError('Type of y values is not single!')
        values_size_y = y.size
        values_size_z = 1
        values_shape = (values_size_x,values_size_y)
    else:
        # 3D case
        number_dimensions = 3
        z = convert_ndarray(z)
        if z.dtype != np.float32:
            raise RuntimeError('Type of z values is not single!')
        values_size_y = y.size
        values_size_z = z.size
        values_shape = (values_size_x, values_size_y, values_size_z)

    # preallocate spline_values array
    values = np.zeros(values_shape, np.float32)

    # call to calculate_values_xd in the library
    if number_dimensions == 1:
        status = lib.calculate_values_1d_func(
            coefficients.ctypes.data_as(POINTER(c_float)),
            c_size_t(coefficients.shape[1]),
            c_size_t(values_size_x),
            x.ctypes.data_as(POINTER(c_float)),
            values.ctypes.data_as(POINTER(c_float))
        )
    elif number_dimensions == 2:
        status = lib.calculate_values_2d(
            coefficients.ctypes.data_as(POINTER(c_float)),
            c_size_t(coefficients.shape[2]),
            c_size_t(coefficients.shape[1]),
            c_size_t(values_size_y),
            c_size_t(values_size_x),
            y.ctypes.data_as(POINTER(c_float)),
            x.ctypes.data_as(POINTER(c_float)),
            values.ctypes.data_as(POINTER(c_float))
        )
    elif number_dimensions == 3:
        status = lib.calculate_values_3d(
            coefficients.ctypes.data_as(POINTER(c_float)),
            c_size_t(coefficients.shape[3]),
            c_size_t(coefficients.shape[2]),
            c_size_t(coefficients.shape[1]),
            c_size_t(values_size_z),
            c_size_t(values_size_y),
            c_size_t(values_size_x),
            z.ctypes.data_as(POINTER(c_float)),
            y.ctypes.data_as(POINTER(c_float)),
            x.ctypes.data_as(POINTER(c_float)),
            values.ctypes.data_as(POINTER(c_float))
        )

    # check status
    if status != 0:
        raise RuntimeError('Call to Gpuspline library failed, exit code {}'.format(status))

    # return spline values
    return values