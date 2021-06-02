.. _api-description:

===============
API description
===============

The Gpuspline source code compiles to a dynamic-link library (DLL), providing C interfaces.
In the sections below, the C interfaces and their arguments are described in detail.

.. _c-interface:

C Interfaces
------------

The C interfaces are defined in the header file: spline.h_.

calculate_coefficients_1d()
+++++++++++++++++++++++++++

This function calculates 1D cubic spline coefficients for intervals between data points of the
input data. Every interval is represented by 4 coefficients.

.. code-block:: cpp

    int calculate_coefficients_1d(
        REAL * data,
        size_t data_size_x,
        REAL * coefficients);
        
.. _api-c-1d-input-parameters:

Description of input parameters
...............................

:data: Pointer to data values

    :type: REAL *
    :length: data_size_x

:data_size_x: Number of data points

    :type: size_t

.. _api-c-1d-output-parameters:

Description of output parameters
................................

:coefficients: calculated spline coefficients

    :type: REAL *
    :length: 4 * (data_size_x - 1)

calculate_coefficients_2d()
+++++++++++++++++++++++++++

This function calculates 2D cubic spline coefficients for intervals between data points of the
input data. Every interval is represented by 16 coefficients.

.. code-block:: cpp
    
    int calculate_coefficients_2d(
        REAL * data,
        size_t data_size_x,
        size_t data_size_y,
        REAL * coefficients);

.. _api-c-2d-input-parameters:

Description of input parameters
...............................

:data: Pointer to data values

    :type: REAL *
    :length: data_size_x * data_size_y

:data_size_x: Data dimension x

    :type: size_t

:data_size_y: Data dimension y

    :type: size_t

.. _api-c-2d-output-parameters:

Description of output parameters
................................

:coefficients: calculated spline coefficients

    :type: REAL *
    :length: 16 * (data_size_x - 1) * (data_size_y - 1)
        
calculate_coefficients_3d()
+++++++++++++++++++++++++++

This function calculates 3D cubic spline coefficients for intervals between data points of the
input data. Every interval is represented by 64 coefficients.

.. code-block:: cpp

    int calculate_coefficients_3d(
        REAL * data,
        size_t data_size_x,
        size_t data_size_y,
        size_t data_size_z,
        REAL * coefficients);

.. _api-c-3d-input-parameters:

Description of input parameters
...............................

:data: Pointer to data values

    :type: REAL *
    :length: data_size_x * data_size_y * data_size_z

:data_size_x: Data dimension x

    :type: size_t

:data_size_y: Data dimension y

    :type: size_t

:data_size_z: Data dimension z

    :type: size_t

.. _api-c-3d-output-parameters:

Description of output parameters
................................

:coefficients: calculated spline coefficients

    :type: REAL *
    :length: 64 * (data_size_x - 1) * (data_size_y - 1) * (data_size_z - 1)

interpolate_1d()
++++++++++++++++

This function performs a 1D data interpolation based on the cubic spline interpolation method.

.. code-block:: cpp

    int interpolate_1d(
        REAL * data,
        size_t data_size_x,
        size_t new_size_x,
        REAL * x_values,
        REAL * interpolated_data);

.. _api-i-1d-input-parameters:

Description of input parameters
...............................

:data: Pointer to data values

    :type: REAL *
    :length: data_size_x

:data_size_x: number of input data points

    :type: size_t

:new_size_x: number of output data points

    :type: size_t

:x_values: Pointer to independent variable values

    :type: REAL *
    :length: new_size_x

.. _api-i-1d-output-parameters:

Description of output parameters
................................

:interpolated_data: Pointer to output data values

    :type: REAL *
    :length: new_size_x

interpolate_2d()
++++++++++++++++

This function performs a 2D data interpolation based on the cubic spline interpolation method.

.. code-block:: cpp

    int interpolate_2d(
        REAL * data,
        size_t data_size_x,
        size_t data_size_y,
        size_t new_size_x,
        size_t new_size_y,
        REAL * x_values,
        REAL * y_values,
        REAL * interpolated_data);

.. _api-i-2d-input-parameters:

Description of input parameters
...............................

:data: Pointer to data values

    :type: REAL *
    :length: data_size_x * data_size_y

:data_size_x: Input data dimension x

    :type: size_t

:data_size_y: Input data dimension y

    :type: size_t

:new_size_x: Output data dimension x

    :type: size_t

:new_size_y: Output data dimension y

    :type: size_t

:x_values: Pointer to independent variable x values

    :type: REAL *
    :length: new_size_x

:y_values: Pointer to independent variable y values

    :type: REAL *
    :length: new_size_y

.. _api-i-2d-output-parameters:

Description of output parameters
................................

:interpolated_data: Pointer to output data values

    :type: REAL *
    :length: new_size_x * new_size_y

interpolate_3d()
++++++++++++++++

This function performs a 3D data interpolation based on the cubic spline interpolation method.

.. code-block:: cpp

    int interpolate_3d(
        REAL * data,
        size_t data_size_x,
        size_t data_size_y,
        size_t data_size_z,
        size_t new_size_x,
        size_t new_size_y,
        size_t new_size_z,
        REAL * x_values,
        REAL * y_values,
        REAL * z_values,
        REAL * interpolated_data);

.. _api-i-3d-input-parameters:

Description of input parameters
...............................

:data: Pointer to data values

    :type: REAL *
    :length: data_size_x * data_size_y * data_size_z

:data_size_x: Input data dimension x

    :type: size_t

:data_size_y: Input data dimension y

    :type: size_t

:data_size_z: Input data dimension z

    :type: size_t

:new_size_x: Output data dimension x

    :type: size_t

:new_size_y: Output data dimension y

    :type: size_t

:new_size_z: Output data dimension z

    :type: size_t

:x_values: Pointer to independent variable x values

    :type: REAL *
    :length: new_size_x

:y_values: Pointer to independent variable y values

    :type: REAL *
    :length: new_size_y

:z_values: Pointer to independent variable z values

    :type: REAL *
    :length: new_size_z

.. _api-i-3d-output-parameters:

Description of output parameters
................................

:interpolated_data: Pointer to output data values

    :type: REAL *
    :length: new_size_x * new_size_y * new_size_z

calculate_values_1d()
+++++++++++++++++++++

This function calculates 1D function values based on provided spline coefficients and independent
variable values.

.. code-block:: cpp

    int calculate_values_1d(
        REAL * coefficients,
        size_t const n_intervals_x,
        size_t const values_size_x,
        REAL * x_values,
        REAL * spline_values);

.. _api-v-1d-input-parameters:

Description of input parameters
...............................

:coefficients: Pointer to spline coefficients

    :type: REAL *
    :length: 4 * n_intervals_x

:n_intervals_x: Number of spline intervals

    :type: size_t

:values_size_x: Number of output data points

    :type: size_t

:x_values: Pointer to independent variable values

    :type: REAL *
    :length: values_size_x

.. _api-v-1d-output-parameters:

Description of output parameters
................................

:spline_values: Pointer to output data values

    :type: REAL *
    :length: values_size_x

calculate_values_2d()
+++++++++++++++++++++

This function calculates function values based on provided spline coefficients and independent
variable values.

.. code-block:: cpp

    int calculate_values_2d(
        REAL * coefficients,
        size_t const n_intervals_x,
        size_t const n_intervals_y,
        size_t const values_size_x,
        size_t const values_size_y,
        REAL * x_values,
        REAL * y_values,
        REAL * spline_values);

.. _api-v-2d-input-parameters:

Description of input parameters
...............................

:coefficients: Pointer to spline coefficients

    :type: REAL *
    :length: 16 * n_intervals_x * n_intervals_y

:n_intervals_x: Number of spline intervals in x

    :type: size_t

:n_intervals_y: Number of spline intervals in y

    :type: size_t

:values_size_x: Output data dimension x

    :type: size_t

:values_size_y: Output data dimension y

    :type: size_t

:x_values: Pointer to independent variable x values

    :type: REAL *
    :length: values_size_x

:y_values: Pointer to independent variable y values

    :type: REAL *
    :length: values_size_y

.. _api-v-2d-output-parameters:

Description of output parameters
................................

:spline_values: Pointer to output data values

    :type: REAL *
    :length: values_size_x * values_size_y

calculate_values_3d()
+++++++++++++++++++++

This function calculates function values based on provided spline coefficients and independent
variable values.

.. code-block:: cpp

    int calculate_values_3d(
        REAL * coefficients,
        size_t const n_intervals_x,
        size_t const n_intervals_y,
        size_t const n_intervals_z,
        size_t const values_size_x,
        size_t const values_size_y,
        size_t const values_size_z,
        REAL * x_values,
        REAL * y_values,
        REAL * z_values,
        REAL * spline_values);

.. _api-v-3d-input-parameters:

Description of input parameters
...............................

:coefficients: Pointer to spline coefficients

    :type: REAL *
    :length: 64 * n_intervals_x * n_intervals_y * n_intervals_z

:n_intervals_x: Number of spline intervals in x

    :type: size_t

:n_intervals_y: Number of spline intervals in y

    :type: size_t

:n_intervals_z: Number of spline intervals in z

    :type: size_t

:values_size_x: Output data dimension x

    :type: size_t

:values_size_y: Output data dimension y

    :type: size_t

:values_size_y: Output data dimension z

    :type: size_t

:x_values: Pointer to independent variable x values

    :type: REAL *
    :length: values_size_x

:y_values: Pointer to independent variable y values

    :type: REAL *
    :length: values_size_y

:z_values: Pointer to independent variable z values

    :type: REAL *
    :length: values_size_z

.. _api-v-3d-output-parameters:

Description of output parameters
................................

:spline_values: Pointer to output data values

    :type: REAL *
    :length: values_size_x * values_size_y * values_size_z

calculate_coefficients_1d_portable()
++++++++++++++++++++++++++++++++++++

This function is a simple wrapper around the :code:`calculate_coefficients_1d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int calculate_coefficients_1d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *calculate_coefficients_1d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Data

        :type: REAL *

    :argv[1]: Number of data points

        :type: size_t *

    :argv[2]: Spline coefficients

        :type: REAL *

calculate_coefficients_2d_portable()
++++++++++++++++++++++++++++++++++++

This function is a simple wrapper around the :code:`calculate_coefficients_2d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int calculate_coefficients_2d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *calculate_coefficients_2d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Data

        :type: REAL *

    :argv[1]: Data dimension x

        :type: size_t *

    :argv[2]: Data dimension y

        :type: size_t *

    :argv[3]: Spline coefficients

        :type: REAL *

calculate_coefficients_3d_portable()
++++++++++++++++++++++++++++++++++++

This function is a simple wrapper around the :code:`calculate_coefficients_3d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int calculate_coefficients_3d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *calculate_coefficients_3d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Data

        :type: REAL *

    :argv[1]: Data dimension x

        :type: size_t *

    :argv[2]: Data dimension y

        :type: size_t *

    :argv[3]: Data dimension z

        :type: size_t *

    :argv[4]: Spline coefficients

        :type: REAL *

interpolate_1d_portable()
+++++++++++++++++++++++++

This function is a simple wrapper around the :code:`interpolate_1d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int interpolate_1d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *interpolate_1d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Input data

        :type: REAL *

    :argv[1]: Input number of data points

        :type: size_t *

    :argv[2]: Output number of data points

        :type: size_t *

    :argv[3]: Independent variable values

        :type: REAL *

    :argv[4]: Output data

        :type: REAL *

interpolate_2d_portable()
+++++++++++++++++++++++++

This function is a simple wrapper around the :code:`interpolate_2d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int interpolate_2d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *interpolate_2d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Input data

        :type: REAL *

    :argv[1]: Input data dimension x

        :type: size_t *

    :argv[2]: Input data dimension y

        :type: size_t *

    :argv[3]: Output data dimension x

        :type: size_t *

    :argv[4]: Output data dimension y

        :type: size_t *

    :argv[5]: Independent variable x values

        :type: REAL *

    :argv[6]: Independent variable y values

        :type: REAL *

    :argv[7]: Output data

        :type: REAL *

interpolate_3d_portable()
+++++++++++++++++++++++++

This function is a simple wrapper around the :code:`interpolate_3d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int interpolate_3d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *interpolate_3d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Input data

        :type: REAL *

    :argv[1]: Input data dimension x

        :type: size_t *

    :argv[2]: Input data dimension y

        :type: size_t *

    :argv[3]: Input data dimension z

        :type: size_t *

    :argv[4]: Output data dimension x

        :type: size_t *

    :argv[5]: Output data dimension y

        :type: size_t *

    :argv[6]: Output data dimension z

        :type: size_t *

    :argv[7]: Independent variable x values

        :type: REAL *

    :argv[8]: Independent variable y values

        :type: REAL *

    :argv[9]: Independent variable z values

        :type: REAL *

    :argv[10]: Output data

        :type: REAL *

calculate_values_1d_portable()
++++++++++++++++++++++++++++++

This function is a simple wrapper around the :code:`calculate_values_1d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int calculate_values_1d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *calculate_values_1d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Spline coefficients

        :type: REAL *

    :argv[1]: Number of spline intervals

        :type: size_t *

    :argv[2]: Number of output data points

        :type: size_t *

    :argv[3]: Independent variable values

        :type: REAL *

    :argv[4]: Output data values

        :type: REAL *

calculate_values_2d_portable()
++++++++++++++++++++++++++++++

This function is a simple wrapper around the :code:`calculate_values_2d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int calculate_values_2d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *calculate_values_2d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Spline coefficients

        :type: REAL *

    :argv[1]: Number of spline intervals in x

        :type: size_t *

    :argv[2]: Number of spline intervals in y

        :type: size_t *

    :argv[3]: Output data dimension x

        :type: size_t *

    :argv[4]: Output data dimension y

        :type: size_t *

    :argv[5]: Independent variable x values

        :type: REAL *

    :argv[6]: Independent variable y values

        :type: REAL *

    :argv[7]: Output data values

        :type: REAL *

calculate_values_3d_portable()
++++++++++++++++++++++++++++++

This function is a simple wrapper around the :code:`calculate_values_3d()` function,
providing an alternative means of passing the function parameters.

.. code-block:: cpp

    int calculate_values_3d_portable(int argc, void *argv[]);

Description of parameters
.........................

:argc: The length of the argv pointer array

:argv: Array of pointers to *calculate_values_3d* parameters, as defined above.
    For reference, the type of each element of the *argv* array is listed below.

    :argv[0]: Spline coefficients

        :type: REAL *

    :argv[1]: Number of spline intervals in x

        :type: size_t *

    :argv[2]: Number of spline intervals in y

        :type: size_t *

    :argv[3]: Number of spline intervals in z

        :type: size_t *

    :argv[4]: Output data dimension x

        :type: size_t *

    :argv[5]: Output data dimension y

        :type: size_t *

    :argv[6]: Output data dimension z

        :type: size_t *

    :argv[7]: Independent variable x values

        :type: REAL *

    :argv[8]: Independent variable y values

        :type: REAL *

    :argv[9]: Independent variable z values

        :type: REAL *

    :argv[10]: Output data values

        :type: REAL *