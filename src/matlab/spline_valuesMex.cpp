#include "src/spline.h"
#include <mex.h>
#include <vector>
#include <algorithm>
#include <cmath>

/*
* Given some data, it calculates the coefficients
*/

#if (defined _MSC_VER && _MSC_VER <= 1800)
#define PRINT_MSG _snprintf_s
#else
#define PRINT_MSG std::snprintf
#endif

#ifdef GPUSPLINE_DOUBLE
#define MX_REAL mxDOUBLE_CLASS
#define TOLERANCE_PRECISION_MESSAGE()\
    mexErrMsgIdAndTxt("Gpuspline:Mex", "tolerance is not a double");
#else
#define MX_REAL mxSINGLE_CLASS
#define TOLERANCE_PRECISION_MESSAGE()\
    mexErrMsgIdAndTxt("Gpuspline:Mex", "tolerance is not a single");
#endif // GPUFIT_DOUBLE

// entry point for Matlab
void mexFunction(
    int          nlhs,
    mxArray      *plhs[],
    int          nrhs,
    mxArray const *prhs[])
{
    // expected number of input arguments
    int const min_nrhs = 6;
    int const max_nrhs = 8;
    // expected number of output arguments
    int const expected_nlhs = 1;

    if (nrhs < min_nrhs || nrhs > max_nrhs)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "wrong number of input arguments.");
        mexErrMsgIdAndTxt("spline_values:Mex", msg);
    }

    if (nlhs != expected_nlhs)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "%d output arguments required.", expected_nlhs);
        mexErrMsgIdAndTxt("spline_values:Mex", msg);
    }

    REAL * coefficients = (REAL*)mxGetPr(prhs[0]);
    std::size_t const n_dimensions = (std::size_t)*mxGetPr(prhs[1]);

    if (n_dimensions < 1 || n_dimensions > 3)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "wrong number of dimensions.");
        mexErrMsgIdAndTxt("spline_values:Mex", msg);
    }

    int const expected_nrhs = min_nrhs + static_cast<int>(n_dimensions) - 1;
    if (nrhs != expected_nrhs)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "%d input arguments required.", expected_nrhs);
        mexErrMsgIdAndTxt("spline_values:Mex", msg);
    }

    // constants
    std::size_t const n_coefficients_per_point = static_cast<std::size_t>(std::pow(4, n_dimensions));

    std::size_t const n_intervals_x = (std::size_t)*mxGetPr(prhs[2]); // note this only works if the type is double of prhs[2-4]!
    std::size_t const n_intervals_y = (std::size_t)*mxGetPr(prhs[3]);
    std::size_t const n_intervals_z = (std::size_t)*mxGetPr(prhs[4]);

    std::size_t n_intervals;
    if (n_dimensions == 1)
        n_intervals = n_intervals_x;
    else if (n_dimensions == 2)
        n_intervals = n_intervals_x *n_intervals_y;
    else if (n_dimensions == 3)
        n_intervals = n_intervals_x * n_intervals_y * n_intervals_z;

    // x values
    REAL * x_values = (REAL*)mxGetPr(prhs[5]);
    std::size_t values_size_x
        = std::max(mxGetDimensions(prhs[5])[0], mxGetDimensions(prhs[5])[1]);

    // y values
    REAL * y_values;
    std::size_t values_size_y = 1;
    if (n_dimensions > 1)
    {
        y_values = (REAL*)mxGetPr(prhs[6]);
        values_size_y
            = std::max(mxGetDimensions(prhs[6])[0], mxGetDimensions(prhs[6])[1]);
    }

    // z values
    REAL * z_values;
    std::size_t values_size_z = 1;
    if (n_dimensions > 2)
    {
        z_values = (REAL*)mxGetPr(prhs[7]);
        values_size_z
            = std::max(mxGetDimensions(prhs[7])[0], mxGetDimensions(prhs[7])[1]);
    }

    // size of spline function
    std::size_t const values_size
        = values_size_x * values_size_y * values_size_z;

    // output parameter: preallocate spline_values array
    REAL * spline_values;
    mxArray * mx_spline_values;
    mx_spline_values = mxCreateNumericMatrix(values_size, 1, MX_REAL, mxREAL);
    spline_values = (REAL*)mxGetData(mx_spline_values);
    plhs[0] = mx_spline_values;

    // TODO as above this could be detected automatically
    switch (n_dimensions)
    {
    case 1:
    {        
        calculate_values_1d(
            coefficients,
            n_intervals_x,
            values_size_x,
            x_values,
            spline_values);

        break;
    }
    case 2:
    {
        calculate_values_2d(
            coefficients,
            n_intervals_x,
            n_intervals_y,
            values_size_x,
            values_size_y,
            x_values,
            y_values,
            spline_values);

        break;
    }
    case 3:
    {
        calculate_values_3d(
            coefficients,
            n_intervals_x,
            n_intervals_y,
            n_intervals_z,
            values_size_x,
            values_size_y,
            values_size_z,
            x_values,
            y_values,
            z_values,
            spline_values);

        break;
    }
    default:
    {
        char msg[50];
        PRINT_MSG(msg, 50, "Wrong number of dimensions!");
        mexErrMsgIdAndTxt("Gpuspline:Mex", msg);
        break;
    }
    }
}