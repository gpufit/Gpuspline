#include "src/spline.h"
#include <mex.h>
#include <cstring>
#include <string>
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
#endif // GPUSPLINE_DOUBLE

// entry point for Matlab
void mexFunction(
    int          nlhs,
    mxArray      *plhs[],
    int          nrhs,
    mxArray const *prhs[])
{
    // expected number of input arguments
    int const expected_nrhs = 5;
    // expected number of output arguments
    int const expected_nlhs = 1;

    if (nrhs != expected_nrhs)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "%d input arguments required.", expected_nrhs);
        mexErrMsgIdAndTxt("Gpuspline:Mex", msg);
    }

    if (nlhs != expected_nlhs)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "%d output arguments required.", expected_nlhs);
        mexErrMsgIdAndTxt("Gpuspline:Mex", msg);
    }

    // TODO data is a matlab array, we can read the dimensions inside here, no need to explicitly write them down?)
    // input parameters (data, nx, ny, nz, ndims)
    REAL * data = (REAL*)mxGetPr(prhs[0]);
    std::size_t n_points_x = (std::size_t)*mxGetPr(prhs[1]);
    std::size_t n_points_y = (std::size_t)*mxGetPr(prhs[2]);
    std::size_t n_points_z = (std::size_t)*mxGetPr(prhs[3]);
    std::size_t n_dimensions = (std::size_t)*mxGetPr(prhs[4]);

    // constants
    std::size_t const n_points = n_points_x * n_points_y * n_points_z;
    std::size_t const n_coefficients_per_point = static_cast<std::size_t>(std::pow(4, n_dimensions));
    std::size_t const n_spline_points_x = n_points_x - 1;
    std::size_t const n_spline_points_y = n_points_y - 1;
    std::size_t const n_spline_points_z = n_points_z - 1;

    std::size_t n_spline_points = 0;
    if (n_dimensions == 1)
        n_spline_points = n_spline_points_x;
    else if (n_dimensions == 2)
        n_spline_points = n_spline_points_x * n_spline_points_y;
    else if (n_dimensions == 3)
        n_spline_points = n_spline_points_x * n_spline_points_y * n_spline_points_z;

    // output parameter: preallocate coefficients array
    REAL * coefficients;
    mxArray * mx_coefficients;
    mx_coefficients = mxCreateNumericMatrix(1, n_spline_points * n_coefficients_per_point, MX_REAL, mxREAL);
    coefficients = (REAL*)mxGetData(mx_coefficients);
    plhs[0] = mx_coefficients;

    // TODO as above this could be detected automatically
    switch (n_dimensions)
    {
    case 1:
    {
        calculate_coefficients_1d(
            data,
            n_points,
            coefficients);

        break;
    }
    case 2:
    {
        calculate_coefficients_2d(
            data,
            n_points_x,
            n_points_y,
            coefficients);

        break;
    }
    case 3:
    {
        calculate_coefficients_3d(
            data,
            n_points_x,
            n_points_y,
            n_points_z,
            coefficients);

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