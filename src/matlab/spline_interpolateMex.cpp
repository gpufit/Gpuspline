#include "src/spline.h"
#include <algorithm>
#include <mex.h>
#include <vector>

#if (defined _MSC_VER && _MSC_VER <= 1800)
#define PRINT_MSG _snprintf_s
#else
#define PRINT_MSG std::snprintf
#define PRINT_MSG std::snprintf
#endif

#ifdef SPLINE_DOUBLE
#define MX_REAL mxDOUBLE_CLASS
#define TOLERANCE_PRECISION_MESSAGE()\
    mexErrMsgIdAndTxt("Gpuspline:Mex", "tolerance is not a double");
#else
#define MX_REAL mxSINGLE_CLASS
#define TOLERANCE_PRECISION_MESSAGE()\
    mexErrMsgIdAndTxt("Gpuspline:Mex", "tolerance is not a single");
#endif // GPUSPLINE_DOUBLE

void mexFunction(
    int          nlhs,
    mxArray      *plhs[],
    int          nrhs,
    mxArray const *prhs[])
{
    // expected number of input arguments
    int const min_nrhs = 3;
    int const max_nrhs = 5;
    // expected number of output arguments
    int const expected_nlhs = 1;

    if (nrhs < min_nrhs || nrhs > max_nrhs)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "wrong number of input arguments.");
        mexErrMsgIdAndTxt("spline_interpolate:Mex", msg);
    }

    if (nlhs != expected_nlhs)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "%d output arguments required.", expected_nlhs);
        mexErrMsgIdAndTxt("spline_interpolate:Mex", msg);
    }

    // input parameters

    // input data
    REAL * data = (REAL*)mxGetPr(prhs[0]);

    // dimensions of input data
    int n_dimensions = (int)*mxGetPr(prhs[1]);

    if (n_dimensions < 1 || n_dimensions > 3)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "wrong number of dimensions.");
        mexErrMsgIdAndTxt("spline_interpolate:Mex", msg);
    }

    int const expected_nrhs = 2 + n_dimensions;
    if (nrhs != expected_nrhs)
    {
        char msg[50];
        PRINT_MSG(msg, 50, "%d input arguments required.", expected_nrhs);
        mexErrMsgIdAndTxt("spline_interpolate:Mex", msg);
    }

    std::vector<std::size_t> data_dimensions(n_dimensions);
    for (std::size_t i = 0; i < data_dimensions.size(); i++)
    {
        data_dimensions[i] = mxGetDimensions(prhs[0])[i];
    }

    // x values
    REAL * x_values = (REAL*)mxGetPr(prhs[2]);
    std::size_t new_size_x
        = std::max(mxGetDimensions(prhs[2])[0], mxGetDimensions(prhs[2])[1]);

    // y values
    REAL * y_values;
    std::size_t new_size_y = 1;
    if (n_dimensions > 1)
    {
        y_values = (REAL*)mxGetPr(prhs[3]);
        new_size_y
            = std::max(mxGetDimensions(prhs[3])[0], mxGetDimensions(prhs[3])[1]);
    }

    // z values
    REAL * z_values;
    std::size_t new_size_z = 1;
    if (n_dimensions > 2)
    {
        z_values = (REAL*)mxGetPr(prhs[4]);
        new_size_z 
            = std::max(mxGetDimensions(prhs[4])[0], mxGetDimensions(prhs[4])[1]);
    }

    // number of points after interpolation
    std::size_t const new_n_points
        = new_size_x * new_size_y * new_size_z;

    // output parameter: preallocate interpolated_data array
    REAL * interpolated_data;
    mxArray * mx_interpolated_data;
    mx_interpolated_data = mxCreateNumericMatrix(new_n_points, 1, MX_REAL, mxREAL);
    interpolated_data = (REAL*)mxGetData(mx_interpolated_data);
    plhs[0] = mx_interpolated_data;

    switch (n_dimensions)
    {
    case 1:
    {
        interpolate_1d(
            data,
            data_dimensions[0],
            new_size_x,
            x_values,
            interpolated_data);
        
        break;
    }
    case 2:
    {
        interpolate_2d(
            data,
            data_dimensions[0],
            data_dimensions[1],
            new_size_x,
            new_size_y,
            x_values,
            y_values,
            interpolated_data);

        break;
    }
    case 3:
    {
        interpolate_3d(
            data,
            data_dimensions[0],
            data_dimensions[1],
            data_dimensions[2],
            new_size_x,
            new_size_y,
            new_size_z,
            x_values,
            y_values,
            z_values,
            interpolated_data);
        
        break;
    }
    }
}