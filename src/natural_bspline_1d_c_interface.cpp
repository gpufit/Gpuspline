#include "spline.h"
#include "natural_bspline_1d.h"

int calculate_coefficients_natural_bspline_1d(
    int num_data_points,
    REAL* data,
    REAL* coefficients_out)
{
    int num_control_points = num_data_points + 2;

   // Create spline and compute coefficients
    Natural_BSpline_1D spline(num_data_points, data);

    // Copy coefficients
    const REAL* coeffs = spline.coefficients();

    for (int i = 0; i < num_control_points; i++) {
        coefficients_out[i] = coeffs[i];
    }

    return 0;
}


int calculate_coefficients_natural_bspline_1d_portable(int argc, void* argv[])
{
    return calculate_coefficients_natural_bspline_1d(
        *((int*) argv[0]),
        (REAL*) argv[1],
        (REAL*) argv[2]);
}


int calculate_values_natural_bspline_1d(
    int num_data_points,
    const REAL* coefficients,
    int n_input_coords,
    const REAL* input_coords,
    int flag_fast_evaluate, 
    REAL* output_values)
{
    int num_control_points = num_data_points + 2;

    // Create spline
    Natural_BSpline_1D spline(num_data_points, num_control_points, coefficients);

    // Set fast evaluation
    spline.set_fast_evaluation(flag_fast_evaluate != 0);

    // Evaluate
    spline.calculate_values(input_coords, n_input_coords, output_values);

    return 0;
}


int calculate_values_natural_bspline_1d_portable(int argc, void* argv[])
{
    return calculate_values_natural_bspline_1d(
        *((int*) argv[0]),
        (REAL*) argv[1],
        *((int*) argv[2]),
        (REAL*) argv[3],
        *((int*)argv[4]),
        (REAL*) argv[5]);
}
