#include "spline.h"
#include "spline_classes.h"

int calculate_coefficients_1d(
    REAL* data,
    std::size_t data_size_x,
    REAL* coefficients)
{
    Spline1D spline_1d(data_size_x);
    spline_1d.initialize(data);
    spline_1d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_1d_portable(int argc, void* argv[])
{
    return calculate_coefficients_1d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        (REAL*)argv[2]);
}

int convert_csaps_coefficients_1d(
    REAL* csaps_coefficients,
    std::size_t n_spline_intervals,
    REAL* grid_spacing_array,
    REAL* reordered_coefficients)
{

    Spline1D spline_1d(n_spline_intervals + 1);
    spline_1d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_1d_portable(int argc, void* argv[])
{
    return convert_csaps_coefficients_1d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        (REAL*)argv[2],
        (REAL*)argv[3]);
}

int interpolate_1d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t new_size_x,
    REAL* x_values,
    REAL* interpolated_data)
{
    Spline1D spline_1d(data_size_x);

    spline_1d.interpolate(
        interpolated_data,
        data,
        x_values,
        new_size_x);

    return 0;
}

int interpolate_1d_portable(int argc, void* argv[])
{
    return interpolate_1d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        (REAL*)argv[3],
        (REAL*)argv[4]);
}

int calculate_values_1d(
    REAL* coefficients,
    std::size_t const n_intervals_x,
    std::size_t const values_size_x,
    REAL* x_values,
    REAL* spline_values)
{
    Spline1D spline_1d(n_intervals_x, coefficients);
    spline_1d.calculate_values(spline_values, x_values, values_size_x);

    return 0;
}

int calculate_values_1d_portable(int argc, void* argv[])
{
    return calculate_values_1d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        (REAL*)argv[3],
        (REAL*)argv[4]);
}

