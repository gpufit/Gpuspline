#include "spline.h"
#include "spline_classes.h"

int calculate_coefficients_2d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    REAL* coefficients)
{
    Spline2D spline_2d(data_size_x, data_size_y);
    spline_2d.initialize(data);
    spline_2d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_2d_portable(int argc, void* argv[])
{
    return calculate_coefficients_2d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        (REAL*)argv[3]);
}

int convert_csaps_coefficients_2d(
    REAL* csaps_coefficients,
    std::size_t n_spline_intervals_x,
    std::size_t n_spline_intervals_y,
    REAL* grid_spacing_array,
    REAL* reordered_coefficients)
{

    Spline2D spline_2d(n_spline_intervals_x + 1, n_spline_intervals_y + 1);
    spline_2d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals_x, n_spline_intervals_y, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_2d_portable(int argc, void* argv[])
{
    return convert_csaps_coefficients_2d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        (REAL*)argv[3],
        (REAL*)argv[4]);
}

int interpolate_2d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t new_size_x,
    std::size_t new_size_y,
    REAL* x_values,
    REAL* y_values,
    REAL* interpolated_data)
{
    Spline2D spline_2d(data_size_x, data_size_y);

    spline_2d.interpolate(
        interpolated_data,
        data,
        x_values, y_values,
        new_size_x, new_size_y);

    return 0;
}

int interpolate_2d_portable(int argc, void* argv[])
{
    return interpolate_2d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        (REAL*)argv[5],
        (REAL*)argv[6],
        (REAL*)argv[7]);
}

int calculate_values_2d(
    REAL* coefficients,
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const values_size_x,
    std::size_t const values_size_y,
    REAL* x_values,
    REAL* y_values,
    REAL* spline_values)
{
    Spline2D spline_2d(n_intervals_x, n_intervals_y, coefficients);
    spline_2d.calculate_values(
        spline_values,
        x_values,
        y_values,
        values_size_x,
        values_size_y);

    return 0;
}

int calculate_values_2d_portable(int argc, void* argv[])
{
    return calculate_values_2d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        (REAL*)argv[5],
        (REAL*)argv[6],
        (REAL*)argv[7]);
}

