#include "spline.h"
#include "spline_classes.h"

int calculate_coefficients_4d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t data_size_w,
    REAL* coefficients)
{
    Spline4D spline_4d(data_size_x, data_size_y, data_size_z, data_size_w);
    spline_4d.initialize(data);
    spline_4d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_4d_portable(int argc, void* argv[])
{
    return calculate_coefficients_4d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        (REAL*)argv[5]);
}

int convert_csaps_coefficients_4d(
    REAL* csaps_coefficients,
    std::size_t n_spline_intervals_x,
    std::size_t n_spline_intervals_y,
    std::size_t n_spline_intervals_z,
    std::size_t n_spline_intervals_w,
    REAL* grid_spacing_array,
    REAL* reordered_coefficients)
{

    Spline4D spline_4d(n_spline_intervals_x + 1, n_spline_intervals_y + 1, n_spline_intervals_z + 1, n_spline_intervals_w + 1);
    spline_4d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals_x, n_spline_intervals_y, n_spline_intervals_z, n_spline_intervals_w, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_4d_portable(int argc, void* argv[])
{
    return convert_csaps_coefficients_4d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        (REAL*)argv[5],
        (REAL*)argv[6]);
}

int interpolate_4d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t data_size_w,
    std::size_t new_size_x,
    std::size_t new_size_y,
    std::size_t new_size_z,
    std::size_t new_size_w,
    REAL* x_values,
    REAL* y_values,
    REAL* z_values,
    REAL* w_values,
    REAL* interpolated_data)
{
    Spline4D spline_4d(data_size_x, data_size_y, data_size_z, data_size_w);

    spline_4d.interpolate(
        interpolated_data,
        data,
        x_values, y_values, z_values, w_values,
        new_size_x, new_size_y, new_size_z, new_size_w);

    return 0;
}

int interpolate_4d_portable(int argc, void* argv[])
{
    return interpolate_4d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        *((std::size_t*)argv[5]),
        *((std::size_t*)argv[6]),
        *((std::size_t*)argv[7]),
        *((std::size_t*)argv[8]),
        (REAL*)argv[9],
        (REAL*)argv[10],
        (REAL*)argv[11],
        (REAL*)argv[12],
        (REAL*)argv[13]);
}

int calculate_values_4d(
    REAL* coefficients,
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const n_intervals_z,
    std::size_t const n_intervals_w,
    std::size_t const values_size_x,
    std::size_t const values_size_y,
    std::size_t const values_size_z,
    std::size_t const values_size_w,
    REAL* x_values,
    REAL* y_values,
    REAL* z_values,
    REAL* w_values,
    REAL* spline_values)
{
    Spline4D spline_4d(n_intervals_x, n_intervals_y, n_intervals_z, n_intervals_w, coefficients);
    spline_4d.calculate_values(
        spline_values,
        x_values,
        y_values,
        z_values,
        w_values,
        values_size_x,
        values_size_y,
        values_size_z,
        values_size_w);

    return 0;
}

int calculate_values_4d_portable(int argc, void* argv[])
{
    return calculate_values_4d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        *((std::size_t*)argv[5]),
        *((std::size_t*)argv[6]),
        *((std::size_t*)argv[7]),
        *((std::size_t*)argv[8]),
        (REAL*)argv[9],
        (REAL*)argv[10],
        (REAL*)argv[11],
        (REAL*)argv[12],
        (REAL*)argv[13]);
}