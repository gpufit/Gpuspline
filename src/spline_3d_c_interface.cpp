#include "spline.h"
#include "spline_classes.h"

int calculate_coefficients_3d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    REAL* coefficients)
{
    Spline3D spline_3d(data_size_x, data_size_y, data_size_z);
    spline_3d.initialize(data);
    spline_3d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_3d_portable(int argc, void* argv[])
{
    return calculate_coefficients_3d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        (REAL*)argv[4]);
}

int convert_csaps_coefficients_3d(
    REAL* csaps_coefficients,
    std::size_t n_spline_intervals_x,
    std::size_t n_spline_intervals_y,
    std::size_t n_spline_intervals_z,
    REAL* grid_spacing_array,
    REAL* reordered_coefficients)
{

    Spline3D spline_3d(n_spline_intervals_x + 1, n_spline_intervals_y + 1, n_spline_intervals_z + 1);
    spline_3d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals_x, n_spline_intervals_y, n_spline_intervals_z, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_3d_portable(int argc, void* argv[])
{
    return convert_csaps_coefficients_3d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        (REAL*)argv[4],
        (REAL*)argv[5]);
}

int interpolate_3d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t new_size_x,
    std::size_t new_size_y,
    std::size_t new_size_z,
    REAL* x_values,
    REAL* y_values,
    REAL* z_values,
    REAL* interpolated_data)
{
    Spline3D spline_3d(data_size_x, data_size_y, data_size_z);

    spline_3d.interpolate(
        interpolated_data,
        data,
        x_values, y_values, z_values,
        new_size_x, new_size_y, new_size_z);

    return 0;
}

int interpolate_3d_portable(int argc, void* argv[])
{
    return interpolate_3d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        *((std::size_t*)argv[5]),
        *((std::size_t*)argv[6]),
        (REAL*)argv[7],
        (REAL*)argv[8],
        (REAL*)argv[9],
        (REAL*)argv[10]);
}

int calculate_values_3d(
    REAL* coefficients,
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const n_intervals_z,
    std::size_t const values_size_x,
    std::size_t const values_size_y,
    std::size_t const values_size_z,
    REAL* x_values,
    REAL* y_values,
    REAL* z_values,
    REAL* spline_values)
{
    Spline3D spline_3d(n_intervals_x, n_intervals_y, n_intervals_z, coefficients);
    spline_3d.calculate_values(
        spline_values,
        x_values,
        y_values,
        z_values,
        values_size_x,
        values_size_y,
        values_size_z);

    return 0;
}

int calculate_values_3d_portable(int argc, void* argv[])
{
    return calculate_values_3d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        *((std::size_t*)argv[5]),
        *((std::size_t*)argv[6]),
        (REAL*)argv[7],
        (REAL*)argv[8],
        (REAL*)argv[9],
        (REAL*)argv[10]);
}