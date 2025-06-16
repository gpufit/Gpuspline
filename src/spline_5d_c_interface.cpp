#include "spline.h"
#include "spline_classes.h"

int calculate_coefficients_5d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t data_size_w,
    std::size_t data_size_v,
    REAL* coefficients)
{
    Spline5D spline_5d(data_size_x, data_size_y, data_size_z, data_size_w, data_size_v);
    spline_5d.initialize(data);
    spline_5d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_5d_portable(int argc, void* argv[])
{
    return calculate_coefficients_5d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        *((std::size_t*)argv[5]),
        (REAL*)argv[6]);
}

int convert_csaps_coefficients_5d(
    REAL* csaps_coefficients,
    std::size_t n_spline_intervals_x,
    std::size_t n_spline_intervals_y,
    std::size_t n_spline_intervals_z,
    std::size_t n_spline_intervals_w,
    std::size_t n_spline_intervals_v,
    REAL* grid_spacing_array,
    REAL* reordered_coefficients)
{

    Spline5D spline_5d(n_spline_intervals_x + 1, n_spline_intervals_y + 1, n_spline_intervals_z + 1, n_spline_intervals_w + 1, n_spline_intervals_v + 1);
    spline_5d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals_x, n_spline_intervals_y, n_spline_intervals_z, n_spline_intervals_w, n_spline_intervals_v, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_5d_portable(int argc, void* argv[])
{
    return convert_csaps_coefficients_5d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        *((std::size_t*)argv[5]),
        (REAL*)argv[6],
        (REAL*)argv[7]);
}

int interpolate_5d(
    REAL* data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t data_size_w,
    std::size_t data_size_v,
    std::size_t new_size_x,
    std::size_t new_size_y,
    std::size_t new_size_z,
    std::size_t new_size_w,
    std::size_t new_size_v,
    REAL* x_values,
    REAL* y_values,
    REAL* z_values,
    REAL* w_values,
    REAL* v_values,
    REAL* interpolated_data)
{
    Spline5D spline_5d(data_size_x, data_size_y, data_size_z, data_size_w, data_size_v);

    spline_5d.interpolate(
        interpolated_data,
        data,
        x_values, y_values, z_values, w_values, v_values,
        new_size_x, new_size_y, new_size_z, new_size_w, new_size_v);

    return 0;
}

int interpolate_5d_portable(int argc, void* argv[])
{
    return interpolate_5d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        *((std::size_t*)argv[5]),
        *((std::size_t*)argv[6]),
        *((std::size_t*)argv[7]),
        *((std::size_t*)argv[8]),
        *((std::size_t*)argv[9]),
        *((std::size_t*)argv[10]),
        (REAL*)argv[11],
        (REAL*)argv[12],
        (REAL*)argv[13],
        (REAL*)argv[14],
        (REAL*)argv[15],
        (REAL*)argv[16]);
}

int calculate_values_5d(
    REAL* coefficients,
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const n_intervals_z,
    std::size_t const n_intervals_w,
    std::size_t const n_intervals_v,
    std::size_t const values_size_x,
    std::size_t const values_size_y,
    std::size_t const values_size_z,
    std::size_t const values_size_w,
    std::size_t const values_size_v,
    REAL* x_values,
    REAL* y_values,
    REAL* z_values,
    REAL* w_values,
    REAL* v_values,
    REAL* spline_values)
{
    Spline5D spline_5d(n_intervals_x, n_intervals_y, n_intervals_z, n_intervals_w, n_intervals_v, coefficients);
    spline_5d.calculate_values(
        spline_values,
        x_values,
        y_values,
        z_values,
        w_values,
        v_values,
        values_size_x,
        values_size_y,
        values_size_z,
        values_size_w,
        values_size_v);

    return 0;
}

int calculate_values_5d_portable(int argc, void* argv[])
{
    return calculate_values_5d(
        (REAL*)argv[0],
        *((std::size_t*)argv[1]),
        *((std::size_t*)argv[2]),
        *((std::size_t*)argv[3]),
        *((std::size_t*)argv[4]),
        *((std::size_t*)argv[5]),
        *((std::size_t*)argv[6]),
        *((std::size_t*)argv[7]),
        *((std::size_t*)argv[8]),
        *((std::size_t*)argv[9]),
        *((std::size_t*)argv[10]),
        (REAL*)argv[11],
        (REAL*)argv[12],
        (REAL*)argv[13],
        (REAL*)argv[14],
        (REAL*)argv[15],
        (REAL*)argv[16]);
}