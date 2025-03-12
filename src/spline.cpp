#include "spline.h"
#include "spline_classes.h"

int calculate_coefficients_1d(
    REAL * data,
    std::size_t data_size_x,
    REAL * coefficients)
{
    Spline1D spline_1d(data_size_x);
    spline_1d.initialize(data);
    spline_1d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_2d(
    REAL * data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    REAL * coefficients)
{
    Spline2D spline_2d(data_size_x, data_size_y);
    spline_2d.initialize(data);
    spline_2d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_3d(
    REAL * data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    REAL * coefficients)
{
    Spline3D spline_3d(data_size_x, data_size_y, data_size_z);
    spline_3d.initialize(data);
    spline_3d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_4d(
    REAL * data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t data_size_w,
    REAL * coefficients)
{
    Spline4D spline_4d(data_size_x, data_size_y, data_size_z, data_size_w);
    spline_4d.initialize(data);
    spline_4d.calculate_coefficients(coefficients);

    return 0;
}

int calculate_coefficients_5d(
    REAL * data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t data_size_w,
    std::size_t data_size_v,
    REAL * coefficients)
{
    Spline5D spline_5d(data_size_x, data_size_y, data_size_z, data_size_w, data_size_v);
    spline_5d.initialize(data);
    spline_5d.calculate_coefficients(coefficients);

    return 0;
}

int convert_csaps_coefficients_1d(
    REAL * csaps_coefficients,
    std::size_t n_spline_intervals,
    REAL * grid_spacing_array,
    REAL * reordered_coefficients)
{

    Spline1D spline_1d(n_spline_intervals + 1);
    spline_1d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_2d(
    REAL * csaps_coefficients,
    std::size_t n_spline_intervals_x,
    std::size_t n_spline_intervals_y,
    REAL * grid_spacing_array,
    REAL * reordered_coefficients)
{

    Spline2D spline_2d(n_spline_intervals_x + 1, n_spline_intervals_y + 1);
    spline_2d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals_x, n_spline_intervals_y, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_3d(
    REAL * csaps_coefficients,
    std::size_t n_spline_intervals_x,
    std::size_t n_spline_intervals_y,
    std::size_t n_spline_intervals_z,
    REAL * grid_spacing_array,
    REAL * reordered_coefficients)
{

    Spline3D spline_3d(n_spline_intervals_x + 1, n_spline_intervals_y + 1, n_spline_intervals_z + 1);
    spline_3d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals_x, n_spline_intervals_y, n_spline_intervals_z, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_4d(
    REAL * csaps_coefficients,
    std::size_t n_spline_intervals_x,
    std::size_t n_spline_intervals_y,
    std::size_t n_spline_intervals_z,
    std::size_t n_spline_intervals_w,
    REAL * grid_spacing_array,
    REAL * reordered_coefficients)
{

    Spline4D spline_4d(n_spline_intervals_x + 1, n_spline_intervals_y + 1, n_spline_intervals_z + 1, n_spline_intervals_w + 1);
    spline_4d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals_x, n_spline_intervals_y, n_spline_intervals_z, n_spline_intervals_w, grid_spacing_array, reordered_coefficients);

    return 0;
}

int convert_csaps_coefficients_5d(
    REAL * csaps_coefficients,
    std::size_t n_spline_intervals_x,
    std::size_t n_spline_intervals_y,
    std::size_t n_spline_intervals_z,
    std::size_t n_spline_intervals_w,
    std::size_t n_spline_intervals_v,
    REAL * grid_spacing_array,
    REAL * reordered_coefficients)
{

    Spline5D spline_5d(n_spline_intervals_x + 1, n_spline_intervals_y + 1, n_spline_intervals_z + 1, n_spline_intervals_w + 1, n_spline_intervals_v + 1);
    spline_5d.convert_csaps_coefficients(csaps_coefficients, n_spline_intervals_x, n_spline_intervals_y, n_spline_intervals_z, n_spline_intervals_w, n_spline_intervals_v, grid_spacing_array, reordered_coefficients);

    return 0;
}

int interpolate_1d(
    REAL * data,
    std::size_t data_size_x,
    std::size_t new_size_x,
    REAL * x_values,
    REAL * interpolated_data)
{
    Spline1D spline_1d(data_size_x);

    spline_1d.interpolate(
        interpolated_data,
        data,
        x_values,
        new_size_x);

    return 0;
}

int interpolate_2d(
    REAL * data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t new_size_x,
    std::size_t new_size_y,
    REAL * x_values,
    REAL * y_values,
    REAL * interpolated_data)
{
    Spline2D spline_2d(data_size_x, data_size_y);

    spline_2d.interpolate(
        interpolated_data,
        data,
        x_values, y_values,
        new_size_x, new_size_y);

    return 0;
}

int interpolate_3d(
    REAL * data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t new_size_x,
    std::size_t new_size_y,
    std::size_t new_size_z,
    REAL * x_values,
    REAL * y_values,
    REAL * z_values,
    REAL * interpolated_data)
{
    Spline3D spline_3d(data_size_x, data_size_y, data_size_z);

    spline_3d.interpolate(
        interpolated_data,
        data,
        x_values, y_values, z_values,
        new_size_x, new_size_y, new_size_z);

    return 0;
}

int interpolate_4d(
    REAL * data,
    std::size_t data_size_x,
    std::size_t data_size_y,
    std::size_t data_size_z,
    std::size_t data_size_w,
    std::size_t new_size_x,
    std::size_t new_size_y,
    std::size_t new_size_z,
    std::size_t new_size_w,
    REAL * x_values,
    REAL * y_values,
    REAL * z_values,
    REAL * w_values,
    REAL * interpolated_data)
{
    Spline4D spline_4d(data_size_x, data_size_y, data_size_z, data_size_w);

    spline_4d.interpolate(
        interpolated_data,
        data,
        x_values, y_values, z_values, w_values,
        new_size_x, new_size_y, new_size_z, new_size_w);

    return 0;
}

int interpolate_5d(
    REAL * data,
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
    REAL * x_values,
    REAL * y_values,
    REAL * z_values,
    REAL * w_values,
    REAL * v_values,
    REAL * interpolated_data)
{
    Spline5D spline_5d(data_size_x, data_size_y, data_size_z, data_size_w, data_size_v);

    spline_5d.interpolate(
        interpolated_data,
        data,
        x_values, y_values, z_values, w_values, v_values,
        new_size_x, new_size_y, new_size_z, new_size_w, new_size_v);

    return 0;
}

int calculate_values_1d(
    REAL * coefficients,
    std::size_t const n_intervals_x,
    std::size_t const values_size_x,
    REAL * x_values,
    REAL* spline_values)
{
    Spline1D spline_1d(n_intervals_x, coefficients);
    spline_1d.calculate_values(spline_values, x_values, values_size_x);

    return 0;
}

int calculate_values_2d(
    REAL * coefficients,
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const values_size_x,
    std::size_t const values_size_y,
    REAL * x_values,
    REAL * y_values,
    REAL * spline_values)
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

int calculate_values_3d(
    REAL * coefficients,
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const n_intervals_z,
    std::size_t const values_size_x,
    std::size_t const values_size_y,
    std::size_t const values_size_z,
    REAL * x_values,
    REAL * y_values,
    REAL * z_values,
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

int calculate_values_4d(
    REAL * coefficients,
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const n_intervals_z,
    std::size_t const n_intervals_w,
    std::size_t const values_size_x,
    std::size_t const values_size_y,
    std::size_t const values_size_z,
    std::size_t const values_size_w,
    REAL * x_values,
    REAL * y_values,
    REAL * z_values,
    REAL * w_values,
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

int calculate_values_5d(
    REAL * coefficients,
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
    REAL * x_values,
    REAL * y_values,
    REAL * z_values,
    REAL * w_values,
    REAL * v_values,
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

int calculate_coefficients_1d_portable(int argc, void *argv[])
{
    return calculate_coefficients_1d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        (REAL *)argv[2]);
}

int calculate_coefficients_2d_portable(int argc, void *argv[])
{
    return calculate_coefficients_2d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        (REAL *)argv[3]);
}

int calculate_coefficients_3d_portable(int argc, void *argv[])
{
    return calculate_coefficients_3d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        (REAL *)argv[4]);
}

int calculate_coefficients_4d_portable(int argc, void *argv[])
{
    return calculate_coefficients_4d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        (REAL *)argv[5]);
}

int calculate_coefficients_5d_portable(int argc, void *argv[])
{
    return calculate_coefficients_5d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        *((std::size_t *)argv[5]),
        (REAL *)argv[6]);
}

int convert_csaps_coefficients_1d_portable(int argc, void *argv[])
{
    return convert_csaps_coefficients_1d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        (REAL *)argv[2],
        (REAL *)argv[3]);
}

int convert_csaps_coefficients_2d_portable(int argc, void *argv[])
{
    return convert_csaps_coefficients_2d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        (REAL *)argv[3],
        (REAL *)argv[4]);
}

int convert_csaps_coefficients_3d_portable(int argc, void *argv[])
{
    return convert_csaps_coefficients_3d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        (REAL *)argv[4],
        (REAL *)argv[5]);
}

int convert_csaps_coefficients_4d_portable(int argc, void *argv[])
{
    return convert_csaps_coefficients_4d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        (REAL *)argv[5],
        (REAL *)argv[6]);
}

int convert_csaps_coefficients_5d_portable(int argc, void *argv[])
{
    return convert_csaps_coefficients_5d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        *((std::size_t *)argv[5]),
        (REAL *)argv[6],
        (REAL *)argv[7]);
}

int interpolate_1d_portable(int argc, void *argv[])
{
    return interpolate_1d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        (REAL *)argv[3],
        (REAL *)argv[4]);
}

int interpolate_2d_portable(int argc, void *argv[])
{
    return interpolate_2d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        (REAL *)argv[5],
        (REAL *)argv[6],
        (REAL *)argv[7]);
}

int interpolate_3d_portable(int argc, void *argv[])
{
    return interpolate_3d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        *((std::size_t *)argv[5]),
        *((std::size_t *)argv[6]),
        (REAL *)argv[7],
        (REAL *)argv[8],
        (REAL *)argv[9],
        (REAL *)argv[10]);
}

int interpolate_4d_portable(int argc, void *argv[])
{
    return interpolate_4d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        *((std::size_t *)argv[5]),
        *((std::size_t *)argv[6]),
        *((std::size_t *)argv[7]),
        *((std::size_t *)argv[8]),
        (REAL *)argv[9],
        (REAL *)argv[10],
        (REAL *)argv[11],
        (REAL *)argv[12],
        (REAL *)argv[13]);
}

int interpolate_5d_portable(int argc, void *argv[])
{
    return interpolate_5d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        *((std::size_t *)argv[5]),
        *((std::size_t *)argv[6]),
        *((std::size_t *)argv[7]),
        *((std::size_t *)argv[8]),
        *((std::size_t *)argv[9]),
        *((std::size_t *)argv[10]),
        (REAL *)argv[11],
        (REAL *)argv[12],
        (REAL *)argv[13],
        (REAL *)argv[14],
        (REAL *)argv[15],
        (REAL *)argv[16]);
}

int calculate_values_1d_portable(int argc, void *argv[])
{
    return calculate_values_1d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        (REAL *)argv[3],
        (REAL *)argv[4]);
}

int calculate_values_2d_portable(int argc, void *argv[])
{
    return calculate_values_2d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        (REAL *)argv[5],
        (REAL *)argv[6],
        (REAL *)argv[7]);
}

int calculate_values_3d_portable(int argc, void *argv[])
{
    return calculate_values_3d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        *((std::size_t *)argv[5]),
        *((std::size_t *)argv[6]),
        (REAL *)argv[7],
        (REAL *)argv[8],
        (REAL *)argv[9],
        (REAL *)argv[10]);
}

int calculate_values_4d_portable(int argc, void *argv[])
{
    return calculate_values_4d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        *((std::size_t *)argv[5]),
        *((std::size_t *)argv[6]),
        *((std::size_t *)argv[7]),
        *((std::size_t *)argv[8]),
        (REAL *)argv[9],
        (REAL *)argv[10],
        (REAL *)argv[11],
        (REAL *)argv[12],
        (REAL *)argv[13]);
}

int calculate_values_5d_portable(int argc, void *argv[])
{
    return calculate_values_5d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        *((std::size_t *)argv[3]),
        *((std::size_t *)argv[4]),
        *((std::size_t *)argv[5]),
        *((std::size_t *)argv[6]),
        *((std::size_t *)argv[7]),
        *((std::size_t *)argv[8]),
        *((std::size_t *)argv[9]),
        *((std::size_t *)argv[10]),
        (REAL *)argv[11],
        (REAL *)argv[12],
        (REAL *)argv[13],
        (REAL *)argv[14],
        (REAL *)argv[15],
        (REAL *)argv[16]);
}


int bspline_calculate_coefficients_1d(
    REAL * data,
    int data_size_x,
    REAL * coefficients)
{
    BSpline1D bspline_1d(data_size_x);
    bspline_1d.initialize(data);
    bspline_1d.calculate_coefficients(coefficients);

    return 0;
}


int bspline_calculate_coefficients_1d_portable(int argc, void *argv[])
{
    return bspline_calculate_coefficients_1d(
        (REAL *)argv[0],
        *((int *)argv[1]),
        (REAL *)argv[2]);
}


int bspline_calculate_values_1d(
    REAL * coefficients,
    std::size_t const n_intervals_x,
    std::size_t const values_size_x,
    REAL * x_values,
    REAL* bspline_values)
{
    BSpline1D bspline_1d(n_intervals_x, coefficients);
    bspline_1d.calculate_values(bspline_values, x_values, values_size_x);

    return 0;
}


int bspline_calculate_values_1d_portable(int argc, void *argv[])
{
    return bspline_calculate_values_1d(
        (REAL *)argv[0],
        *((std::size_t *)argv[1]),
        *((std::size_t *)argv[2]),
        (REAL *)argv[3],
        (REAL *)argv[4]);
}