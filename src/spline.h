#ifndef SPLINE_H_INCLUDED
#define SPLINE_H_INCLUDED

#ifdef __linux__
#define VISIBLE __attribute__((visibility("default")))
#endif

#ifdef _WIN32
#define VISIBLE
#endif

#include <cstddef>
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    VISIBLE int calculate_coefficients_1d(
        REAL * data,
        std::size_t data_size_x,
        REAL * coefficients);

    VISIBLE int calculate_coefficients_2d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t data_size_y,
        REAL * coefficients);

    VISIBLE int calculate_coefficients_3d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t data_size_y,
        std::size_t data_size_z,
        REAL * coefficients);

    VISIBLE int calculate_coefficients_4d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t data_size_y,
        std::size_t data_size_z,
        std::size_t data_size_w,
        REAL * coefficients);

    VISIBLE int interpolate_1d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t new_size_x,
        REAL * x_values,
        REAL * interpolated_data);

    VISIBLE int interpolate_2d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t data_size_y,
        std::size_t new_size_x,
        std::size_t new_size_y,
        REAL * x_values,
        REAL * y_values,
        REAL * interpolated_data);

    VISIBLE int interpolate_3d(
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
        REAL * interpolated_data);

    VISIBLE int interpolate_4d(
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
        REAL * interpolated_data);

    VISIBLE int calculate_values_1d(
        REAL * coefficients,
        std::size_t const n_intervals_x,
        std::size_t const values_size_x,
        REAL * x_values,
        REAL * spline_values);

    VISIBLE int calculate_values_2d(
        REAL * coefficients,
        std::size_t const n_intervals_x,
        std::size_t const n_intervals_y,
        std::size_t const values_size_x,
        std::size_t const values_size_y,
        REAL * x_values,
        REAL * y_values,
        REAL * spline_values);

    VISIBLE int calculate_values_3d(
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
        REAL * spline_values);

    VISIBLE int calculate_values_4d(
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
        REAL * spline_values);

    VISIBLE int convert_csaps_coefficients_1d(
        REAL * csaps_coefficients,
        std::size_t n_spline_intervals,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

    VISIBLE int convert_csaps_coefficients_2d(
        REAL * csaps_coefficients,
        std::size_t n_spline_intervals_x,
        std::size_t n_spline_intervals_y,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

    VISIBLE int convert_csaps_coefficients_3d(
        REAL * csaps_coefficients,
        std::size_t n_spline_intervals_x,
        std::size_t n_spline_intervals_y,
        std::size_t n_spline_intervals_z,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

    VISIBLE int convert_csaps_coefficients_4d(
        REAL * csaps_coefficients,
        std::size_t n_spline_intervals_x,
        std::size_t n_spline_intervals_y,
        std::size_t n_spline_intervals_z,
        std::size_t n_spline_intervals_w,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

    VISIBLE int calculate_coefficients_1d_portable(int argc, void *argv[]);

    VISIBLE int calculate_coefficients_2d_portable(int argc, void *argv[]);

    VISIBLE int calculate_coefficients_3d_portable(int argc, void *argv[]);

    VISIBLE int calculate_coefficients_4d_portable(int argc, void *argv[]);

    VISIBLE int interpolate_1d_portable(int argc, void *argv[]);

    VISIBLE int interpolate_2d_portable(int argc, void *argv[]);
    
    VISIBLE int interpolate_3d_portable(int argc, void *argv[]);

    VISIBLE int interpolate_4d_portable(int argc, void *argv[]);

    VISIBLE int calculate_values_1d_portable(int argc, void *argv[]);

    VISIBLE int calculate_values_2d_portable(int argc, void *argv[]);

    VISIBLE int calculate_values_3d_portable(int argc, void *argv[]);

    VISIBLE int calculate_values_4d_portable(int argc, void *argv[]);

    VISIBLE int convert_csaps_coefficients_1d_portable(int argc, void *argv[]);

    VISIBLE int convert_csaps_coefficients_2d_portable(int argc, void *argv[]);

    VISIBLE int convert_csaps_coefficients_3d_portable(int argc, void *argv[]);

    VISIBLE int convert_csaps_coefficients_4d_portable(int argc, void *argv[]);

#ifdef __cplusplus
}
#endif

#endif // SPLINE_H_INCLUDED