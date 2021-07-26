#ifndef SPLINE_H_INCLUDED
#define SPLINE_H_INCLUDED

#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    int calculate_coefficients_1d(
        REAL * data,
        std::size_t data_size_x,
        REAL * coefficients);

    int calculate_coefficients_2d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t data_size_y,
        REAL * coefficients);

    int calculate_coefficients_3d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t data_size_y,
        std::size_t data_size_z,
        REAL * coefficients);

    int interpolate_1d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t new_size_x,
        REAL * x_values,
        REAL * interpolated_data);

    int interpolate_2d(
        REAL * data,
        std::size_t data_size_x,
        std::size_t data_size_y,
        std::size_t new_size_x,
        std::size_t new_size_y,
        REAL * x_values,
        REAL * y_values,
        REAL * interpolated_data);

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
        REAL * interpolated_data);

    int calculate_values_1d(
        REAL * coefficients,
        std::size_t const n_intervals_x,
        std::size_t const values_size_x,
        REAL * x_values,
        REAL * spline_values);

    int calculate_values_2d(
        REAL * coefficients,
        std::size_t const n_intervals_x,
        std::size_t const n_intervals_y,
        std::size_t const values_size_x,
        std::size_t const values_size_y,
        REAL * x_values,
        REAL * y_values,
        REAL * spline_values);

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
        REAL * spline_values);

    int calculate_coefficients_1d_portable(int argc, void *argv[]);

    int calculate_coefficients_2d_portable(int argc, void *argv[]);

    int calculate_coefficients_3d_portable(int argc, void *argv[]);

    int interpolate_1d_portable(int argc, void *argv[]);

    int interpolate_2d_portable(int argc, void *argv[]);
    
    int interpolate_3d_portable(int argc, void *argv[]);

    int calculate_values_1d_portable(int argc, void *argv[]);

    int calculate_values_2d_portable(int argc, void *argv[]);

    int calculate_values_3d_portable(int argc, void *argv[]);

#ifdef __cplusplus
}
#endif

#endif // !SPLINE_H_INCLUDED