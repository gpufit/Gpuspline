#include "spline_classes.h"

EquationSystem_1D equation_system_1d;

Spline1D::Spline1D(std::size_t const data_size) :
    data_size_(data_size),
    n_intervals_(data_size - 1),
    coefficients_calculated_(false),
    data_initialized_(false)
{}

Spline1D::Spline1D(std::size_t const n_intervals, REAL * coefficients) :
    data_size_(n_intervals + 1),
    n_intervals_(n_intervals),
    coefficients_calculated_(true),
    data_initialized_(false),
    coefficients_(coefficients)
{}

Spline1D::~Spline1D()
{}

void Spline1D::initialize(REAL const * data)
{
    data_ = data;
    data_initialized_ = true;
}

void Spline1D::set_coefficients(REAL * const coefficients)
{
    coefficients_ = coefficients;
}

// this follows http://mathworld.wolfram.com/CubicSpline.html which sets the second derivatives at the endpoints to zero
void Spline1D::calculate_coefficients(
    REAL * coefficients)
{
    coefficients_ = coefficients;

    equation_system_1d.resize(data_size_);

    equation_system_1d.set_vector(data_);

    equation_system_1d.solve();

    for (std::size_t i = 0; i < equation_system_1d.N_ - 1; i++)
    {
        // coefficients c[0]+c[1]*t+c[2]*t^2+c[3]*t^3
        REAL * c = coefficients + Spline1D::n_coefficients_per_point * i;
        // y = data
        REAL const * y = data_;
        REAL const * D = equation_system_1d.solution_.data();

        c[0] = y[i];
        c[1] = D[i];
        c[2] = 3 * (y[i + 1] - y[i]) - 2 * D[i] - D[i + 1];
        c[3] = 2 * (y[i] - y[i + 1]) + D[i] + D[i + 1];
    }

    coefficients_calculated_ = true;
}

REAL Spline1D::calculate_value(REAL const x)
{
    // determine interval index
    int i = static_cast<std::size_t>(std::floor(x));

    // adjust i to its bounds TODO this behavior should be documented (constant extrapolation)
    i = i >= 0 ? i : 0;
    i = i < static_cast<int>(n_intervals_) ? i : static_cast<int>(n_intervals_) - 1;

    // determine fraction in interval
    REAL const x_diff = x - static_cast<REAL>(i);

    // calculate spline value in interval
    REAL power_factor = 1;
    REAL spline_value = 0;
    for (std::size_t order = 0; order < n_coefficients_per_point; order++)
    {
        spline_value
            += coefficients_[i * n_coefficients_per_point + order]
            * power_factor;
        power_factor *= x_diff;
    }

    return spline_value;
}

// calls calculate_value() iteratively
void Spline1D::calculate_values(
    REAL * spline_values,
    REAL const * x_values,
    std::size_t const size_x)
{
    for (std::size_t x_index = 0; x_index < size_x; x_index++)
    {

        if (x_index == (size_x - 1)){

            int a = 0;

        }

        if (x_index == (size_x - 2)){

            int b = 0;

        }

        spline_values[x_index] = calculate_value(x_values[x_index]);
    }
}

void Spline1D::interpolate(
    REAL * interpolated_data,
    REAL const * data,
    REAL const * x_values,
    std::size_t const size_x)
{
    initialize(data);

    std::vector<REAL> coefficients(n_intervals_ * n_coefficients_per_point);

    calculate_coefficients(coefficients.data());

    calculate_values(
        interpolated_data,
        x_values,
        size_x);
}

// change the ordering of a coefficients array
void Spline1D::convert_csaps_coefficients(
    REAL * csaps_coefficients,
    std::size_t const n_spline_intervals,
    REAL * grid_spacing_array,
    REAL * reordered_coefficients)
{

    REAL dx = grid_spacing_array[0];

    REAL dx_scale_factors[4];
    
    dx_scale_factors[0] = 1.0;
    dx_scale_factors[1] = dx;
    dx_scale_factors[2] = dx * dx;
    dx_scale_factors[3] = dx * dx * dx;

    for (std::size_t x_index = 0; x_index < n_spline_intervals; x_index++)
    {
        for (std::size_t order_index = 0; order_index < n_coefficients_per_point; order_index++){

            std::size_t csaps_index = ((n_coefficients_per_point - 1) - order_index)*n_spline_intervals + x_index;
            std::size_t std_index = x_index * n_coefficients_per_point + order_index;

            reordered_coefficients[std_index] = csaps_coefficients[csaps_index] * dx_scale_factors[order_index];

        }
    }
}