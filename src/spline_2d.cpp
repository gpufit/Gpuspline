#include "spline_classes.h"

EquationSystem_2D equation_system_2d;

Spline2D::Spline2D(std::size_t const data_size_x, std::size_t const data_size_y) :
    data_size_x_(data_size_x),
    data_size_y_(data_size_y),
    n_intervals_x_(data_size_x - 1),
    n_intervals_y_(data_size_y - 1),
    n_intervals_((data_size_x - 1) * (data_size_y - 1)),
    coefficients_calculated_(false),
    data_initialized_(false)
{}

Spline2D::Spline2D(
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    REAL * coefficients)
    :
    data_size_x_(n_intervals_x + 1),
    data_size_y_(n_intervals_y + 1),
    n_intervals_x_(n_intervals_x),
    n_intervals_y_(n_intervals_y),
    n_intervals_(n_intervals_x_ * n_intervals_y_),
    coefficients_calculated_(true),
    data_initialized_(false),
    coefficients_(coefficients)
{}

Spline2D::~Spline2D()
{}

void Spline2D::initialize(REAL const * data)
{
    data_ = data;
    data_initialized_ = true;
}

// follows https://github.com/ZhuangLab/storm-analysis/blob/master/storm_analysis/spliner/spline2D.py
void Spline2D::calculate_coefficients(REAL * coefficients)
{
    coefficients_ = coefficients;

    // create splines along x axis for each y value = y_splines
    std::vector<Spline1D> y_splines(data_size_y_, Spline1D(data_size_x_));
    std::vector<REAL> y_coefficients(data_size_y_ * n_intervals_x_ * Spline1D::n_coefficients_per_point);

    for (std::size_t y_index = 0; y_index < data_size_y_; y_index++)
    {
        REAL const * y_data = data_ + y_index * data_size_x_;
        REAL * yc = y_coefficients.data() + y_index * Spline1D::n_coefficients_per_point * n_intervals_x_;
        
        y_splines[y_index].initialize(y_data);
        y_splines[y_index].calculate_coefficients(yc);
    }

    // x_splines

    // use y_splines to creates splines along the y-axis with sub-integer spacing along x-axis (factor of 3 times more
    // to get 16 values per 2D interval in total) = x_splines
    std::vector<Spline1D> x_splines(3 * n_intervals_x_ + 1, Spline1D(data_size_y_));
    // splines are stored in x_splines along x-axis first
    std::vector<REAL> x_coefficients(x_splines.size() * n_intervals_y_ * Spline1D::n_coefficients_per_point);
    std::vector<REAL> x_data(data_size_y_);

    for (std::size_t x_index = 0; x_index < x_splines.size(); x_index++)
    {
        REAL x = static_cast<REAL>(x_index) / 3;

        // for each point along the y-axis, calculate spline value from y_splines and use as new data value
        for (std::size_t y_index = 0; y_index < data_size_y_; y_index++)
        {
            x_data[y_index] = y_splines[y_index].calculate_value(x);
        }

        // calculate coefficient for 1D splines
        REAL * xc = x_coefficients.data() + x_index * Spline1D::n_coefficients_per_point * n_intervals_y_;

        x_splines[x_index].initialize(x_data.data());
        x_splines[x_index].calculate_coefficients(xc);
    }

    // Compute spline coefficients using the x axis splines to generate 16 values per grid cell
    // and then solving for the coefficients.
    for (std::size_t i = 0; i < n_intervals_x_; i++)
    {
        for (std::size_t j = 0; j < n_intervals_y_; j++)
        {
            REAL * current_coefficients
                = coefficients
                + (i * n_intervals_y_ + j)
                * Spline2D::n_coefficients_per_point;
            
            equation_system_2d.set_vector(x_splines, i, j);

            equation_system_2d.solve(current_coefficients);
        }
    }

    coefficients_calculated_ = true;
}

REAL Spline2D::calculate_value(REAL const x, REAL const y)
{
    // determine interval index
    int i = static_cast<std::size_t>(std::floor(x));
    int j = static_cast<std::size_t>(std::floor(y));

    // adjust i, j and k to their bounds
    i = i >= 0 ? i : 0;
    i = i < static_cast<int>(n_intervals_x_) ? i : static_cast<int>(n_intervals_x_) - 1;
    j = j >= 0 ? j : 0;
    j = j < static_cast<int>(n_intervals_y_) ? j : static_cast<int>(n_intervals_y_) - 1;

    // get fraction in interval
    REAL const x_diff = x - static_cast<REAL>(i);
    REAL const y_diff = y - static_cast<REAL>(j);

    REAL spline_value = 0;
    REAL power_factor_i = 1;
    for (std::size_t order_i = 0; order_i < Spline1D::n_coefficients_per_point; order_i++)
    {
        REAL power_factor_j = 1;
        for (std::size_t order_j = 0; order_j < Spline1D::n_coefficients_per_point; order_j++)
        {
            std::size_t const point_index = i * n_intervals_y_ + j;
            std::size_t const coeff_index = order_i * Spline1D::n_coefficients_per_point + order_j;

            spline_value
                += coefficients_[point_index * n_coefficients_per_point + coeff_index]
                * power_factor_i
                * power_factor_j;

            power_factor_j *= y_diff;
        }
        power_factor_i *= x_diff;
    }

    return spline_value;
}

void Spline2D::calculate_values(
    REAL * spline_values,
    REAL const * x_values,
    REAL const * y_values,
    std::size_t const size_x,
    std::size_t const size_y)
{
    for (std::size_t y_index = 0; y_index < size_y; y_index++)
    {
        for (std::size_t x_index = 0; x_index < size_x; x_index++)
        {
            std::size_t point_index = y_index * size_x + x_index;
            spline_values[point_index] = calculate_value(x_values[x_index], y_values[y_index]);
        }
    }
}

void Spline2D::interpolate(
    REAL * interpolated_data,
    REAL const * data,
    REAL const * x_values,
    REAL const * y_values,
    std::size_t const size_x,
    std::size_t const size_y)
{
    initialize(data);

    std::vector<REAL> coefficients(n_intervals_ * n_coefficients_per_point);

    calculate_coefficients(coefficients.data());

    calculate_values(
        interpolated_data,
        x_values,
        y_values,
        size_x,
        size_y);
}