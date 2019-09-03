#include "spline_classes.h"

EquationSystem_3D equation_system_3d;

Spline3D::Spline3D(std::size_t const data_size_x, std::size_t const data_size_y, std::size_t const data_size_z) :
    data_size_x_(data_size_x),
    data_size_y_(data_size_y),
    data_size_z_(data_size_z),
    n_intervals_x_(data_size_x - 1),
    n_intervals_y_(data_size_y - 1),
    n_intervals_z_(data_size_z - 1),
    n_intervals_(n_intervals_x_ * n_intervals_y_ * n_intervals_z_),
    coefficients_calculated_(false),
    data_initialized_(false)
{}

Spline3D::Spline3D(
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const n_intervals_z,
    REAL * coefficients)
    :
    data_size_x_(n_intervals_x + 1),
    data_size_y_(n_intervals_y + 1),
    data_size_z_(n_intervals_z + 1),
    n_intervals_x_(n_intervals_x),
    n_intervals_y_(n_intervals_y),
    n_intervals_z_(n_intervals_z),
    n_intervals_(n_intervals_x_ * n_intervals_y_ * n_intervals_z_),
    coefficients_calculated_(true),
    data_initialized_(false),
    coefficients_(coefficients)
{}

Spline3D::~Spline3D()
{}

void Spline3D::initialize(REAL const * data)
{
    data_ = data;
    data_initialized_ = true;
}

// follows https://github.com/ZhuangLab/storm-analysis/blob/master/storm_analysis/spliner/spline3D.py
void Spline3D::calculate_coefficients(REAL * coefficients)
{
    coefficients_ = coefficients;

    // create 2D splines along z-axis for each xy plane = xy_splines
    std::vector<Spline2D> xy_splines(data_size_z_, Spline2D(data_size_x_, data_size_y_));

    std::size_t const n_spline_intervals_xy = xy_splines[0].n_intervals_;
    std::vector<REAL> xy_coefficients(data_size_z_ * n_spline_intervals_xy * Spline2D::n_coefficients_per_point);

    for (std::size_t z_index = 0; z_index < data_size_z_; z_index++)
    {
        REAL const * xy_data = data_ + z_index * data_size_x_ * data_size_y_;

        xy_splines[z_index].initialize(xy_data);
        xy_splines[z_index].calculate_coefficients(
            xy_coefficients.data() + z_index * n_spline_intervals_xy * Spline2D::n_coefficients_per_point);
    }

    // use xy_splines to create splines along the z-axis with sub-integer spacing withing the xy-plane (factor of 3 times more to get 64 values per 3D interval in total) = z_splines
    std::size_t const interpolated_size_x = 3 * n_intervals_x_ + 1;
    std::size_t const interpolated_size_y = 3 * n_intervals_y_ + 1;
    std::size_t const interpolated_size_xy = interpolated_size_x * interpolated_size_y;
    std::vector<Spline1D> z_splines(interpolated_size_xy, Spline1D(data_size_z_));

    // order: z coefficients, y-axis, x-axis
    std::vector<REAL> z_coefficients(interpolated_size_xy * n_intervals_z_ * Spline1D::n_coefficients_per_point);

    std::vector<REAL> z_data(data_size_z_);

    // for each xy spline
    for (std::size_t x_index = 0; x_index < interpolated_size_x; x_index++)
    {
        // 0, 1/3, 2/3, 1, ...
        REAL const x = static_cast<REAL>(x_index) / 3;

        for (std::size_t y_index = 0; y_index < interpolated_size_y; y_index++)
        {
            // 0, 1/3, 2/3, 1, ...
            REAL const y = static_cast<REAL>(y_index) / 3;

            // for each point along the z-axis, calculate spline values from xy_splines and use as new data values
            for (std::size_t z_index = 0; z_index < data_size_z_; z_index++)
            {
                z_data[z_index] = xy_splines[z_index].calculate_value(x, y);
            }

            // order in z_splines is dy first, then dx
            std::size_t xy_index = x_index * interpolated_size_y + y_index;
            // initialize with spline interpolated values and calculate coefficients
            REAL * zc = z_coefficients.data() + xy_index * n_intervals_z_ * Spline1D::n_coefficients_per_point;

            z_splines[xy_index].initialize(z_data.data());
            z_splines[xy_index].calculate_coefficients(zc);
        }
    }

    // compute spline coefficients using the z-axis splines to generate 64 values per 3D interval and then
    // solving for the coefficients
    for (std::size_t i = 0; i < n_intervals_x_; i++)
    {
        for (std::size_t j = 0; j < n_intervals_y_; j++)
        {
            for (std::size_t k = 0; k < n_intervals_z_; k++)
            {
                REAL * current_coefficients
                    = coefficients
                    + (i * n_intervals_y_ * n_intervals_z_ + j * n_intervals_z_ + k)
                    * n_coefficients_per_point;
                
                equation_system_3d.set_vector(z_splines, interpolated_size_x, interpolated_size_y, i, j, k);

                equation_system_3d.solve(current_coefficients);
            }
        }
    }

    coefficients_calculated_ = true;
}

REAL Spline3D::calculate_value(REAL const x, REAL const y, REAL const z)
{
    // determine interval index
    int i = static_cast<std::size_t>(std::floor(x));
    int j = static_cast<std::size_t>(std::floor(y));
    int k = static_cast<std::size_t>(std::floor(z));

    // adjust i, j and k to their bounds
    i = i >= 0 ? i : 0;
    i = i < static_cast<int>(n_intervals_x_) ? i : static_cast<int>(n_intervals_x_) - 1;
    j = j >= 0 ? j : 0;
    j = j < static_cast<int>(n_intervals_y_) ? j : static_cast<int>(n_intervals_y_) - 1;
    k = k >= 0 ? k : 0;
    k = k < static_cast<int>(n_intervals_z_) ? k : static_cast<int>(n_intervals_z_) - 1;

    // get fraction in interval
    REAL const x_diff = x - static_cast<REAL>(i);
    REAL const y_diff = y - static_cast<REAL>(j);
    REAL const z_diff = z - static_cast<REAL>(k);

    REAL spline_value = 0;
    REAL power_factor_i = 1;
    for (std::size_t order_i = 0; order_i < Spline1D::n_coefficients_per_point; order_i++)
    {
        REAL power_factor_j = 1;
        for (std::size_t order_j = 0; order_j < Spline1D::n_coefficients_per_point; order_j++)
        {
            REAL power_factor_k = 1;
            for (std::size_t order_k = 0; order_k < Spline1D::n_coefficients_per_point; order_k++)
            {
                std::size_t const point_index
                    = i * n_intervals_y_ * n_intervals_z_ + j * n_intervals_z_ + k;
                
                std::size_t const coeff_index = order_i * 16 + order_j * 4 + order_k;

                spline_value
                    += coefficients_[point_index * n_coefficients_per_point + coeff_index]
                    * power_factor_i
                    * power_factor_j
                    * power_factor_k;

                power_factor_k *= z_diff;
            }
            power_factor_j *= y_diff;
        }
        power_factor_i *= x_diff;
    }

    return spline_value;
}

void Spline3D::calculate_values(
    REAL * spline_values,
    REAL const * x_values,
    REAL const * y_values,
    REAL const * z_values,
    std::size_t const size_x,
    std::size_t const size_y,
    std::size_t const size_z)
{
    for (std::size_t z_index = 0; z_index < size_z; z_index++)
    {
        for (std::size_t y_index = 0; y_index < size_y; y_index++)
        {
            for (std::size_t x_index = 0; x_index < size_x; x_index++)
            {
                std::size_t point_index
                    = z_index * size_x * size_y + y_index * size_x + x_index;
                
                spline_values[point_index]
                    = calculate_value(x_values[x_index], y_values[y_index], z_values[z_index]);
            }
        }
    }
}

void Spline3D::interpolate(
    REAL * interpolated_data,
    REAL const * data,
    REAL const * x_values,
    REAL const * y_values,
    REAL const * z_values,
    std::size_t const size_x,
    std::size_t const size_y,
    std::size_t const size_z)
{
    initialize(data);

    std::vector<REAL> coefficients(n_intervals_ * n_coefficients_per_point);

    calculate_coefficients(coefficients.data());

    calculate_values(
        interpolated_data,
        x_values,
        y_values,
        z_values,
        size_x,
        size_y,
        size_z);
}
