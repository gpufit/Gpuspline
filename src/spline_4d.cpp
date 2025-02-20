#include "spline_classes.h"

EquationSystem_4D equation_system_4d;

Spline4D::Spline4D(std::size_t const data_size_x, std::size_t const data_size_y, std::size_t const data_size_z, std::size_t const data_size_w) :
    data_size_x_(data_size_x),
    data_size_y_(data_size_y),
    data_size_z_(data_size_z),
    data_size_w_(data_size_w),
    n_intervals_x_(data_size_x - 1),
    n_intervals_y_(data_size_y - 1),
    n_intervals_z_(data_size_z - 1),
    n_intervals_w_(data_size_w - 1),
    n_intervals_((data_size_x - 1) * (data_size_y - 1) * (data_size_z - 1) * (data_size_w - 1)),
    coefficients_calculated_(false),
    data_initialized_(false)
{}

Spline4D::Spline4D(
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const n_intervals_z,
    std::size_t const n_intervals_w,
    REAL * coefficients)
    :
    data_size_x_(n_intervals_x + 1),
    data_size_y_(n_intervals_y + 1),
    data_size_z_(n_intervals_z + 1),
    data_size_w_(n_intervals_w + 1),
    n_intervals_x_(n_intervals_x),
    n_intervals_y_(n_intervals_y),
    n_intervals_z_(n_intervals_z),
    n_intervals_w_(n_intervals_w),
    n_intervals_(n_intervals_x * n_intervals_y * n_intervals_z * n_intervals_w),
    coefficients_calculated_(true),
    data_initialized_(false),
    coefficients_(coefficients)
{}

Spline4D::~Spline4D()
{}

void Spline4D::initialize(REAL const * data)
{
    data_ = data;
    data_initialized_ = true;
}

//
//void Spline4D::calculate_coefficients(REAL * coefficients)
//{
//    coefficients_ = coefficients;
//
//    // Step 1: Create 3D splines for each w slice
//    std::vector<Spline3D> xyz_splines(data_size_w_, Spline3D(data_size_x_, data_size_y_, data_size_z_));
//
//    std::size_t const n_spline_intervals_xyz = xyz_splines[0].n_intervals_;
//    std::vector<REAL> xyz_coefficients(data_size_w_ * n_spline_intervals_xyz * Spline3D::n_coefficients_per_point);
//
//    for (std::size_t w_index = 0; w_index < data_size_w_; w_index++)
//    {
//        REAL const * xyz_data = data_ + w_index * data_size_x_ * data_size_y_ * data_size_z_;
//
//        xyz_splines[w_index].initialize(xyz_data);
//        xyz_splines[w_index].calculate_coefficients(
//            xyz_coefficients.data() + w_index * n_spline_intervals_xyz * Spline3D::n_coefficients_per_point);
//    }
//
//    // Step 2: Create 1D splines along the w-axis
//    std::size_t const interpolated_size_x = 3 * n_intervals_x_ + 1;
//    std::size_t const interpolated_size_y = 3 * n_intervals_y_ + 1;
//    std::size_t const interpolated_size_z = 3 * n_intervals_z_ + 1;
//    std::size_t const interpolated_size_xyz = interpolated_size_x * interpolated_size_y * interpolated_size_z;
//    std::vector<Spline1D> w_splines(interpolated_size_xyz, Spline1D(data_size_w_));
//
//    // Order: w coefficients, z-axis, y-axis, x-axis
//    std::vector<REAL> w_coefficients(interpolated_size_xyz * n_intervals_w_ * Spline1D::n_coefficients_per_point);
//    std::vector<REAL> w_data(data_size_w_);
//
//    // Step 3: Compute interpolated values along w-axis
//    for (std::size_t x_index = 0; x_index < interpolated_size_x; x_index++)
//    {
//        REAL const x = static_cast<REAL>(x_index) / 3;
//
//        for (std::size_t y_index = 0; y_index < interpolated_size_y; y_index++)
//        {
//            REAL const y = static_cast<REAL>(y_index) / 3;
//
//            for (std::size_t z_index = 0; z_index < interpolated_size_z; z_index++)
//            {
//                REAL const z = static_cast<REAL>(z_index) / 3;
//
//                // Compute new values along the w-axis using xyz_splines
//                for (std::size_t w_index = 0; w_index < data_size_w_; w_index++)
//                {
//                    w_data[w_index] = xyz_splines[w_index].calculate_value(x, y, z);
//                }
//
//                // Corrected calculation for xyz_index
//                std::size_t xyz_index = (x_index * interpolated_size_y * interpolated_size_z)
//                    + (y_index * interpolated_size_z)
//                    + z_index;
//
//                // Initialize with spline interpolated values and calculate coefficients
//                REAL * wc = w_coefficients.data() + xyz_index * n_intervals_w_ * Spline1D::n_coefficients_per_point;
//
//                w_splines[xyz_index].initialize(w_data.data());
//                w_splines[xyz_index].calculate_coefficients(wc);
//            }
//        }
//    }
//
//    // Step 4: Solve for 4D spline coefficients
//    for (std::size_t i = 0; i < n_intervals_x_; i++)
//    {
//        for (std::size_t j = 0; j < n_intervals_y_; j++)
//        {
//            for (std::size_t k = 0; k < n_intervals_z_; k++)
//            {
//                for (std::size_t m = 0; m < n_intervals_w_; m++)  // Renamed from "l" to "m"
//                {
//                    REAL * current_coefficients
//                        = coefficients
//                        + ((i * n_intervals_y_ * n_intervals_z_ * n_intervals_w_)
//                        + (j * n_intervals_z_ * n_intervals_w_)
//                        + (k * n_intervals_w_) + m)  // Updated indexing
//                        * n_coefficients_per_point;
//
//                    equation_system_4d.set_vector(w_splines, interpolated_size_x, interpolated_size_y, interpolated_size_z, i, j, k, m);
//
//                    equation_system_4d.solve(current_coefficients);
//                }
//            }
//        }
//    }
//
//    coefficients_calculated_ = true;
//}
//

void Spline4D::calculate_coefficients(REAL * coefficients)
{
    coefficients_ = coefficients;

    // create 3D splines along w-axis for each xyz volume = xyz_splines
    std::vector<Spline3D> xyz_splines(data_size_w_, Spline3D(data_size_x_, data_size_y_, data_size_z_));

    std::size_t const n_spline_intervals_xyz = xyz_splines[0].n_intervals_;
    std::vector<REAL> xyz_coefficients(data_size_w_ * n_spline_intervals_xyz * Spline3D::n_coefficients_per_point);

    for (std::size_t w_index = 0; w_index < data_size_w_; w_index++)
    {
        REAL const * xyz_data = data_ + w_index * data_size_x_ * data_size_y_ * data_size_z_;

        xyz_splines[w_index].initialize(xyz_data);
        xyz_splines[w_index].calculate_coefficients(
            xyz_coefficients.data() + w_index * n_spline_intervals_xyz * Spline3D::n_coefficients_per_point);
    }

    // use xyz_splines to create splines along the w-axis with sub-integer spacing within the xyz-volume (factor of 3 times more to get 256 values per 4D interval in total) = w_splines
    std::size_t const interpolated_size_x = 3 * n_intervals_x_ + 1;
    std::size_t const interpolated_size_y = 3 * n_intervals_y_ + 1;
    std::size_t const interpolated_size_z = 3 * n_intervals_z_ + 1;
    std::size_t const interpolated_size_xyz = interpolated_size_x * interpolated_size_y * interpolated_size_z;
    std::vector<Spline1D> w_splines(interpolated_size_xyz, Spline1D(data_size_w_));

    // order: w coefficients, z-axis, y-axis, x-axis
    std::vector<REAL> t_coefficients(interpolated_size_xyz * n_intervals_w_ * Spline1D::n_coefficients_per_point);

    std::vector<REAL> w_data(data_size_w_);

    // for each xyz spline
    for (std::size_t x_index = 0; x_index < interpolated_size_x; x_index++)
    {
        // 0, 1/3, 2/3, 1, ...
        REAL const x = static_cast<REAL>(x_index) / 3;

        for (std::size_t y_index = 0; y_index < interpolated_size_y; y_index++)
        {
            // 0, 1/3, 2/3, 1, ...
            REAL const y = static_cast<REAL>(y_index) / 3;

            for (std::size_t z_index = 0; z_index < interpolated_size_z; z_index++)
            {
                // 0, 1/3, 2/3, 1, ...
                REAL const z = static_cast<REAL>(z_index) / 3;

                // for each point along the z-axis, calculate spline values from xy_splines and use as new data values
                for (std::size_t w_index = 0; w_index < data_size_w_; w_index++)
                {
                    w_data[w_index] = xyz_splines[w_index].calculate_value(x, y, z);
                }

                // order in w_splines is dz, dy, dx
                std::size_t xyz_index = x_index * interpolated_size_y * interpolated_size_z + y_index * interpolated_size_z + z_index;
                // initialize with spline interpolated values and calculate coefficients
                REAL * tc = t_coefficients.data() + xyz_index * n_intervals_w_ * Spline1D::n_coefficients_per_point;

                w_splines[xyz_index].initialize(w_data.data());
                w_splines[xyz_index].calculate_coefficients(tc);

            }
        }
    }

    // compute spline coefficients using the w-axis splines to generate 256 values per 4D interval and then
    // solving for the coefficients
    for (std::size_t i = 0; i < n_intervals_x_; i++)
    {
        for (std::size_t j = 0; j < n_intervals_y_; j++)
        {
            for (std::size_t k = 0; k < n_intervals_z_; k++)
            {
                for (std::size_t m = 0; m < n_intervals_w_; m++)
                {

                    REAL * current_coefficients
                        = coefficients
                        + (i * n_intervals_y_ * n_intervals_z_ * n_intervals_w_ + j * n_intervals_z_ * n_intervals_w_ + k * n_intervals_w_ + m)
                        * n_coefficients_per_point;
                
                    equation_system_4d.set_vector(w_splines, interpolated_size_x, interpolated_size_y, interpolated_size_z, i, j, k, m);

                    equation_system_4d.solve(current_coefficients);
                }
            }
        }
    }

    coefficients_calculated_ = true;
}


REAL Spline4D::calculate_value(REAL const x, REAL const y, REAL const z, REAL const w)
{
    // determine interval index
    int i = static_cast<std::size_t>(std::floor(x));
    int j = static_cast<std::size_t>(std::floor(y));
    int k = static_cast<std::size_t>(std::floor(z));
    int m = static_cast<std::size_t>(std::floor(w));

    // adjust i, j, k, m to their bounds
    i = i >= 0 ? i : 0;
    i = i < static_cast<int>(n_intervals_x_) ? i : static_cast<int>(n_intervals_x_) - 1;
    j = j >= 0 ? j : 0;
    j = j < static_cast<int>(n_intervals_y_) ? j : static_cast<int>(n_intervals_y_) - 1;
    k = k >= 0 ? k : 0;
    k = k < static_cast<int>(n_intervals_z_) ? k : static_cast<int>(n_intervals_z_) - 1;
    m = m >= 0 ? m : 0;
    m = m < static_cast<int>(n_intervals_w_) ? m : static_cast<int>(n_intervals_w_) - 1;

    // get fraction in interval
    REAL const x_diff = x - static_cast<REAL>(i);
    REAL const y_diff = y - static_cast<REAL>(j);
    REAL const z_diff = z - static_cast<REAL>(k);
    REAL const w_diff = w - static_cast<REAL>(m);

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
                REAL power_factor_m = 1;
                for (std::size_t order_m = 0; order_m < Spline1D::n_coefficients_per_point; order_m++)
                {

                    std::size_t const point_index
                        = i * n_intervals_y_ * n_intervals_z_ * n_intervals_w_ + j * n_intervals_z_ * n_intervals_w_ + k * n_intervals_w_ + m;
                
                    std::size_t const coeff_index = order_i * 64 + order_j * 16 + order_k * 4 + order_m;

                    spline_value
                        += coefficients_[point_index * n_coefficients_per_point + coeff_index]
                        * power_factor_i
                        * power_factor_j
                        * power_factor_k
                        * power_factor_m;

                    power_factor_m *= w_diff;
                }
                power_factor_k *= z_diff;
            }
            power_factor_j *= y_diff;
        }
        power_factor_i *= x_diff;
    }

    return spline_value;
}

void Spline4D::calculate_values(
    REAL * spline_values,
    REAL const * x_values,
    REAL const * y_values,
    REAL const * z_values,
    REAL const * w_values,
    std::size_t const size_x,
    std::size_t const size_y,
    std::size_t const size_z,
    std::size_t const size_w)
{
    for (std::size_t w_index = 0; w_index < size_w; w_index++)
    {
        for (std::size_t z_index = 0; z_index < size_z; z_index++)
        {
            for (std::size_t y_index = 0; y_index < size_y; y_index++)
            {
                for (std::size_t x_index = 0; x_index < size_x; x_index++)
                {
                    std::size_t point_index
                        = w_index * size_x * size_y * size_z + z_index * size_x * size_y + y_index * size_x + x_index;
                
                    spline_values[point_index]
                        = calculate_value(x_values[x_index], y_values[y_index], z_values[z_index], w_values[w_index]);
                }
            }
        }
    }
}

void Spline4D::interpolate(
    REAL * interpolated_data,
    REAL const * data,
    REAL const * x_values,
    REAL const * y_values,
    REAL const * z_values,
    REAL const * w_values,
    std::size_t const size_x,
    std::size_t const size_y,
    std::size_t const size_z,
    std::size_t const size_w)
{

    initialize(data);

    std::vector<REAL> coefficients(n_intervals_ * n_coefficients_per_point);

    calculate_coefficients(coefficients.data());

    calculate_values(
        interpolated_data,
        x_values,
        y_values,
        z_values,
        w_values,
        size_x,
        size_y,
        size_z,
        size_w);

}


// change the ordering of a coefficients array
void Spline4D::convert_csaps_coefficients(
    REAL * csaps_coefficients,
    std::size_t const n_spline_intervals_x,
    std::size_t const n_spline_intervals_y,
    std::size_t const n_spline_intervals_z,
    std::size_t const n_spline_intervals_w,
    REAL * grid_spacing_array,
    REAL * reordered_coefficients)
{

    REAL dx = grid_spacing_array[0];
    REAL dy = grid_spacing_array[1];
    REAL dz = grid_spacing_array[2];
    REAL dw = grid_spacing_array[3];

    REAL dx_scale_factors[4], dy_scale_factors[4], dz_scale_factors[4], dw_scale_factors[4];

    dx_scale_factors[0] = 1.0;
    dx_scale_factors[1] = dx;
    dx_scale_factors[2] = dx * dx;
    dx_scale_factors[3] = dx * dx * dx;

    dy_scale_factors[0] = 1.0;
    dy_scale_factors[1] = dy;
    dy_scale_factors[2] = dy * dy;
    dy_scale_factors[3] = dy * dy * dy;

    dz_scale_factors[0] = 1.0;
    dz_scale_factors[1] = dz;
    dz_scale_factors[2] = dz * dz;
    dz_scale_factors[3] = dz * dz * dz;

    dw_scale_factors[0] = 1.0;
    dw_scale_factors[1] = dw;
    dw_scale_factors[2] = dw * dw;
    dw_scale_factors[3] = dw * dw * dw;

    std::size_t ncpp1d = Spline1D::n_coefficients_per_point;

    std::size_t ncpp1d_2 = ncpp1d * ncpp1d;
    std::size_t ncpp1d_3 = ncpp1d * ncpp1d * ncpp1d;

    std::size_t nzw = n_spline_intervals_z * n_spline_intervals_w;
    std::size_t nyzw = n_spline_intervals_y * n_spline_intervals_z * n_spline_intervals_w;
    std::size_t nxyzw = n_spline_intervals_x * n_spline_intervals_y * n_spline_intervals_z * n_spline_intervals_w;

    for (std::size_t w_index = 0; w_index < n_spline_intervals_w; w_index++)
    {
        for (std::size_t z_index = 0; z_index < n_spline_intervals_z; z_index++)
        {
            for (std::size_t y_index = 0; y_index < n_spline_intervals_y; y_index++)
            {
                for (std::size_t x_index = 0; x_index < n_spline_intervals_x; x_index++)
                {
                    for (std::size_t order_w_index = 0; order_w_index < Spline1D::n_coefficients_per_point; order_w_index++)
                    {
                        for (std::size_t order_z_index = 0; order_z_index < Spline1D::n_coefficients_per_point; order_z_index++)
                        {
                            for (std::size_t order_y_index = 0; order_y_index < Spline1D::n_coefficients_per_point; order_y_index++)
                            {
                                for (std::size_t order_x_index = 0; order_x_index < Spline1D::n_coefficients_per_point; order_x_index++)
                                {

                                    std::size_t const std_point_index = (x_index * nyzw) + (y_index * nzw) + z_index * n_spline_intervals_w + w_index;
                                    std::size_t const std_coeff_index = order_x_index * ncpp1d_3 + order_y_index * ncpp1d_2 + order_z_index * ncpp1d + order_w_index;
                                    std::size_t const combined_std_coeff_index = std_point_index * n_coefficients_per_point + std_coeff_index;


                                    std::size_t const combined_csaps_coeff_index = ((ncpp1d - 1) - order_x_index) * ncpp1d_3 * nxyzw +
                                                                                   ((ncpp1d - 1) - order_y_index) * ncpp1d_2 * nxyzw +
                                                                                   ((ncpp1d - 1) - order_z_index) * ncpp1d * nxyzw +
                                                                                   ((ncpp1d - 1) - order_w_index) * nxyzw +
                                                                                   (x_index) * nyzw +
                                                                                   (y_index) * nzw +
                                                                                   (z_index) * n_spline_intervals_w + 
                                                                                   (w_index);

                                    reordered_coefficients[combined_std_coeff_index] = csaps_coefficients[combined_csaps_coeff_index] *
                                                                                           dx_scale_factors[order_x_index] *
                                                                                           dy_scale_factors[order_y_index] *
                                                                                           dz_scale_factors[order_z_index] *
                                                                                           dw_scale_factors[order_w_index];

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
