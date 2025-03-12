#include "spline_classes.h"

EquationSystem_5D equation_system_5d;

Spline5D::Spline5D(std::size_t const data_size_x, std::size_t const data_size_y, std::size_t const data_size_z, std::size_t const data_size_w, std::size_t const data_size_v) :
    data_size_x_(data_size_x),
    data_size_y_(data_size_y),
    data_size_z_(data_size_z),
    data_size_w_(data_size_w),
	data_size_v_(data_size_v),
    n_intervals_x_(data_size_x - 1),
    n_intervals_y_(data_size_y - 1),
    n_intervals_z_(data_size_z - 1),
    n_intervals_w_(data_size_w - 1),
	n_intervals_v_(data_size_v - 1),
    n_intervals_((data_size_x - 1) * (data_size_y - 1) * (data_size_z - 1) * (data_size_w - 1) * (data_size_v - 1)),
    coefficients_calculated_(false),
    data_initialized_(false)
{}

Spline5D::Spline5D(
    std::size_t const n_intervals_x,
    std::size_t const n_intervals_y,
    std::size_t const n_intervals_z,
    std::size_t const n_intervals_w,
    std::size_t const n_intervals_v,
    REAL * coefficients)
    :
    data_size_x_(n_intervals_x + 1),
    data_size_y_(n_intervals_y + 1),
    data_size_z_(n_intervals_z + 1),
    data_size_w_(n_intervals_w + 1),
    data_size_v_(n_intervals_v + 1),
    n_intervals_x_(n_intervals_x),
    n_intervals_y_(n_intervals_y),
    n_intervals_z_(n_intervals_z),
    n_intervals_w_(n_intervals_w),
    n_intervals_v_(n_intervals_v),
    n_intervals_(n_intervals_x * n_intervals_y * n_intervals_z * n_intervals_w * n_intervals_v),
    coefficients_calculated_(true),
    data_initialized_(false),
    coefficients_(coefficients)
{}

Spline5D::~Spline5D()
{}

void Spline5D::initialize(REAL const * data)
{
    data_ = data;
    data_initialized_ = true;
}

void Spline5D::calculate_coefficients(REAL * coefficients) {
    coefficients_ = coefficients;

    // Step 1: Compute 4D splines for each v slice
    std::vector<Spline4D> xyzw_splines(data_size_v_, Spline4D(data_size_x_, data_size_y_, data_size_z_, data_size_w_));

    std::size_t n_spline_intervals_xyzw = xyzw_splines[0].n_intervals_;
    std::vector<REAL> xyzw_coefficients(data_size_v_ * n_spline_intervals_xyzw * Spline4D::n_coefficients_per_point);

    for (std::size_t v_index = 0; v_index < data_size_v_; v_index++) {
        REAL const* xyzw_data = data_ + v_index * data_size_x_ * data_size_y_ * data_size_z_ * data_size_w_;

        xyzw_splines[v_index].initialize(xyzw_data);
        xyzw_splines[v_index].calculate_coefficients(
            xyzw_coefficients.data() + v_index * n_spline_intervals_xyzw * Spline4D::n_coefficients_per_point);
    }

    // Step 2: Create 1D splines along the v-axis
    std::size_t const interpolated_size_x = 3 * n_intervals_x_ + 1;
    std::size_t const interpolated_size_y = 3 * n_intervals_y_ + 1;
    std::size_t const interpolated_size_z = 3 * n_intervals_z_ + 1;
    std::size_t const interpolated_size_w = 3 * n_intervals_w_ + 1;
    std::size_t const interpolated_size_xyzw = interpolated_size_x * interpolated_size_y * interpolated_size_z * interpolated_size_w;

    std::vector<Spline1D> v_splines(interpolated_size_xyzw, Spline1D(data_size_v_));
    std::vector<REAL> v_coefficients(interpolated_size_xyzw * n_intervals_v_ * Spline1D::n_coefficients_per_point);
    std::vector<REAL> v_data(data_size_v_); 

    // Interpolate along v-dimension
	for (std::size_t x_index = 0; x_index < interpolated_size_x; x_index++) {
		REAL const x = static_cast<REAL>(x_index) / 3;

		for (std::size_t y_index = 0; y_index < interpolated_size_y; y_index++) {
			REAL const y = static_cast<REAL>(y_index) / 3;

			for (std::size_t z_index = 0; z_index < interpolated_size_z; z_index++) {
				REAL const z = static_cast<REAL>(z_index) / 3;

				for (std::size_t w_index = 0; w_index < interpolated_size_w; w_index++) {
					REAL const w = static_cast<REAL>(w_index) / 3;

					std::size_t xyzw_index = (x_index * interpolated_size_y * interpolated_size_z * interpolated_size_w) +
											 (y_index * interpolated_size_z * interpolated_size_w) +
											 (z_index * interpolated_size_w) + w_index;

					// Compute values along v-axis
					for (std::size_t v_index = 0; v_index < data_size_v_; v_index++) {
						v_data[v_index] = xyzw_splines[v_index].calculate_value(x, y, z, w);
					}

					// Compute v-spline coefficients (only once per xyzw index)
					REAL* vc = v_coefficients.data() + xyzw_index * n_intervals_v_ * Spline1D::n_coefficients_per_point;

					v_splines[xyzw_index].initialize(v_data.data());
					v_splines[xyzw_index].calculate_coefficients(vc);
				}
			}
		}
	}


    // Step 3: Solve for final 5D spline coefficients
    for (std::size_t i = 0; i < n_intervals_x_; i++) {
        for (std::size_t j = 0; j < n_intervals_y_; j++) {
            for (std::size_t k = 0; k < n_intervals_z_; k++) {
                for (std::size_t q = 0; q < n_intervals_w_; q++) {
                    for (std::size_t r = 0; r < n_intervals_v_; r++) {
                        REAL* current_coefficients = coefficients
                            + ((i * n_intervals_y_ * n_intervals_z_ * n_intervals_w_ * n_intervals_v_)
                            +  (j * n_intervals_z_ * n_intervals_w_ * n_intervals_v_)
                            +  (k * n_intervals_w_ * n_intervals_v_)
                            +  (q * n_intervals_v_) + r)
                            * n_coefficients_per_point;

                        equation_system_5d.set_vector(v_splines, interpolated_size_x, interpolated_size_y, interpolated_size_z, interpolated_size_w, i, j, k, q, r);
                        equation_system_5d.solve(current_coefficients);
                    }
                }
            }
        }
    }

    coefficients_calculated_ = true;
}

REAL Spline5D::calculate_value(REAL const x, REAL const y, REAL const z, REAL const w, REAL const v) 
{
    // Determine interval indices (floor operation to find the correct cell)
    int i = static_cast<int>(std::floor(x));
    int j = static_cast<int>(std::floor(y));
    int k = static_cast<int>(std::floor(z));
    int m = static_cast<int>(std::floor(w));
    int n = static_cast<int>(std::floor(v));

    // Clamp indices to valid range
    i = clamp(i, 0, static_cast<int>(n_intervals_x_ - 1));
    j = clamp(j, 0, static_cast<int>(n_intervals_y_ - 1));
    k = clamp(k, 0, static_cast<int>(n_intervals_z_ - 1));
    m = clamp(m, 0, static_cast<int>(n_intervals_w_ - 1));
    n = clamp(n, 0, static_cast<int>(n_intervals_v_ - 1));

    // Compute fractional positions within each interval
    REAL const x_diff = x - static_cast<REAL>(i);
    REAL const y_diff = y - static_cast<REAL>(j);
    REAL const z_diff = z - static_cast<REAL>(k);
    REAL const w_diff = w - static_cast<REAL>(m);
    REAL const v_diff = v - static_cast<REAL>(n);

    // Compute the spline value using a nested polynomial expansion
    REAL spline_value = 0;
    REAL power_factor_i = 1;
    for (std::size_t order_i = 0; order_i < Spline1D::n_coefficients_per_point; order_i++) {
        REAL power_factor_j = 1;
        for (std::size_t order_j = 0; order_j < Spline1D::n_coefficients_per_point; order_j++) {
            REAL power_factor_k = 1;
            for (std::size_t order_k = 0; order_k < Spline1D::n_coefficients_per_point; order_k++) {
                REAL power_factor_m = 1;
                for (std::size_t order_m = 0; order_m < Spline1D::n_coefficients_per_point; order_m++) {
                    REAL power_factor_n = 1;
                    for (std::size_t order_n = 0; order_n < Spline1D::n_coefficients_per_point; order_n++) {
                        
                        // Compute the correct coefficient index
                        std::size_t point_index = (i * n_intervals_y_ * n_intervals_z_ * n_intervals_w_ * n_intervals_v_) +
                                                  (j * n_intervals_z_ * n_intervals_w_ * n_intervals_v_) +
                                                  (k * n_intervals_w_ * n_intervals_v_) +
                                                  (m * n_intervals_v_) + n;

                        std::size_t coeff_index = (order_i * 256) + (order_j * 64) + (order_k * 16) + (order_m * 4) + order_n;

                        // Accumulate the weighted coefficient contributions
                        spline_value += coefficients_[point_index * n_coefficients_per_point + coeff_index] *
                                        power_factor_i * power_factor_j * power_factor_k * power_factor_m * power_factor_n;

                        power_factor_n *= v_diff;
                    }
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

void Spline5D::calculate_values(
    REAL * spline_values,
    REAL const * x_values,
    REAL const * y_values,
    REAL const * z_values,
    REAL const * w_values,
    REAL const * v_values,
    std::size_t const size_x,
    std::size_t const size_y,
    std::size_t const size_z,
    std::size_t const size_w,
    std::size_t const size_v)
{
    for (std::size_t v_index = 0; v_index < size_v; v_index++) {
        REAL const v = v_values[v_index];

        for (std::size_t w_index = 0; w_index < size_w; w_index++) {
            REAL const w = w_values[w_index];

            for (std::size_t z_index = 0; z_index < size_z; z_index++) {
                REAL const z = z_values[z_index];

                for (std::size_t y_index = 0; y_index < size_y; y_index++) {
                    REAL const y = y_values[y_index];

                    for (std::size_t x_index = 0; x_index < size_x; x_index++) {
                        REAL const x = x_values[x_index];

                        // Compute 1D flattened index
                        std::size_t point_index =
                            (v_index * size_x * size_y * size_z * size_w) +
                            (w_index * size_x * size_y * size_z) +
                            (z_index * size_x * size_y) +
                            (y_index * size_x) +
                            x_index;

                        // Evaluate spline function
                        spline_values[point_index] = calculate_value(x, y, z, w, v);
                    }
                }
            }
        }
    }
}

void Spline5D::interpolate(
    REAL * interpolated_data,
    REAL const * data,
    REAL const * x_values,
    REAL const * y_values,
    REAL const * z_values,
    REAL const * w_values,
    REAL const * v_values,
    std::size_t const size_x,
    std::size_t const size_y,
    std::size_t const size_z,
    std::size_t const size_w,
    std::size_t const size_v)
{
    // Step 1: Initialize the spline with the given data
    initialize(data);

    // Step 2: Compute the 5D spline coefficients
    std::vector<REAL> coefficients(n_intervals_ * n_coefficients_per_point);
    calculate_coefficients(coefficients.data());

    // Step 3: Evaluate the spline at the requested points
    calculate_values(
        interpolated_data,
        x_values,
        y_values,
        z_values,
        w_values,
        v_values,
        size_x,
        size_y,
        size_z,
        size_w,
        size_v);
}

void Spline5D::convert_csaps_coefficients(
    REAL * csaps_coefficients,
    std::size_t const n_spline_intervals_x,
    std::size_t const n_spline_intervals_y,
    std::size_t const n_spline_intervals_z,
    std::size_t const n_spline_intervals_w,
    std::size_t const n_spline_intervals_v,
    REAL * grid_spacing_array,
    REAL * reordered_coefficients)
{
    // Extract grid spacings
    REAL dx = grid_spacing_array[0];
    REAL dy = grid_spacing_array[1];
    REAL dz = grid_spacing_array[2];
    REAL dw = grid_spacing_array[3];
    REAL dv = grid_spacing_array[4];

    // Compute scale factors for each dimension
    REAL dx_scale_factors[4], dy_scale_factors[4], dz_scale_factors[4], dw_scale_factors[4], dv_scale_factors[4];

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

    dv_scale_factors[0] = 1.0;
    dv_scale_factors[1] = dv;
    dv_scale_factors[2] = dv * dv;
    dv_scale_factors[3] = dv * dv * dv;

    // Compute coefficients per point in 1D splines
    std::size_t ncpp1d = Spline1D::n_coefficients_per_point;
    std::size_t ncpp1d_2 = ncpp1d * ncpp1d;
    std::size_t ncpp1d_3 = ncpp1d * ncpp1d * ncpp1d;
    std::size_t ncpp1d_4 = ncpp1d * ncpp1d * ncpp1d * ncpp1d;

    // Compute array strides
    std::size_t nwv = n_spline_intervals_w * n_spline_intervals_v;
    std::size_t nzwv = n_spline_intervals_z * nwv;
    std::size_t nyzwv = n_spline_intervals_y * nzwv;
    std::size_t nxyzwv = n_spline_intervals_x * nyzwv;

    // Loop through all spline intervals
    for (std::size_t v_index = 0; v_index < n_spline_intervals_v; v_index++) {
        for (std::size_t w_index = 0; w_index < n_spline_intervals_w; w_index++) {
            for (std::size_t z_index = 0; z_index < n_spline_intervals_z; z_index++) {
                for (std::size_t y_index = 0; y_index < n_spline_intervals_y; y_index++) {
                    for (std::size_t x_index = 0; x_index < n_spline_intervals_x; x_index++) {
                        for (std::size_t order_v = 0; order_v < Spline1D::n_coefficients_per_point; order_v++) {
                            for (std::size_t order_w = 0; order_w < Spline1D::n_coefficients_per_point; order_w++) {
                                for (std::size_t order_z = 0; order_z < Spline1D::n_coefficients_per_point; order_z++) {
                                    for (std::size_t order_y = 0; order_y < Spline1D::n_coefficients_per_point; order_y++) {
                                        for (std::size_t order_x = 0; order_x < Spline1D::n_coefficients_per_point; order_x++) {

                                            // Compute standard index (our internal storage)
                                            std::size_t std_point_index =
                                                (x_index * nyzwv) + (y_index * nzwv) + (z_index * nwv) + (w_index * n_spline_intervals_v) + v_index;
                                            
                                            std::size_t std_coeff_index =
                                                (order_x * ncpp1d_4) + (order_y * ncpp1d_3) + (order_z * ncpp1d_2) + (order_w * ncpp1d) + order_v;
                                            
                                            std::size_t combined_std_coeff_index =
                                                std_point_index * n_coefficients_per_point + std_coeff_index;

                                            // Compute CSAPS index (external format)
                                            std::size_t combined_csaps_coeff_index =
                                                ((ncpp1d - 1) - order_x) * ncpp1d_4 * nxyzwv +
                                                ((ncpp1d - 1) - order_y) * ncpp1d_3 * nxyzwv +
                                                ((ncpp1d - 1) - order_z) * ncpp1d_2 * nxyzwv +
                                                ((ncpp1d - 1) - order_w) * ncpp1d * nxyzwv +
                                                ((ncpp1d - 1) - order_v) * nxyzwv +
                                                (x_index * nyzwv) +
                                                (y_index * nzwv) +
                                                (z_index * nwv) +
                                                (w_index * n_spline_intervals_v) +
                                                v_index;

                                            // Apply scaling and copy coefficients
                                            reordered_coefficients[combined_std_coeff_index] =
                                                csaps_coefficients[combined_csaps_coeff_index] *
                                                dx_scale_factors[order_x] *
                                                dy_scale_factors[order_y] *
                                                dz_scale_factors[order_z] *
                                                dw_scale_factors[order_w] *
                                                dv_scale_factors[order_v];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
