#include "spline_classes.h"



// Constructor (initializes an empty coefficient array)
BSpline5D::BSpline5D(int Nx, int Ny, int Nz, int Nw, int Nv)
    : Nx_(Nx), Ny_(Ny), Nz_(Nz), Nw_(Nw), Nv_(Nv),
    coefficients_(nullptr), coefficients_calculated_(false), data_initialized_(false) {
    int total_size = 4 * Nx_ * Ny_ * Nz_ * Nw_ * Nv_; // 4 coefficients per dimension
    coefficients_ = new REAL[total_size];
    std::memset(coefficients_, 0, total_size * sizeof(REAL));
}

// Constructor (accepts precomputed coefficients)
BSpline5D::BSpline5D(int Nx, int Ny, int Nz, int Nw, int Nv, REAL* coefficients)
    : Nx_(Nx), Ny_(Ny), Nz_(Nz), Nw_(Nw), Nv_(Nv),
    coefficients_(coefficients), coefficients_calculated_(true), data_initialized_(true) {}

// Destructor (frees allocated memory if needed)
BSpline5D::~BSpline5D() {
    if (coefficients_ != nullptr) {
        delete[] coefficients_;
    }
}

// Initialize coefficients from raw data
void BSpline5D::initialize(const REAL* data) {
    if (!data) {
        data_initialized_ = false;
        return;
    }

    // Process data to initialize coefficients
    // (Assuming some preprocessing or copying is required before calling calculate_coefficients)
    calculate_coefficients(coefficients_, data);

    data_initialized_ = true;
}




// Compute B-spline coefficients from input data
void BSpline5D::calculate_coefficients(REAL* coefficients, const REAL* data) {
    if (!data_initialized_) {
        return;
    }

    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            for (int z = 0; z < Nz_; z++) {
                for (int w = 0; w < Nw_; w++) {
                    for (int v = 0; v < Nv_; v++) {
                        // Compute the linear index explicitly with clear memory ordering
                        int index = x + Nx_ * (y + Ny_ * (z + Nz_ * (w + Nw_ * v)));
                        int coeff_index = index * 4; // Each grid point has 4 coefficients

                        // Compute coefficients based on data points
                        coefficients[coeff_index + 0] = data[index];  // Example: Direct assignment
                        coefficients[coeff_index + 1] = (x > 0) ? (data[index] - data[index - 1]) : 0;
                        coefficients[coeff_index + 2] = (y > 0) ? (data[index] - data[index - Nx_]) : 0;
                        coefficients[coeff_index + 3] = (z > 0) ? (data[index] - data[index - Nx_ * Ny_]) : 0;
                    }
                }
            }
        }
    }
    coefficients_calculated_ = true;
}



// Evaluate the B-spline function
REAL BSpline5D::evaluate(REAL x, REAL y, REAL z, REAL w, REAL v) const {
    if (!coefficients_calculated_) {
        return 0.0;
    }

    // Compute basis functions for each dimension
    REAL basis_x[4], basis_y[4], basis_z[4], basis_w[4], basis_v[4];
    compute_basis_functions(basis_x, x);
    compute_basis_functions(basis_y, y);
    compute_basis_functions(basis_z, z);
    compute_basis_functions(basis_w, w);
    compute_basis_functions(basis_v, v);

    // Determine the starting index for the given (x,y,z,w,v)
    int i = static_cast<int>(x);
    int j = static_cast<int>(y);
    int k = static_cast<int>(z);
    int m = static_cast<int>(w);
    int n = static_cast<int>(v);

    i = std::max(0, std::min(i, Nx_ - 2));
    j = std::max(0, std::min(j, Ny_ - 2));
    k = std::max(0, std::min(k, Nz_ - 2));
    m = std::max(0, std::min(m, Nw_ - 2));
    n = std::max(0, std::min(n, Nv_ - 2));

    // Compute the function value
    REAL value = 0.0;
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
                for (int d = 0; d < 4; ++d) {
                    for (int e = 0; e < 4; ++e) {
                        int index = (n + e) * (Nw_ * Nz_ * Ny_ * Nx_) +
                            (m + d) * (Nz_ * Ny_ * Nx_) +
                            (k + c) * (Ny_ * Nx_) +
                            (j + b) * Nx_ +
                            (i + a);
                        index *= 4; // Multiply by 4 to get the correct coefficient offset
                        REAL coeff = coefficients_[index];
                        value += coeff * basis_x[a] * basis_y[b] * basis_z[c] * basis_w[d] * basis_v[e];
                    }
                }
            }
        }
    }

    return value;
}

// Compute cubic B-spline basis functions for a given position
void BSpline5D::compute_basis_functions(REAL* basis_values, REAL x) const {
    REAL t = x - std::floor(x);

    basis_values[0] = (1.0 / 6.0) * (1 - t) * (1 - t) * (1 - t);
    basis_values[1] = (1.0 / 6.0) * (3 * t * t * t - 6 * t * t + 4);
    basis_values[2] = (1.0 / 6.0) * (-3 * t * t * t + 3 * t * t + 3 * t + 1);
    basis_values[3] = (1.0 / 6.0) * (t * t * t);
}


// Compute cubic B-spline derivative basis functions for a given position
void BSpline5D::compute_derivative_basis_functions(REAL* basis_derivatives, REAL x) const {
    REAL t = x - std::floor(x);

    basis_derivatives[0] = (-0.5) * (1 - t) * (1 - t);
    basis_derivatives[1] = (1.5 * t * t - 2 * t);
    basis_derivatives[2] = (-1.5 * t * t + t + 0.5);
    basis_derivatives[3] = (0.5 * t * t);
}


// Evaluate the B-spline function and its first derivatives
void BSpline5D::evaluate_with_derivatives(REAL x, REAL y, REAL z, REAL w, REAL v,
    REAL& function_value,
    REAL& dx, REAL& dy, REAL& dz, REAL& dw, REAL& dv) const {
    if (!coefficients_calculated_) {
        function_value = dx = dy = dz = dw = dv = 0.0;
        return;
    }

    // Compute basis functions and their derivatives
    REAL basis_x[4], basis_y[4], basis_z[4], basis_w[4], basis_v[4];
    REAL deriv_x[4], deriv_y[4], deriv_z[4], deriv_w[4], deriv_v[4];
    compute_basis_functions(basis_x, x);
    compute_basis_functions(basis_y, y);
    compute_basis_functions(basis_z, z);
    compute_basis_functions(basis_w, w);
    compute_basis_functions(basis_v, v);
    compute_derivative_basis_functions(deriv_x, x);
    compute_derivative_basis_functions(deriv_y, y);
    compute_derivative_basis_functions(deriv_z, z);
    compute_derivative_basis_functions(deriv_w, w);
    compute_derivative_basis_functions(deriv_v, v);

    // Determine the starting index for the given (x,y,z,w,v)
    int i = static_cast<int>(x);
    int j = static_cast<int>(y);
    int k = static_cast<int>(z);
    int m = static_cast<int>(w);
    int n = static_cast<int>(v);

    i = std::max(0, std::min(i, Nx_ - 2));
    j = std::max(0, std::min(j, Ny_ - 2));
    k = std::max(0, std::min(k, Nz_ - 2));
    m = std::max(0, std::min(m, Nw_ - 2));
    n = std::max(0, std::min(n, Nv_ - 2));

    // Compute function value and derivatives
    function_value = dx = dy = dz = dw = dv = 0.0;
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
                for (int d = 0; d < 4; ++d) {
                    for (int e = 0; e < 4; ++e) {
                        int index = (i + a) + Nx_ * ((j + b) + Ny_ * ((k + c) + Nz_ * ((m + d) + Nw_ * (n + e))));
                        index *= 4;
                        REAL coeff = coefficients_[index];
                        function_value += coeff * basis_x[a] * basis_y[b] * basis_z[c] * basis_w[d] * basis_v[e];
                        dx += coeff * deriv_x[a] * basis_y[b] * basis_z[c] * basis_w[d] * basis_v[e];
                        dy += coeff * basis_x[a] * deriv_y[b] * basis_z[c] * basis_w[d] * basis_v[e];
                        dz += coeff * basis_x[a] * basis_y[b] * deriv_z[c] * basis_w[d] * basis_v[e];
                        dw += coeff * basis_x[a] * basis_y[b] * basis_z[c] * deriv_w[d] * basis_v[e];
                        dv += coeff * basis_x[a] * basis_y[b] * basis_z[c] * basis_w[d] * deriv_v[e];
                    }
                }
            }
        }
    }
}
