#include "spline_classes.h"


// Constructor (initializes an empty coefficient array)
BSpline1D::BSpline1D(int Nx) : 
    Nx_(Nx), 
    coefficients_(nullptr), 
    coefficients_calculated_(false), 
    data_initialized_(false) 
{
    allocate_knots();
}

// Constructor (accepts precomputed coefficients)
BSpline1D::BSpline1D(int Nx, REAL* coefficients) : 
    Nx_(Nx), 
    coefficients_(coefficients), 
    coefficients_calculated_(true), 
    data_initialized_(false) 
{
    allocate_knots();
}

// Destructor (frees allocated memory if needed)
BSpline1D::~BSpline1D() 
{
    delete[] knots_;
}

// Initialize spline with raw data
void BSpline1D::initialize(const REAL* data) 
{
    if (!data) {
        data_initialized_ = false;
        return;
    }
    data_ = data;
    data_initialized_ = true;
}


// Allocate and initialize the knot vector
void BSpline1D::allocate_knots() {
    int num_knots = Nx_ + 6;  // Correct number of knots
    knots_ = new REAL[num_knots];

    // First four knots (repeated at 0)
    for (int i = 0; i < 4; i++) {
        knots_[i] = 0;
    }

    // Internal knots (evenly spaced from 1 to Nx-2)
    for (int i = 4; i < Nx_ + 2; i++) {
        knots_[i] = i - 3;
    }

    // Last four knots (repeated at Nx-1)
    for (int i = Nx_ + 2; i < num_knots; i++) {
        knots_[i] = Nx_ - 1;
    }
}


// Compute spline coefficients with natural boundary conditions
void BSpline1D::calculate_coefficients(REAL* coefficients) {
    if (!data_initialized_ || !coefficients) {
        return;
    }

    for (int x = 0; x < Nx_; x++) {
        int index = x * 4;

        // Function value at control point
        coefficients[index + 0] = data_[x];

        // First derivative approximation
        if (x == 0) {
            coefficients[index + 1] = data_[1] - data_[0]; // Forward difference
        }
        else if (x == Nx_ - 1) {
            coefficients[index + 1] = data_[Nx_ - 1] - data_[Nx_ - 2]; // Backward difference
        }
        else {
            coefficients[index + 1] = (data_[x + 1] - data_[x - 1]) / 2.0; // Central difference
        }

        // Second derivative approximation (Natural boundary conditions)
        if (x == 0 || x == Nx_ - 1) {
            coefficients[index + 2] = 0.0; // Enforce natural boundary
        }
        else {
            coefficients[index + 2] = (data_[x + 1] - 2.0 * data_[x] + data_[x - 1]);
        }

        // Third derivative approximation
        if (x < 2 || x > Nx_ - 3) {
            coefficients[index + 3] = 0.0; // Avoid undefined behavior near boundaries
        }
        else {
            coefficients[index + 3] = (data_[x + 2] - 3.0 * data_[x + 1] + 3.0 * data_[x] - data_[x - 1]);
        }
    }

    coefficients_calculated_ = true;
}


// Compute basis functions using explicit knots
void BSpline1D::compute_basis_functions(REAL* basis_values, REAL x) const {
    int i = static_cast<int>(std::floor(x));
    while (i > 0 && x < knots_[i]) {
        i--;
    }
    REAL t = (x - knots_[i]) / (knots_[i + 1] - knots_[i]);
    basis_values[0] = (REAL(1) / 6) * (1 - t) * (1 - t) * (1 - t);
    basis_values[1] = (REAL(1) / 6) * (3 * t * t * t - 6 * t * t + 4);
    basis_values[2] = (REAL(1) / 6) * (-3 * t * t * t + 3 * t * t + 3 * t + 1);
    basis_values[3] = (REAL(1) / 6) * (t * t * t);
}

// Compute cubic B-spline derivative basis functions for a given position
void BSpline1D::compute_derivative_basis_functions(REAL* basis_derivatives, REAL x) const {
    REAL t = x - std::floor(x);
    basis_derivatives[0] = (-0.5f) * (1 - t) * (1 - t);
    basis_derivatives[1] = (1.5f * t * t - 2 * t);
    basis_derivatives[2] = (-1.5f * t * t + t + 0.5f);
    basis_derivatives[3] = (0.5f * t * t);
}

// Evaluate the B-spline function
REAL BSpline1D::evaluate(REAL x) const {
    if (!coefficients_calculated_) {
        return 0.0;
    }

    REAL basis_values[4];
    compute_basis_functions(basis_values, x);

    int i = static_cast<int>(std::floor(x));
    if (i < 0) {
        i = 0;
    }
    else if (i > N_control_ - 4) { // Ensure at least 4 control points exist
        i = N_control_ - 4;
    }

    REAL value = 0.0;
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) { // Iterate over all coefficients
            int index = (i + a) * 4 + b;
            value += coefficients_[index] * basis_values[a];
        }
    }

    return value;
}


// Evaluate the B-spline function and its first derivative
void BSpline1D::evaluate_with_derivative(REAL x, REAL& function_value, REAL& dx) const {
    if (!coefficients_calculated_) {
        function_value = dx = 0.0;
        return;
    }

    REAL basis_values[4], basis_derivatives[4];
    compute_basis_functions(basis_values, x);
    compute_derivative_basis_functions(basis_derivatives, x);

    int i = static_cast<int>(x);
    i = std::max(0, std::min(i, Nx_ - 2));

    function_value = dx = 0.0;
    for (int a = 0; a < 4; ++a) {
        int index = (i + a) * 4;
        function_value += coefficients_[index] * basis_values[a];
        dx += coefficients_[index] * basis_derivatives[a];
    }
}

// calls calculate_value() iteratively
void BSpline1D::calculate_values(
    REAL * bspline_values,
    REAL const * x_values,
    std::size_t const size_x)
{
    for (std::size_t x_index = 0; x_index < size_x; x_index++)
    {
        bspline_values[x_index] = evaluate(x_values[x_index]);
    }
}
