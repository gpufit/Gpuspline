
#include "natural_bspline_1d.h"
#include "libs/fitpack/gpuspline_fitpack_functions.h"
#include "bspline_fast_cubic_basis_evaluate.h"

Natural_BSpline_1D::Natural_BSpline_1D(int num_data_points)
    : num_data_points_(num_data_points), 
    num_control_points_(num_data_points + 2), 
    num_knots_(num_control_points_ + SPLINE_DEGREE + 1), 
    coefficients_(nullptr), 
    knots_(nullptr) {

    coefficients_ = new REAL[num_control_points_]{};
    knots_ = new REAL[num_knots_]{};

    use_fast_evaluation_ = true;

    knots_ready_ = false;
    coefficients_ready_ = false;

    init_knot_vector();          // Creates 3+N+3 = N+6 knots for a cubic spline

}

Natural_BSpline_1D::Natural_BSpline_1D(int num_data_points, const REAL* values)
    : num_data_points_(num_data_points),
    num_control_points_(num_data_points + 2),
    num_knots_(num_control_points_ + SPLINE_DEGREE + 1) {

    coefficients_ = new REAL[num_control_points_]{};
    knots_ = new REAL[num_knots_]{};

    use_fast_evaluation_ = true;

    knots_ready_ = false;
    coefficients_ready_ = false;

    init_knot_vector();          // Creates 3+N+3 = N+6 knots for a cubic spline
    calculate_coefficients(values);
}

Natural_BSpline_1D::~Natural_BSpline_1D() {
    delete[] coefficients_;
    delete[] knots_;
}

bool Natural_BSpline_1D::get_fast_evaluation() {
    return use_fast_evaluation_;
}

void Natural_BSpline_1D::set_fast_evaluation(bool flag) {
    use_fast_evaluation_ = flag;
}

void Natural_BSpline_1D::init_knot_vector() {
    // Create knot vector with _augknt logic: [x[0]]*3 + x + [x[-1]]*3
    // For cubic spline, this means 3 left and 3 right repeats

    // This assumes data points are uniformly spaced from 0 to N-1.
    // For more general spacing, adapt this loop accordingly.

    int i;
    int k = SPLINE_DEGREE;
    int N = num_data_points_;

    for (i = 0; i < k; ++i) {
        knots_[i] = 0.0f;  // Repeat first x-value (assumed 0.0)
    }

    for (i = 0; i < N; ++i) {
        knots_[k + i] = static_cast<REAL>(i);
    }

    for (i = 0; i < k; ++i) {
        knots_[k + N + i] = static_cast<REAL>(N - 1);  // Repeat last x-value
    }

    knots_ready_ = true;

}

int Natural_BSpline_1D::find_knot_interval(REAL x) const {
    int k = SPLINE_DEGREE;
    int N = num_data_points_;

    if (x <= 0.0)
        return k;
    if (x >= static_cast<REAL>(N - 1))
        return num_control_points_ - 1;

    return static_cast<int>(x) + k;
}


void Natural_BSpline_1D::evaluate_basis(const REAL x, const int span, const int derivative_order, REAL* result) const {

    int i;
    constexpr int degree = SPLINE_DEGREE;

    // Convert knot vector to double for FITPACK
    std::vector<double> double_knots(num_knots_);
    for (i = 0; i < num_knots_; ++i) {
        double_knots[i] = static_cast<double>(knots_[i]);
    }

    double input_x = static_cast<double>(x);

    // FITPACK requires a workspace of size 2k + 2
    double wrk[2 * degree + 2] = {};

    fitpack::_deBoor_D(
        double_knots.data(),
        input_x,
        degree,
        static_cast<int>(span),
        derivative_order,
        wrk
    );

    for (i = 0; i <= degree; ++i) {
        result[i] = static_cast<REAL>(wrk[i]);
    }
}


void Natural_BSpline_1D::evaluate_fast_cubic_basis(const REAL x, const int span, REAL* basis) const {
 
    REAL t = x - knots_[span];

    if (span == SPLINE_DEGREE) {
        // boundary case: first support interval
        evaluate_fast_cubic_basis_first_span(t, basis);
    }
    else if (span == SPLINE_DEGREE + 1) {
        // boundary case: second support interval
        evaluate_fast_cubic_basis_second_span(t, basis);
    }
    else if (span >= SPLINE_DEGREE + 2 && span <= num_control_points_ - 3) {
        // interior case
        evaluate_fast_cubic_basis_interior_span(t, basis);
    }
    else if (span == num_control_points_ - 2) {
        // boundary case near the end
        evaluate_fast_cubic_basis_second_last_span(t, basis);
    }
    else if (span == num_control_points_ - 1) {
        // last span
        evaluate_fast_cubic_basis_last_span(t, basis);
    }
    else {
        // handle error or unexpected span
        assert(false && "Invalid span in evaluate_fast_cubic_basis()");
    }
}


inline void Natural_BSpline_1D::evaluate_fast_cubic_basis_derivative(const REAL x, const int span, REAL* dbasis) const {

    REAL t = x - knots_[span];

    if (span == SPLINE_DEGREE) {
        evaluate_fast_cubic_basis_derivative_first_span(t, dbasis);
    }
    else if (span == SPLINE_DEGREE + 1) {
        evaluate_fast_cubic_basis_derivative_second_span(t, dbasis);
    }
    else if (span >= SPLINE_DEGREE + 2 && span <= num_control_points_ - 3) {
        evaluate_fast_cubic_basis_derivative_interior_span(t, dbasis);
    }
    else if (span == num_control_points_ - 2) {
        evaluate_fast_cubic_basis_derivative_second_last_span(t, dbasis);
    }
    else if (span == num_control_points_ - 1) {
        evaluate_fast_cubic_basis_derivative_last_span(t, dbasis);
    }
    else {
        // handle error or unexpected span
        assert(false && "Invalid span in evaluate_fast_cubic_basis_derivative()");
    }
}


REAL Natural_BSpline_1D::evaluate(REAL x) const {

    int i;
    int idx;

    assert(coefficients_ready_ == true);

    int span = find_knot_interval(x);

    REAL basis[4] = {};
    if (use_fast_evaluation_) {
        evaluate_fast_cubic_basis(x, span, basis);
    }
    else {
        evaluate_basis(x, span, 0, basis);
    }

    REAL result = 0.0;
    for (i = 0; i <= SPLINE_DEGREE; ++i) {
        idx = span - SPLINE_DEGREE + i;
        if (idx >= 0 && idx < num_control_points_) {
            result += basis[i] * coefficients_[idx];
        }
    }
    return result;
}

REAL Natural_BSpline_1D::evaluate_derivative(REAL x) const {

    int i;
    int idx;

    assert(coefficients_ready_ == true);

    int span = find_knot_interval(x);

    REAL dbasis[4] = {};
    if (use_fast_evaluation_) {
        evaluate_fast_cubic_basis_derivative(x, span, dbasis);
    }
    else {
        evaluate_basis(x, span, 1, dbasis);
    }

    REAL result = 0.0;
    for (i = 0; i <= SPLINE_DEGREE; ++i) {
        idx = span - SPLINE_DEGREE + i;
        if (idx >= 0 && idx < num_control_points_) {
            result += dbasis[i] * coefficients_[idx];
        }
    }

    return result;
}

void Natural_BSpline_1D::calculate_coefficients(const REAL* values) {

    int i, j;
    int col;
    int span;
    constexpr int p = SPLINE_DEGREE;
    const int N = num_data_points_;
    const int M = num_control_points_;

    // Resize and zero-out A_ and b_ in double precision
    A_ = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(M), static_cast<Eigen::Index>(M));
    b_ = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(M));

    // Fill interior rows with basis values at data points
    for (i = 0; i < N; ++i) {
        REAL x = static_cast<REAL>(i); // Assume data points at x = 0, 1, ..., N-1
        span = find_knot_interval(x);

        REAL basis[p + 1] = {};
        evaluate_basis(x, span, 0, basis);

        for (j = 0; j <= p; ++j) {
            col = span - p + j;
            if (col >= 0 && col < M) {
                A_(i, col) = static_cast<double>(basis[j]);
            }
        }

        b_(i) = static_cast<double>(values[i]);
    }

    // Add natural boundary condition at start (second derivative = 0)
    {
        REAL dd_basis[p + 1] = {};
        span = find_knot_interval(knots_[p]); // Evaluate at first interior knot
        //evaluate_basis_second_derivative(knots_[p], span, dd_basis);
        evaluate_basis(knots_[p], span, 2, dd_basis);
        for (j = 0; j <= p; ++j) {
            col = span - p + j;
            if (col >= 0 && col < M) {
                A_(N, col) = static_cast<double>(dd_basis[j]);
            }
        }
        b_(N) = 0.0;
    }

    // Add natural boundary condition at end (second derivative = 0)
    {
        REAL dd_basis[p + 1] = {};
        REAL x = static_cast<REAL>(N - 1);
        span = find_knot_interval(x); // Evaluate at last interior knot
        //evaluate_basis_second_derivative(knots_[M - 1], span, dd_basis);
        evaluate_basis(x, span, 2, dd_basis);
        for (j = 0; j <= p; ++j) {
            col = span - p + j;
            if (col >= 0 && col < M) {
                A_(N + 1, col) = static_cast<double>(dd_basis[j]);
            }
        }
        b_(N + 1) = 0.0;
    }

    // Solve linear system
    Eigen::VectorXd x = A_.colPivHouseholderQr().solve(b_);

    for (i = 0; i < M; ++i) {
        coefficients_[i] = static_cast<REAL>(x(i));
    }

    coefficients_ready_ = true;
}


void Natural_BSpline_1D::calculate_values(const REAL* x_values, int num_points, REAL* output_values) const {
    for (int i = 0; i < num_points; ++i) {
        output_values[i] = evaluate(x_values[i]);
    }
}

void Natural_BSpline_1D::calculate_derivatives(const REAL* x_values, int num_points, REAL* output_derivatives) const {
    for (int i = 0; i < num_points; ++i) {
        output_derivatives[i] = evaluate_derivative(x_values[i]);
    }
}

