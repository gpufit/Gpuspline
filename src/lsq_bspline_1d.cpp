#include "spline_classes.h"
#include "libs/Eigen/Dense"

LSQ_BSpline_1D::LSQ_BSpline_1D()
    : num_data_points_(0), num_control_points_(0), num_knots_(0),
    coefficients_(nullptr), knots_(nullptr) {}

LSQ_BSpline_1D::LSQ_BSpline_1D(int num_data_points)
    : num_data_points_(num_data_points) {
    int num_internal_knots = num_data_points_ - 2;
    num_knots_ = 2 * (SPLINE_DEGREE + 1) + num_internal_knots;
    num_control_points_ = num_knots_ - (SPLINE_DEGREE + 1);

    coefficients_ = new REAL[num_control_points_];
    knots_ = new REAL[num_knots_];

    init_zero();
    init_uniform_clamped_knots();
}

LSQ_BSpline_1D::LSQ_BSpline_1D(int num_data_points, const REAL* values)
    : LSQ_BSpline_1D(num_data_points) {
    calculate_coefficients(values);
}

LSQ_BSpline_1D::~LSQ_BSpline_1D() {
    delete[] coefficients_;
    delete[] knots_;
}

void LSQ_BSpline_1D::init_zero() {
    std::memset(coefficients_, 0, sizeof(REAL) * num_control_points_);
}

void LSQ_BSpline_1D::init_uniform_clamped_knots() {
    int left_pad = SPLINE_DEGREE + 1;
    int right_pad = SPLINE_DEGREE + 1;
    int num_internal_knots = num_data_points_ - 2;
    int internal_start = 1;

    // Left clamped knots
    for (int i = 0; i < left_pad; ++i) {
        knots_[i] = static_cast<REAL>(0);
    }

    // Internal uniform knots
    for (int i = 0; i < num_internal_knots; ++i) {
        knots_[left_pad + i] = static_cast<REAL>(internal_start + i);
    }

    // Right clamped knots
    REAL end_value = static_cast<REAL>(num_data_points_ - 1);
    for (int i = num_knots_ - right_pad; i < num_knots_; ++i) {
        knots_[i] = end_value;
    }
}

void LSQ_BSpline_1D::calculate_coefficients(const REAL* values) {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int M = num_data_points_;
    const int N = num_control_points_;

    MatrixXd A = MatrixXd::Zero(M, N);
    VectorXd y = VectorXd::Zero(M);

    // Loop over each data point (x = 0, 1, ..., M-1)
    for (int i = 0; i < M; ++i) {
        REAL x = static_cast<REAL>(i);
        y(i) = static_cast<double>(values[i]);

        // Find the knot span index j such that knots[j] <= x < knots[j+1]
        int span = -1;
        for (int j = SPLINE_DEGREE; j < num_knots_ - SPLINE_DEGREE - 1; ++j) {
            if (x >= knots_[j] && x < knots_[j + 1]) {
                span = j;
                break;
            }
        }
        if (x == knots_[num_knots_ - 1]) {
            span = num_knots_ - SPLINE_DEGREE - 2;  // Clamp to last interval
        }
        assert(span >= SPLINE_DEGREE && span < num_control_points_);

        // Evaluate the nonzero basis functions at x
        REAL t = x - static_cast<int>(x);
        REAL basis[SPLINE_DEGREE + 1];
        //evaluate_fast_basis(t, basis);
        evaluate_basis_for_fitting(x, span, basis, knots_);

        // Fill row i of matrix A
        int start = span - SPLINE_DEGREE;
        for (int k = 0; k <= SPLINE_DEGREE; ++k) {
            int col = start + k;
            if (col >= 0 && col < N) {
                A(i, col) = basis[k];
            }
        }

    }

    std::cout << "\nDesign matrix A (" << A.rows() << " x " << A.cols() << "):\n";
    std::cout << A << "\n";

    std::cout << "Right-hand side y:\n";
    for (int i = 0; i < M; ++i) {
        std::cout << "  y[" << i << "] = " << y(i) << "\n";
    }

    // Solve the least-squares system: A * c = y
    VectorXd c = A.colPivHouseholderQr().solve(y);

    // Copy back to raw coefficient array
    for (int i = 0; i < N; ++i) {
        coefficients_[i] = static_cast<REAL>(c(i));
    }
}

void inline LSQ_BSpline_1D::evaluate_fast_basis(REAL t, REAL* basis) const {
    REAL t2 = t * t;
    REAL t3 = t2 * t;
    REAL omt = 1.0 - t;
    REAL omt2 = omt * omt;
    REAL omt3 = omt2 * omt;

    basis[0] = (1.0 / 6.0) * omt3;
    basis[1] = (1.0 / 6.0) * (3.0 * t3 - 6.0 * t2 + 4.0);
    basis[2] = (1.0 / 6.0) * (-3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0);
    basis[3] = (1.0 / 6.0) * t3;
}

void inline LSQ_BSpline_1D::evaluate_fast_basis_derivative(REAL t, REAL* dbasis) const {
    REAL t2 = t * t;
    REAL omt = 1.0 - t;
    REAL omt2 = omt * omt;

    dbasis[0] = -0.5 * omt2;
    dbasis[1] = 1.5 * t2 - 2.0 * t;
    dbasis[2] = -1.5 * t2 + t + 0.5;
    dbasis[3] = 0.5 * t2;
}

void LSQ_BSpline_1D::evaluate_basis_for_fitting(REAL x, int span, REAL* basis, const REAL* knots) {
    constexpr int DEGREE = 3;
    REAL left[DEGREE + 1];
    REAL right[DEGREE + 1];

    basis[0] = 1.0;

    for (int j = 1; j <= DEGREE; ++j) {
        left[j] = x - knots[span + 1 - j];
        right[j] = knots[span + j] - x;
        REAL saved = 0.0;

        for (int r = 0; r < j; ++r) {
            REAL temp = basis[r] / (right[r + 1] + left[j - r]);
            basis[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        basis[j] = saved;
    }
}

REAL inline LSQ_BSpline_1D::evaluate(REAL x) const {
    int i = static_cast<int>(std::floor(x));
    int span = std::min(i + SPLINE_DEGREE, num_control_points_ - 1);

    REAL t = x - static_cast<REAL>(i);  // local coordinate in [0, 1]
    REAL basis[SPLINE_DEGREE + 1];
    evaluate_fast_basis(t, basis);

    REAL result = 0.0;
    int first_basis = span - SPLINE_DEGREE;

    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        int index = first_basis + j;
        if (index >= 0 && index < num_control_points_) {
            result += coefficients_[index] * basis[j];
        }
    }
    return result;
}

REAL inline LSQ_BSpline_1D::evaluate_derivative(REAL x) const {
    int i = static_cast<int>(std::floor(x));
    int span = std::min(i + SPLINE_DEGREE, num_control_points_ - 1);

    REAL t = x - static_cast<REAL>(i);
    REAL dbasis[SPLINE_DEGREE + 1];
    evaluate_fast_basis_derivative(t, dbasis);

    REAL result = 0.0;
    int first_basis = span - SPLINE_DEGREE;

    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        int index = first_basis + j;
        if (index >= 0 && index < num_control_points_) {
            result += coefficients_[index] * dbasis[j];
        }
    }
    return result;
}

void LSQ_BSpline_1D::calculate_values(const REAL* x_values, int num_points, REAL* output_values) const {
    for (int i = 0; i < num_points; ++i) {
        output_values[i] = evaluate(x_values[i]);
    }
}

void LSQ_BSpline_1D::calculate_derivatives(const REAL* x_values, int num_points, REAL* output_derivatives) const {
    for (int i = 0; i < num_points; ++i) {
        output_derivatives[i] = evaluate_derivative(x_values[i]);
    }
}
