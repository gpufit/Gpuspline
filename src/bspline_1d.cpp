#include "spline_classes.h"
#include "libs/Eigen/Dense"

BSpline1D::BSpline1D(int num_control_points)
    : num_control_points_(num_control_points) {
    coefficients_ = new REAL[num_control_points_];
    init_zero();
    init_not_a_knot();
    //init_uniform_clamped();
}

BSpline1D::BSpline1D(int num_control_points, const REAL* values)
    : num_control_points_(num_control_points) {
    coefficients_ = new REAL[num_control_points_];
    init_zero();
    init_not_a_knot();
    //init_uniform_clamped();
    calculate_coefficients(values);
}


void BSpline1D::init_zero() {
    for (int i = 0; i < num_control_points_; ++i) {
        coefficients_[i] = 0.0;
    }
}


BSpline1D::~BSpline1D() {
    delete[] coefficients_;
}


void BSpline1D::init_not_a_knot() {
    // Number of knots needed for a not-a-knot cubic spline: n + k + 1
    int n_knots = num_control_points_ + SPLINE_DEGREE + 1;
    knots_.resize(n_knots);

    // Step 1: Create uniform grid points: x = [0, 1, ..., n-1]
    std::vector<REAL> x(num_control_points_);
    for (int i = 0; i < num_control_points_; ++i) {
        x[i] = static_cast<REAL>(i);
    }

    // Step 2: Determine how many points to remove from each side (k2 = (k + 1) // 2)
    int k = SPLINE_DEGREE;
    int k2 = (k + 1) / 2; // Only valid for odd k

    // Step 3: Insert repeated knots at each boundary (k + 1 copies)
    for (int i = 0; i <= k; ++i) {
        knots_[i] = x[0];                        // Start of knot vector
        knots_[n_knots - 1 - i] = x.back();      // End of knot vector
    }

    // Step 4: Fill internal knots by skipping the second and second-to-last grid points
    for (int i = 0; i < num_control_points_ - 2 * k2; ++i) {
        knots_[i + k + 1] = x[i + k2];
    }
}



void BSpline1D::init_uniform_clamped() {
    int k = SPLINE_DEGREE;
    int n_knots = num_control_points_ + k + 1;
    knots_.resize(n_knots);

    for (int i = 0; i <= k; ++i) {
        knots_[i] = 0.0;  // Clamped left
    }
    for (int i = k + 1; i < num_control_points_; ++i) {
        knots_[i] = static_cast<REAL>(i - k);
    }
    for (int i = num_control_points_; i < n_knots; ++i) {
        knots_[i] = static_cast<REAL>(num_control_points_ - k);
    }
}

inline bool is_in_fast_region(int i, int num_control_points, int degree) {
    return (i >= degree) && (i < num_control_points - degree - 1);
}

inline void fast_uniform_basis(REAL t, REAL* basis) {
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

inline void fast_uniform_basis_derivative(REAL t, REAL* dbasis) {
    REAL t2 = t * t;
    REAL omt = 1.0 - t;
    REAL omt2 = omt * omt;

    dbasis[0] = -0.5 * omt2;
    dbasis[1] = 0.5 * (9.0 * t2 - 8.0 * t);
    dbasis[2] = -0.5 * (9.0 * t2 - 4.0 * t - 1.0);
    dbasis[3] = 0.5 * t2;
}


void BSpline1D::evaluate_basis(REAL x, int interval, REAL* basis) const {
    REAL left[SPLINE_DEGREE + 1];
    REAL right[SPLINE_DEGREE + 1];
    REAL ndu[SPLINE_DEGREE + 1][SPLINE_DEGREE + 1];

    ndu[0][0] = 1.0;

    for (int j = 1; j <= SPLINE_DEGREE; ++j) {
        left[j] = x - knots_[interval + 1 - j];
        right[j] = knots_[interval + j] - x;
        REAL saved = 0.0;
        for (int r = 0; r < j; ++r) {
            REAL denom = right[r + 1] + left[j - r];
            REAL temp = (denom == 0.0) ? 0.0 : ndu[r][j - 1] / denom;
            ndu[j][r] = denom;
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        basis[j] = ndu[j][SPLINE_DEGREE];
    }
}


void BSpline1D::evaluate_basis_derivative(REAL x, int interval, REAL* dbasis) const {
    const int degree = SPLINE_DEGREE;
    REAL basis_lo[degree];  // degree-2 basis functions

    // Step 1: Compute (degree - 1) basis functions
    for (int j = 0; j < degree; ++j) basis_lo[j] = 0.0;
    int lo_interval = interval;
    REAL left[degree + 1];
    REAL right[degree + 1];
    REAL ndu[degree + 1][degree + 1] = { 0 };
    ndu[0][0] = 1.0;

    for (int j = 1; j <= degree - 1; ++j) {
        left[j] = x - knots_[lo_interval + 1 - j];
        right[j] = knots_[lo_interval + j] - x;
        REAL saved = 0.0;
        for (int r = 0; r < j; ++r) {
            REAL denom = right[r + 1] + left[j - r];
            REAL temp = (denom == 0.0) ? 0.0 : ndu[r][j - 1] / denom;
            ndu[j][r] = denom;
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    for (int j = 0; j < degree; ++j) {
        basis_lo[j] = ndu[j][degree - 1];
    }

    // Step 2: Compute first derivatives
    for (int j = 0; j <= degree; ++j) {
        dbasis[j] = 0.0;
    }

    for (int j = 0; j < degree; ++j) {
        int i = interval - degree + j + 1;
        REAL denom = knots_[i + degree] - knots_[i];
        REAL factor = (denom != 0.0) ? degree / denom : 0.0;

        dbasis[j] -= factor * basis_lo[j];     // subtract lower
        dbasis[j + 1] += factor * basis_lo[j];     // add upper
    }
}




int BSpline1D::find_knot_interval(REAL x) const {
    int k = SPLINE_DEGREE;
    int n_knots = static_cast<int>(knots_.size());
    int n_intervals = n_knots - 1;

    if (x <= knots_[k]) return k;
    if (x >= knots_[n_knots - k - 1]) return n_knots - k - 2;

    // Binary search for t_i <= x < t_{i+1}
    for (int i = k; i < n_knots - k - 1; ++i) {
        if (x >= knots_[i] && x < knots_[i + 1]) {
            return i;
        }
    }

    return n_knots - k - 2;  // Fallback to the last safe interval
}


void BSpline1D::calculate_coefficients(const REAL* values) {
    int n = num_control_points_;
    Eigen::MatrixXd collocation = Eigen::MatrixXd::Zero(n, n);

    // Build collocation matrix: each row corresponds to one x position
    for (int row = 0; row < n; ++row) {
        REAL x = static_cast<REAL>(row);  // regular grid at x = 0, 1, ..., n-1
        int interval = find_knot_interval(x);

        REAL basis[SPLINE_DEGREE + 1];
        evaluate_basis(x, interval, basis);


        REAL row_sum = 0.0;
        for (int j = 0; j <= SPLINE_DEGREE; ++j) {
            int col = interval - SPLINE_DEGREE + j;
            if (col >= 0 && col < n) {
                collocation(row, col) = basis[j];
                row_sum += basis[j];
            }
        }
        std::cout << "Row " << row << " basis sum: " << row_sum << std::endl;
    }

    // Solve collocation * coefficients = values
    Eigen::VectorXd rhs(n);
    for (int i = 0; i < n; ++i) rhs[i] = values[i];

    std::cout << "\nDesign matrix A (" << collocation.rows() << " x " << collocation.cols() << "):\n";
    std::cout << collocation << "\n";

    std::cout << "Right-hand side y:\n";
    for (int i = 0; i < n; ++i) {
        std::cout << "  y[" << i << "] = " << rhs(i) << "\n";
    }

    Eigen::VectorXd coeff = collocation.colPivHouseholderQr().solve(rhs);

    for (int i = 0; i < n; ++i) coefficients_[i] = coeff[i];
}


//REAL BSpline1D::evaluate(REAL x) const {
//    if (x < knots_.front() || x > knots_.back()) {
//        return 0.0; // or handle extrapolation policy
//    }
//
//    // Use fast interval lookup assuming uniform spacing
//    int interval = find_knot_interval(x);
//
//    if (interval < SPLINE_DEGREE || interval >= num_control_points_) {
//        return 0.0; // or NAN, or throw, depending on your policy
//    }
//
//    // Evaluate basis functions
//    REAL basis[SPLINE_DEGREE + 1];
//    evaluate_basis(x, interval, basis);
//
//    // Accumulate weighted sum
//    REAL result = 0.0;
//    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
//        int idx = interval - SPLINE_DEGREE + j;
//        if (idx >= 0 && idx < num_control_points_) {
//            result += coefficients_[idx] * basis[j];
//        }
//    }
//
//    return result;
//}


REAL BSpline1D::evaluate(REAL x) const {
    int i = static_cast<int>(std::floor(x));
    REAL t = x - i;

    if (is_in_fast_region(i, num_control_points_, SPLINE_DEGREE)) {
        REAL basis[4];
        fast_uniform_basis(t, basis);

        REAL result = 0.0;
        for (int j = 0; j < 4; ++j) {
            result += coefficients_[i - 1 + j] * basis[j];
        }
        return result;
    }

    // Fallback to full Cox–de Boor
    int interval = find_knot_interval(x);
    REAL basis[SPLINE_DEGREE + 1];
    evaluate_basis(x, interval, basis);

    REAL result = 0.0;
    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        result += coefficients_[interval - SPLINE_DEGREE + j] * basis[j];
    }
    return result;
}


//REAL BSpline1D::evaluate_derivative(REAL x) const {
//    if (x < knots_.front() || x > knots_.back()) {
//        return 0.0; // or handle extrapolation policy
//    }
//
//    int interval = find_knot_interval(x);
//
//    if (interval < SPLINE_DEGREE || interval >= num_control_points_) {
//        return 0.0; // or NAN, or throw, depending on your policy
//    }
//
//    // Evaluate derivative basis functions
//    REAL dbasis[SPLINE_DEGREE + 1];
//    evaluate_basis_derivative(x, interval, dbasis);
//
//    REAL result = 0.0;
//    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
//        int idx = interval - SPLINE_DEGREE + j;
//        if (idx >= 0 && idx < num_control_points_) {
//            result += coefficients_[idx] * dbasis[j];
//        }
//    }
//
//    return result;
//}

REAL BSpline1D::evaluate_derivative(REAL x) const {
    int i = static_cast<int>(std::floor(x));
    REAL t = x - i;

    if (is_in_fast_region(i, num_control_points_, SPLINE_DEGREE)) {
        REAL dbasis[4];
        fast_uniform_basis_derivative(t, dbasis);

        REAL result = 0.0;
        for (int j = 0; j < 4; ++j) {
            result += coefficients_[i - 1 + j] * dbasis[j];
        }
        return result;
    }

    // Fallback to full derivative logic
    int interval = find_knot_interval(x);
    REAL dbasis[SPLINE_DEGREE + 1];
    evaluate_basis_derivative(x, interval, dbasis);

    REAL result = 0.0;
    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        result += coefficients_[interval - SPLINE_DEGREE + j] * dbasis[j];
    }
    return result;
}

void BSpline1D::calculate_values(const REAL* x_values, int num_points, REAL* output_values) const {
    for (int i = 0; i < num_points; ++i) {
        output_values[i] = evaluate(x_values[i]);
    }
}


void BSpline1D::calculate_derivatives(const REAL* x_values, int num_points, REAL* output_derivatives) const {
    for (int i = 0; i < num_points; ++i) {
        output_derivatives[i] = evaluate_derivative(x_values[i]);
    }
}