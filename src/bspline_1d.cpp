#include "spline_classes.h"

const int BSpline1D::SPLINE_DEGREE = 3;

BSpline1D::BSpline1D(int num_control_points)
    : num_control_points_(num_control_points) {
    coefficients_ = new REAL[num_control_points_];
    init_zero();
    init_not_a_knot();
}

BSpline1D::BSpline1D(int num_control_points, const REAL* values)
    : num_control_points_(num_control_points) {
    coefficients_ = new REAL[num_control_points_];
    init_zero();
    init_not_a_knot();
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


void BSpline1D::evaluate_basis(REAL x, int interval, REAL* basis /* expected size: SPLINE_DEGREE + 1 */) const {

    // Zero-initialize all entries
    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        basis[j] = 0.0;
    }

    // Special case: last knot
    if (x == knots_.back()) {
        basis[SPLINE_DEGREE] = 1.0;
        return;
    }

    // Initialize zeroth-degree functions
    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        int idx = interval - SPLINE_DEGREE + j;
        if (idx >= 0 && idx + 1 < knots_.size() && knots_[idx] <= x && x < knots_[idx + 1]) {
            basis[j] = 1.0;
        }
    }

    // Cox–de Boor recursion
    for (int d = 1; d <= SPLINE_DEGREE; ++d) {
        for (int j = 0; j <= SPLINE_DEGREE - d; ++j) {
            int i = interval - SPLINE_DEGREE + j;

            REAL left = 0.0;
            REAL denom_left = knots_[i + d] - knots_[i];
            if (denom_left != 0.0) {
                left = (x - knots_[i]) / denom_left * basis[j];
            }

            REAL right = 0.0;
            REAL denom_right = knots_[i + d + 1] - knots_[i + 1];
            if (denom_right != 0.0) {
                right = (knots_[i + d + 1] - x) / denom_right * basis[j + 1];
            }

            basis[j] = left + right;
        }
    }
}


void BSpline1D::evaluate_basis_derivative(REAL x, int interval, REAL* dbasis) const {
    REAL left[SPLINE_DEGREE + 1];
    REAL right[SPLINE_DEGREE + 1];
    REAL ndu[SPLINE_DEGREE + 1][SPLINE_DEGREE + 1] = { 0 };
    REAL a[2][SPLINE_DEGREE + 1] = { 0 };

    ndu[0][0] = 1.0;

    for (int j = 1; j <= SPLINE_DEGREE; ++j) {
        left[j] = x - knots_[interval + 1 - j];
        right[j] = knots_[interval + j] - x;
        REAL saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            REAL temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        dbasis[j] = 0.0;
    }

    for (int r = 0; r <= SPLINE_DEGREE; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;

        for (int k = 1; k <= 1; ++k) {  // first derivative only
            REAL d = 0.0;
            int rk = r - k;
            int pk = SPLINE_DEGREE - k;
            int j1 = (rk >= 0) ? 0 : -rk;
            int j2 = (r - 1 <= pk) ? k : SPLINE_DEGREE - r;

            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) /
                    (knots_[interval + r + j + 1 - k] - knots_[interval + r - k + j]);
                d += a[s2][j] * ndu[r - k + j][pk];
            }
            dbasis[r] = d * SPLINE_DEGREE;
        }
    }
}


int BSpline1D::find_knot_interval(REAL x) const {
    int k = SPLINE_DEGREE;
    int n_knots = static_cast<int>(knots_.size());
    int max_interval = n_knots - k - 2;  // Last valid interval for evaluation

    // Estimate interval index using floor(x) + offset for repeated boundary knots
    int interval = static_cast<int>(std::floor(x)) + k;

    // Clamp to valid interval range
    if (interval < k) {
        interval = k;
    }
    else if (interval > max_interval) {
        interval = max_interval;
    }

    return interval;
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

        for (int j = 0; j <= SPLINE_DEGREE; ++j) {
            int col = interval - SPLINE_DEGREE + j;
            if (col >= 0 && col < n) {
                collocation(row, col) = basis[j];
            }
        }
    }

    // Solve collocation * coefficients = values
    Eigen::VectorXd rhs(n);
    for (int i = 0; i < n; ++i) rhs[i] = values[i];

    Eigen::VectorXd coeff = collocation.colPivHouseholderQr().solve(rhs);

    for (int i = 0; i < n; ++i) coefficients_[i] = coeff[i];
}


REAL BSpline1D::evaluate(REAL x) const {
    if (x < knots_.front() || x > knots_.back()) {
        return 0.0; // or handle extrapolation policy
    }

    // Use fast interval lookup assuming uniform spacing
    int interval = find_knot_interval(x);

    // Evaluate basis functions
    REAL basis[SPLINE_DEGREE + 1];
    evaluate_basis(x, interval, basis);

    // Accumulate weighted sum
    REAL result = 0.0;
    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        int idx = interval - SPLINE_DEGREE + j;
        if (idx >= 0 && idx < num_control_points_) {
            result += coefficients_[idx] * basis[j];
        }
    }

    return result;
}


REAL BSpline1D::evaluate_derivative(REAL x) const {
    if (x < knots_.front() || x > knots_.back()) {
        return 0.0; // or handle extrapolation policy
    }

    int interval = find_knot_interval(x);

    // Evaluate derivative basis functions
    REAL dbasis[SPLINE_DEGREE + 1];
    evaluate_basis_derivative(x, interval, dbasis);

    REAL result = 0.0;
    for (int j = 0; j <= SPLINE_DEGREE; ++j) {
        int idx = interval - SPLINE_DEGREE + j;
        if (idx >= 0 && idx < num_control_points_) {
            result += coefficients_[idx] * dbasis[j];
        }
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