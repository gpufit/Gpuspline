
#include <numeric>
#include "natural_bspline_1d.h"
#include "natural_bspline_nd.h"
#include "libs/fitpack/gpuspline_fitpack_functions.h"
#include "bspline_fast_cubic_basis_evaluate.h"

Natural_BSpline_ND::Natural_BSpline_ND(const std::vector<int>& data_dims, const REAL* data_values)
    : num_dims_(static_cast<int>(data_dims.size())), data_dims_(data_dims) {

    assert(num_dims_ >= 1);
    assert(data_values != nullptr);

    use_fast_evaluation_ = false;

    control_point_dims_.resize(num_dims_);
    for (int i = 0; i < num_dims_; ++i) {
        control_point_dims_[i] = data_dims[i] + 2;
    }

    total_data_points_ = std::accumulate(data_dims_.begin(), data_dims_.end(), 1, std::multiplies<int>());
    total_control_points_ = std::accumulate(control_point_dims_.begin(), control_point_dims_.end(), 1, std::multiplies<int>());

    // Initialize stride vectors
    init_stride_vectors();

    // Initialize multi-indices
    build_multi_indices();

    // Initialise knot vectors for each dimension
    init_knot_vectors();

    // Compute spline coefficients from input data
    calculate_coefficients(data_values);
}

Natural_BSpline_ND::Natural_BSpline_ND(const std::vector<int>& data_dims,
                                       const std::vector<REAL>& coefficients,
                                       const std::vector<int>& multi_indices)
    : data_dims_(data_dims),
    num_dims_(static_cast<int>(data_dims.size())),
    multi_indices_flat_(multi_indices),
    coefficients_(coefficients)
{
    constexpr int degree = SPLINE_DEGREE;

    use_fast_evaluation_ = false;

    // Reconstruct control_point_dims_ from total coefficient count
    control_point_dims_ = data_dims_;  // assume initially same size
    int total_expected = 1;
    for (int d = 0; d < num_dims_; ++d) {
        int n_ctrl = data_dims_[d] + degree - 1; // natural cubic spline control points
        control_point_dims_[d] = n_ctrl;
        total_expected *= n_ctrl;
    }

    assert(static_cast<int>(coefficients.size()) == total_expected);
    assert(static_cast<int>(multi_indices_flat_.size()) % num_dims_ == 0);

    total_data_points_ = std::accumulate(data_dims_.begin(), data_dims_.end(), 1, std::multiplies<int>());
    total_control_points_ = std::accumulate(control_point_dims_.begin(), control_point_dims_.end(), 1, std::multiplies<int>());
    num_combinations_ = static_cast<int>(multi_indices_flat_.size()) / num_dims_;

    // Initialize stride vectors
    init_stride_vectors();

    // Initialise knot vectors for each dimension
    init_knot_vectors();
}

Natural_BSpline_ND::~Natural_BSpline_ND() {

}

bool Natural_BSpline_ND::get_fast_evaluation() {
    return use_fast_evaluation_;
}

void Natural_BSpline_ND::set_fast_evaluation(bool flag) {
    use_fast_evaluation_ = flag;
}

void Natural_BSpline_ND::init_stride_vectors() {

    // Compute strides for flattened indexing
    data_strides_.resize(num_dims_);
    data_strides_[0] = 1;
    for (int i = 1; i < num_dims_; ++i)
        data_strides_[i] = data_strides_[i - 1] * data_dims_[i - 1];

    control_point_strides_.resize(num_dims_);
    control_point_strides_[0] = 1;
    for (int i = 1; i < num_dims_; ++i)
        control_point_strides_[i] = control_point_strides_[i - 1] * control_point_dims_[i - 1];

}

void Natural_BSpline_ND::build_multi_indices() {
    constexpr int degree = SPLINE_DEGREE;
    const int num_dims = num_dims_;

    num_combinations_ = static_cast<int>(std::pow(degree + 1, num_dims));
    multi_indices_flat_.resize(num_combinations_ * num_dims);

    for (int i = 0; i < num_combinations_; ++i) {
        int val = i;
        for (int d = 0; d < num_dims_; ++d) {  // now row-major
            multi_indices_flat_[i * num_dims_ + d] = val % (degree + 1);
            val /= (degree + 1);
        }
    }
}

void Natural_BSpline_ND::init_knot_vectors() {
    constexpr int k = SPLINE_DEGREE;

    knot_vectors_.resize(num_dims_);

    for (int d = 0; d < num_dims_; ++d) {
        int N = data_dims_[d]; // number of data points in dimension d
        int num_knots = N + 2 * k; // e.g., 10 + 6 = 16
        std::vector<REAL> knots(num_knots);

        // First k knots = 0.0
        for (int i = 0; i < k; ++i) {
            knots[i] = 0.0;
        }

        // Internal knots = 1, 2, ..., N - 1
        for (int i = 0; i < N; ++i) {
            knots[k + i] = static_cast<REAL>(i);
        }

        // Last k knots = N - 1
        for (int i = 0; i < k; ++i) {
            knots[k + N + i] = static_cast<REAL>(N - 1);
        }

        knot_vectors_[d] = std::move(knots);
    }
}

int Natural_BSpline_ND::find_knot_interval(REAL x, int dimension) const {
    int k = SPLINE_DEGREE;
    int N = data_dims_[dimension];
    int M = control_point_dims_[dimension];

    if (x <= 0.0)
        return k;
    if (x >= static_cast<REAL>(N - 1))
        return M - 1;

    return static_cast<int>(x) + k;
}

inline void Natural_BSpline_ND::evaluate_basis_fitpack(REAL x, const std::vector<REAL>& knots, int span, int derivative_order, REAL* basis) {
    constexpr int degree = SPLINE_DEGREE;

    std::vector<double> d_knots(knots.size());
    for (int i = 0; i < knots.size(); ++i) {
        d_knots[i] = static_cast<double>(knots[i]);
    }

    double wrk[2 * degree + 2] = {};

    fitpack::_deBoor_D(
        d_knots.data(),
        static_cast<double>(x),
        degree,
        span,
        derivative_order,
        wrk
    );

    for (int i = 0; i <= degree; ++i) {
        basis[i] = static_cast<REAL>(wrk[i]);
    }
}

void Natural_BSpline_ND::evaluate_fast_cubic_basis(const REAL x, const int span, const int axis, REAL* basis) const {

    std::vector<REAL> knots = knot_vectors_[axis];
    int num_control_points = control_point_dims_[axis];

    REAL t = x - knots[span];

    if (span == SPLINE_DEGREE) {
        // boundary case: first support interval
        evaluate_fast_cubic_basis_first_span(t, basis);
    }
    else if (span == SPLINE_DEGREE + 1) {
        // boundary case: second support interval
        evaluate_fast_cubic_basis_second_span(t, basis);
    }
    else if (span >= SPLINE_DEGREE + 2 && span <= num_control_points - 3) {
        // interior case
        evaluate_fast_cubic_basis_interior_span(t, basis);
    }
    else if (span == num_control_points - 2) {
        // boundary case near the end
        evaluate_fast_cubic_basis_second_last_span(t, basis);
    }
    else if (span == num_control_points - 1) {
        // last span
        evaluate_fast_cubic_basis_last_span(t, basis);
    }
    else {
        // handle error or unexpected span
        assert(false && "Invalid span in evaluate_fast_cubic_basis()");
    }
}

inline void Natural_BSpline_ND::evaluate_fast_cubic_basis_derivative(const REAL x, const int span, const int axis, REAL* dbasis) const {

    std::vector<REAL> knots = knot_vectors_[axis];
    int num_control_points = control_point_dims_[axis];

    REAL t = x - knots[span];

    if (span == SPLINE_DEGREE) {
        evaluate_fast_cubic_basis_derivative_first_span(t, dbasis);
    }
    else if (span == SPLINE_DEGREE + 1) {
        evaluate_fast_cubic_basis_derivative_second_span(t, dbasis);
    }
    else if (span >= SPLINE_DEGREE + 2 && span <= num_control_points - 3) {
        evaluate_fast_cubic_basis_derivative_interior_span(t, dbasis);
    }
    else if (span == num_control_points - 2) {
        evaluate_fast_cubic_basis_derivative_second_last_span(t, dbasis);
    }
    else if (span == num_control_points - 1) {
        evaluate_fast_cubic_basis_derivative_last_span(t, dbasis);
    }
    else {
        // handle error or unexpected span
        assert(false && "Invalid span in evaluate_fast_cubic_basis_derivative()");
    }
}

void Natural_BSpline_ND::calculate_coefficients(const REAL* values) {
    // Initialize
    coefficients_.assign(total_control_points_, 0.0);

    std::vector<int> current_data_shape = data_dims_;
    std::vector<REAL> temp_values(values, values + total_data_points_);

    int ndim = num_dims_;

    // Local flatten function assuming C-style (x-fastest) layout
    auto local_flatten = [](int i, int slice_idx, int axis, const std::vector<int>& shape) -> int {
        int ndim = shape.size();
        std::vector<int> strides(ndim);
        strides[0] = 1;
        for (int d = 1; d < ndim; ++d) {
            strides[d] = strides[d - 1] * shape[d - 1];
        }

        int flat_idx = i * strides[axis];
        int remaining = slice_idx;
        for (int d = 0; d < ndim; ++d) {
            if (d == axis) continue;
            int coord = remaining % shape[d];
            flat_idx += coord * strides[d];
            remaining /= shape[d];
        }
        return flat_idx;
    };

    // Solve along each axis
    for (int axis = 0; axis < ndim; ++axis) {
        int n_data = current_data_shape[axis];
        int n_control = control_point_dims_[axis];

        //std::cout << "Starting solve on axis " << axis << ", current shape: ";
        //for (auto v : current_data_shape) std::cout << v << " ";
        //std::cout << "\n";

        int num_slices = 1;
        for (int d = 0; d < ndim; ++d) {
            if (d != axis) num_slices *= current_data_shape[d];
        }

        std::vector<int> temp_output_shape = current_data_shape;
        temp_output_shape[axis] = n_control;

        std::vector<REAL> new_coeffs(num_slices * n_control, 0.0);

        Natural_BSpline_1D spline_1d(n_data);

        for (int slice_idx = 0; slice_idx < num_slices; ++slice_idx) {
            // Extract 1D slice
            std::vector<REAL> slice(n_data);
            for (int i = 0; i < n_data; ++i) {
                int idx = local_flatten(i, slice_idx, axis, current_data_shape);
                slice[i] = temp_values[idx];
            }

            //std::cout << "Slice " << slice_idx << " along axis " << axis << ": ";
            //for (auto v : slice) std::cout << v << " ";
            //std::cout << "\n";

            // Solve 1D spline
            spline_1d.calculate_coefficients(slice.data());

            // Copy solved coefficients
            for (int i = 0; i < n_control; ++i) {
                int idx = local_flatten(i, slice_idx, axis, temp_output_shape);
                new_coeffs[idx] = spline_1d.coefficients()[i];
            }
        }

        // Update working state
        temp_values = std::move(new_coeffs);
        current_data_shape = temp_output_shape;
    }

    // Save final coefficients
    coefficients_ = std::move(temp_values);
}

REAL Natural_BSpline_ND::evaluate(const REAL* pt) const {
    if (use_fast_evaluation_) {
        switch (num_dims_) {
        case 1: return evaluate_1d_fast(pt);
        case 2: return evaluate_2d_fast(pt);
        case 3: return evaluate_3d_fast(pt);
        case 4: return evaluate_4d_fast(pt);
        case 5: return evaluate_5d_fast(pt);
        case 6: return evaluate_6d_fast(pt);
        default: return evaluate_general(pt);  // fallback for D > 5
        }
    }
    else {
        return evaluate_general(pt);  // always use slow path
    }
}

REAL Natural_BSpline_ND::evaluate_general(const REAL* pt) const {
    constexpr int degree = SPLINE_DEGREE;
    const int num_dims = num_dims_;

    std::vector<std::array<REAL, degree + 1>> basis_values(num_dims);
    std::vector<int> spans(num_dims);

    // Step 1: Compute spans and basis functions in each dimension
    for (int dim = 0; dim < num_dims; ++dim) {
        const REAL x_dim = pt[dim];
        const auto& knots = knot_vectors_[dim];

        int span = find_knot_interval(x_dim, dim);
        spans[dim] = span;

        if (use_fast_evaluation_) {
            evaluate_fast_cubic_basis(x_dim, span, dim, basis_values[dim].data());
        }
        else {
            evaluate_basis_fitpack(x_dim, knots, span, 0, basis_values[dim].data());
        }

    }

    // Step 2: Tensor-product sum
    REAL result = 0.0;

    std::vector<int> ctrl_offset(num_dims);
    for (int dim = 0; dim < num_dims; ++dim) {
        ctrl_offset[dim] = spans[dim] - degree;
    }

    const std::vector<int>& strides = control_point_strides_;

    for (int i = 0; i < num_combinations_; ++i) {
        REAL weight = 1.0;
        int coeff_idx = 0;

        for (int dim = 0; dim < num_dims; ++dim) {
            int basis_idx = multi_indices_flat_[i * num_dims_ + dim];
            weight *= basis_values[dim][basis_idx];
            int idx = ctrl_offset[dim] + basis_idx;
            coeff_idx += idx * strides[dim];
        }

        result += weight * coefficients_[coeff_idx];
    }

    return result;
}

REAL Natural_BSpline_ND::evaluate_1d_fast(const REAL* pt) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<REAL, span_size> B;
    int span = find_knot_interval(pt[0], 0);
    evaluate_fast_cubic_basis(pt[0], span, 0, B.data());

    int offset = span - degree;
    int stride = control_point_strides_[0];

    REAL result = 0.0;
    for (int i = 0; i < span_size; ++i) {
        int idx = (offset + i) * stride;
        result += B[i] * coefficients_[idx];
    }

    return result;
}

REAL Natural_BSpline_ND::evaluate_2d_fast(const REAL* pt) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 2> B;
    std::array<int, 2> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 2; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;
        evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
    }

    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];

    // Step 2: tensor-product sum
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j) {
            REAL w = B[0][i] * B[1][j];
            int idx = (o0 + i) * s0 + (o1 + j) * s1;
            result += w * coefficients_[idx];
        }

    return result;
}

REAL Natural_BSpline_ND::evaluate_3d_fast(const REAL* pt) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 3> B;
    std::array<int, 3> spans;

    for (int d = 0; d < 3; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;
        evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
    }

    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;
    int o2 = spans[2] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];
    const int s2 = control_point_strides_[2];

    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j)
            for (int k = 0; k < span_size; ++k) {
                REAL w = B[0][i] * B[1][j] * B[2][k];
                int idx = (o0 + i) * s0 + (o1 + j) * s1 + (o2 + k) * s2;
                result += w * coefficients_[idx];
            }

    return result;
}

REAL Natural_BSpline_ND::evaluate_4d_fast(const REAL* pt) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 4> B;
    std::array<int, 4> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 4; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;
        evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
    }

    // Step 2: offset and strides
    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;
    int o2 = spans[2] - degree;
    int o3 = spans[3] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];
    const int s2 = control_point_strides_[2];
    const int s3 = control_point_strides_[3];

    // Step 3: tensor-product sum
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j)
            for (int k = 0; k < span_size; ++k)
                for (int l = 0; l < span_size; ++l) {
                    REAL w = B[0][i] * B[1][j] * B[2][k] * B[3][l];
                    int idx = (o0 + i) * s0 + (o1 + j) * s1 + (o2 + k) * s2 + (o3 + l) * s3;
                    result += w * coefficients_[idx];
                }

    return result;
}

REAL Natural_BSpline_ND::evaluate_5d_fast(const REAL* pt) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 5> B;
    std::array<int, 5> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 5; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;
        evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
    }

    // Step 2: offset and strides
    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;
    int o2 = spans[2] - degree;
    int o3 = spans[3] - degree;
    int o4 = spans[4] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];
    const int s2 = control_point_strides_[2];
    const int s3 = control_point_strides_[3];
    const int s4 = control_point_strides_[4];

    // Step 3: tensor-product sum
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j)
            for (int k = 0; k < span_size; ++k)
                for (int l = 0; l < span_size; ++l)
                    for (int m = 0; m < span_size; ++m) {
                        REAL w = B[0][i] * B[1][j] * B[2][k] * B[3][l] * B[4][m];
                        int idx = (o0 + i) * s0 + (o1 + j) * s1 + (o2 + k) * s2 +
                            (o3 + l) * s3 + (o4 + m) * s4;
                        result += w * coefficients_[idx];
                    }

    return result;
}

REAL Natural_BSpline_ND::evaluate_6d_fast(const REAL* pt) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 6> B;
    std::array<int, 6> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 6; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;
        evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
    }

    // Step 2: offsets and strides
    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;
    int o2 = spans[2] - degree;
    int o3 = spans[3] - degree;
    int o4 = spans[4] - degree;
    int o5 = spans[5] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];
    const int s2 = control_point_strides_[2];
    const int s3 = control_point_strides_[3];
    const int s4 = control_point_strides_[4];
    const int s5 = control_point_strides_[5];

    // Step 3: tensor-product sum
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j)
            for (int k = 0; k < span_size; ++k)
                for (int l = 0; l < span_size; ++l)
                    for (int m = 0; m < span_size; ++m)
                        for (int n = 0; n < span_size; ++n) {
                            REAL w = B[0][i] * B[1][j] * B[2][k] * B[3][l] * B[4][m] * B[5][n];
                            int idx = (o0 + i) * s0 + (o1 + j) * s1 + (o2 + k) * s2 +
                                (o3 + l) * s3 + (o4 + m) * s4 + (o5 + n) * s5;
                            result += w * coefficients_[idx];
                        }

    return result;
}

REAL Natural_BSpline_ND::evaluate_derivative(const REAL* pt, int d_target) const {
    if (use_fast_evaluation_) {
        switch (num_dims_) {
        case 1: return evaluate_derivative_1d_fast(pt, d_target);
        case 2: return evaluate_derivative_2d_fast(pt, d_target);
        case 3: return evaluate_derivative_3d_fast(pt, d_target);
        case 4: return evaluate_derivative_4d_fast(pt, d_target);
        case 5: return evaluate_derivative_5d_fast(pt, d_target);
        case 6: return evaluate_derivative_6d_fast(pt, d_target);
        default: return evaluate_derivative_general(pt, d_target);  // fallback
        }
    }
    else {
        return evaluate_derivative_general(pt, d_target);
    }
}

REAL Natural_BSpline_ND::evaluate_derivative_general(const REAL* pt, int d_target) const {
    constexpr int degree = SPLINE_DEGREE;
    const int num_dims = num_dims_;

    std::vector<std::array<REAL, degree + 1>> basis_values(num_dims);
    std::array<REAL, degree + 1> derivative_values;
    std::vector<int> spans(num_dims);

    // Step 1: Evaluate basis and derivative basis functions
    for (int dim = 0; dim < num_dims; ++dim) {

        const REAL x_dim = pt[dim];
        const auto& knots = knot_vectors_[dim];

        int span = find_knot_interval(x_dim, dim);
        spans[dim] = span;

        if (dim == d_target) {
            if (use_fast_evaluation_) {
                evaluate_fast_cubic_basis_derivative(x_dim, span, dim, derivative_values.data());
            }
            else {
                evaluate_basis_fitpack(x_dim, knots, span, 1, derivative_values.data());
            }
        }
        else {
            if (use_fast_evaluation_) {
                evaluate_fast_cubic_basis(x_dim, span, dim, basis_values[dim].data());
            }
            else {
                evaluate_basis_fitpack(x_dim, knots, span, 0, basis_values[dim].data());
            }
        }
    }

    // Step 2: Tensor-product sum
    REAL result = 0.0;

    std::vector<int> ctrl_offset(num_dims);
    for (int dim = 0; dim < num_dims; ++dim) {
        ctrl_offset[dim] = spans[dim] - degree;
    }

    const std::vector<int>& strides = control_point_strides_;

    for (int i = 0; i < num_combinations_; ++i) {
        REAL weight = 1.0;
        int coeff_idx = 0;

        for (int dim = 0; dim < num_dims; ++dim) {
            int basis_idx = multi_indices_flat_[i * num_dims_ + dim];
            int offset = ctrl_offset[dim] + basis_idx;
            coeff_idx += offset * strides[dim];

            if (dim == d_target) {
                weight *= derivative_values[basis_idx];
            }
            else {
                weight *= basis_values[dim][basis_idx];
            }
        }

        result += weight * coefficients_[coeff_idx];
    }

    return result;
}

REAL Natural_BSpline_ND::evaluate_derivative_1d_fast(const REAL* pt, int d_target) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<REAL, span_size> dB;
    int span = find_knot_interval(pt[0], 0);
    evaluate_fast_cubic_basis_derivative(pt[0], span, 0, dB.data());

    int offset = span - degree;
    int stride = control_point_strides_[0];

    REAL result = 0.0;
    for (int i = 0; i < span_size; ++i) {
        int idx = (offset + i) * stride;
        result += dB[i] * coefficients_[idx];
    }

    return result;
}

REAL Natural_BSpline_ND::evaluate_derivative_2d_fast(const REAL* pt, int d_target) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 2> B;
    std::array<REAL, span_size> dB;
    std::array<int, 2> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 2; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;

        if (d == d_target) {
            evaluate_fast_cubic_basis_derivative(pt[d], span, d, dB.data());
        }
        else {
            evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
        }
    }

    // Step 2: offsets and strides
    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];

    // Step 3: tensor-product sum
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j) {
            REAL w = 1.0;
            w *= (d_target == 0) ? dB[i] : B[0][i];
            w *= (d_target == 1) ? dB[j] : B[1][j];

            int idx = (o0 + i) * s0 + (o1 + j) * s1;
            result += w * coefficients_[idx];
        }

    return result;
}

REAL Natural_BSpline_ND::evaluate_derivative_3d_fast(const REAL* pt, int d_target) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 3> B;
    std::array<REAL, span_size> dB;
    std::array<int, 3> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 3; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;

        if (d == d_target) {
            evaluate_fast_cubic_basis_derivative(pt[d], span, d, dB.data());
        }
        else {
            evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
        }
    }

    // Step 2: offsets and strides
    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;
    int o2 = spans[2] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];
    const int s2 = control_point_strides_[2];

    // Step 3: tensor-product sum
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j)
            for (int k = 0; k < span_size; ++k) {
                REAL w = 1.0;
                w *= (d_target == 0) ? dB[i] : B[0][i];
                w *= (d_target == 1) ? dB[j] : B[1][j];
                w *= (d_target == 2) ? dB[k] : B[2][k];

                int idx = (o0 + i) * s0 + (o1 + j) * s1 + (o2 + k) * s2;
                result += w * coefficients_[idx];
            }

    return result;
}

REAL Natural_BSpline_ND::evaluate_derivative_4d_fast(const REAL* pt, int d_target) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 4> B;
    std::array<REAL, span_size> dB;
    std::array<int, 4> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 4; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;

        if (d == d_target) {
            evaluate_fast_cubic_basis_derivative(pt[d], span, d, dB.data());
        }
        else {
            evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
        }
    }

    // Step 2: offsets and strides
    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;
    int o2 = spans[2] - degree;
    int o3 = spans[3] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];
    const int s2 = control_point_strides_[2];
    const int s3 = control_point_strides_[3];

    // Step 3: tensor-product sum
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j)
            for (int k = 0; k < span_size; ++k)
                for (int l = 0; l < span_size; ++l) {
                    REAL w = 1.0;
                    w *= (d_target == 0) ? dB[i] : B[0][i];
                    w *= (d_target == 1) ? dB[j] : B[1][j];
                    w *= (d_target == 2) ? dB[k] : B[2][k];
                    w *= (d_target == 3) ? dB[l] : B[3][l];

                    int idx = (o0 + i) * s0 + (o1 + j) * s1 +
                        (o2 + k) * s2 + (o3 + l) * s3;

                    result += w * coefficients_[idx];
                }

    return result;
}

REAL Natural_BSpline_ND::evaluate_derivative_5d_fast(const REAL* pt, int d_target) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 5> B;
    std::array<REAL, span_size> dB;
    std::array<int, 5> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 5; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;

        if (d == d_target) {
            evaluate_fast_cubic_basis_derivative(pt[d], span, d, dB.data());
        }
        else {
            evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
        }
    }

    // Step 2: offsets and strides
    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;
    int o2 = spans[2] - degree;
    int o3 = spans[3] - degree;
    int o4 = spans[4] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];
    const int s2 = control_point_strides_[2];
    const int s3 = control_point_strides_[3];
    const int s4 = control_point_strides_[4];

    // Step 3: tensor-product sum
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j)
            for (int k = 0; k < span_size; ++k)
                for (int l = 0; l < span_size; ++l)
                    for (int m = 0; m < span_size; ++m) {
                        REAL w = 1.0;
                        w *= (d_target == 0) ? dB[i] : B[0][i];
                        w *= (d_target == 1) ? dB[j] : B[1][j];
                        w *= (d_target == 2) ? dB[k] : B[2][k];
                        w *= (d_target == 3) ? dB[l] : B[3][l];
                        w *= (d_target == 4) ? dB[m] : B[4][m];

                        int idx = (o0 + i) * s0 + (o1 + j) * s1 + (o2 + k) * s2 +
                            (o3 + l) * s3 + (o4 + m) * s4;

                        result += w * coefficients_[idx];
                    }

    return result;
}

REAL Natural_BSpline_ND::evaluate_derivative_6d_fast(const REAL* pt, int d_target) const {
    constexpr int degree = SPLINE_DEGREE;
    constexpr int span_size = degree + 1;

    std::array<std::array<REAL, span_size>, 6> B;
    std::array<REAL, span_size> dB;
    std::array<int, 6> spans;

    // Step 1: compute spans and basis values
    for (int d = 0; d < 6; ++d) {
        int span = find_knot_interval(pt[d], d);
        spans[d] = span;

        if (d == d_target) {
            evaluate_fast_cubic_basis_derivative(pt[d], span, d, dB.data());
        }
        else {
            evaluate_fast_cubic_basis(pt[d], span, d, B[d].data());
        }
    }

    // Step 2: offsets and strides
    int o0 = spans[0] - degree;
    int o1 = spans[1] - degree;
    int o2 = spans[2] - degree;
    int o3 = spans[3] - degree;
    int o4 = spans[4] - degree;
    int o5 = spans[5] - degree;

    const int s0 = control_point_strides_[0];
    const int s1 = control_point_strides_[1];
    const int s2 = control_point_strides_[2];
    const int s3 = control_point_strides_[3];
    const int s4 = control_point_strides_[4];
    const int s5 = control_point_strides_[5];

    // Step 3: tensor-product sum with substitution for derivative dimension
    REAL result = 0.0;

    for (int i = 0; i < span_size; ++i)
        for (int j = 0; j < span_size; ++j)
            for (int k = 0; k < span_size; ++k)
                for (int l = 0; l < span_size; ++l)
                    for (int m = 0; m < span_size; ++m)
                        for (int n = 0; n < span_size; ++n) {
                            REAL w = 1.0;
                            w *= (d_target == 0) ? dB[i] : B[0][i];
                            w *= (d_target == 1) ? dB[j] : B[1][j];
                            w *= (d_target == 2) ? dB[k] : B[2][k];
                            w *= (d_target == 3) ? dB[l] : B[3][l];
                            w *= (d_target == 4) ? dB[m] : B[4][m];
                            w *= (d_target == 5) ? dB[n] : B[5][n];

                            int idx = (o0 + i) * s0 + (o1 + j) * s1 + (o2 + k) * s2 +
                                (o3 + l) * s3 + (o4 + m) * s4 + (o5 + n) * s5;

                            result += w * coefficients_[idx];
                        }

    return result;
}

void Natural_BSpline_ND::evaluate_batch(int n_points, const REAL* input_coords, REAL* output_values) const {
    for (int i = 0; i < n_points; ++i) {
        const REAL* pt = input_coords + i * num_dims_;
        output_values[i] = evaluate(pt);  
    }
}

void Natural_BSpline_ND::evaluate_batch_derivatives(int n_points,
    const REAL* input_coords,
    int d_target,
    REAL* output_values) const
{
    for (int i = 0; i < n_points; ++i) {
        const REAL* pt = input_coords + i * num_dims_;
        output_values[i] = evaluate_derivative(pt, d_target);
    }
}
