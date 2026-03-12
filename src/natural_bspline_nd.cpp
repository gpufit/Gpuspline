
#include <numeric>
#include "natural_bspline_1d.h"
#include "natural_bspline_nd.h"
#include "libs/fitpack/gpuspline_fitpack_functions.h"
#include "bspline_fast_cubic_basis_evaluate.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#include <algorithm>
#endif


struct axis_lut_t {
    std::vector<int> o;               // span - degree, length = output_grid_dims[d]
    std::vector<REAL> b0, b1, b2, b3; // basis weights, length = output_grid_dims[d]
};


Natural_BSpline_ND::Natural_BSpline_ND(const std::vector<int>& data_dims, const REAL* data_values)
    : num_dims_(static_cast<int>(data_dims.size())), 
    data_dims_(data_dims) 
{

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
                                       const std::vector<REAL>& coefficients)
    : num_dims_(static_cast<int>(data_dims.size())),
    data_dims_(data_dims),
    coefficients_(coefficients)
{
    assert(num_dims_ >= 1);

    use_fast_evaluation_ = false;

    control_point_dims_.resize(num_dims_);
    int total_expected = 1;
    for (int i = 0; i < num_dims_; ++i) {
        int n_ctrl = data_dims[i] + 2;          // natural cubic spline # control points
        control_point_dims_[i] = n_ctrl;
        total_expected *= n_ctrl;
    }

    assert(static_cast<int>(coefficients.size()) == total_expected);

    total_data_points_ = std::accumulate(data_dims_.begin(), data_dims_.end(), 1, std::multiplies<int>());
    total_control_points_ = std::accumulate(control_point_dims_.begin(), control_point_dims_.end(), 1, std::multiplies<int>());

    // Initialize stride vectors
    init_stride_vectors();

    // Initialize multi-indices
    build_multi_indices();

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


static inline void evaluate_regular_grid_1d_lut(const Natural_BSpline_ND&, const std::vector<axis_lut_t>&, const int*, REAL*);
static inline void evaluate_regular_grid_2d_lut(const Natural_BSpline_ND&, const std::vector<axis_lut_t>&, const int*, REAL*);
static inline void evaluate_regular_grid_3d_lut(const Natural_BSpline_ND&, const std::vector<axis_lut_t>&, const int*, REAL*);
static inline void evaluate_regular_grid_4d_lut(const Natural_BSpline_ND&, const std::vector<axis_lut_t>&, const int*, REAL*);
static inline void evaluate_regular_grid_5d_lut(const Natural_BSpline_ND&, const std::vector<axis_lut_t>&, const int*, REAL*);
static inline void evaluate_regular_grid_6d_lut(const Natural_BSpline_ND&, const std::vector<axis_lut_t>&, const int*, REAL*);


void Natural_BSpline_ND::evaluate_regular_grid(
    const int* output_grid_dims,
    const REAL* output_origin,
    const REAL* output_step_size,
    REAL* output_values) const
{

#ifdef HAVE_OPENMP
    omp_set_dynamic(0);
    const int nproc = omp_get_num_procs();
    omp_set_num_threads(std::max(1, nproc - 1));
#endif

    assert(output_grid_dims && output_origin && output_step_size && output_values);

    for (int d = 0; d < num_dims_; ++d) {
        if (output_grid_dims[d] <= 0) return; // or throw / assert
    }

    auto fallback = [&]() {
        std::size_t total_pts = 1;
        for (int d = 0; d < num_dims_; ++d)
            total_pts *= static_cast<std::size_t>(output_grid_dims[d]);

        std::vector<int> idx(num_dims_, 0);
        std::vector<REAL> pt(num_dims_);
        for (int d = 0; d < num_dims_; ++d) pt[d] = output_origin[d];

        for (std::size_t lin = 0; lin < total_pts; ++lin) {
            output_values[lin] = evaluate(pt.data());

            int dd = 0;
            while (dd < num_dims_) {
                idx[dd] += 1;
                pt[dd] += output_step_size[dd];
                if (idx[dd] < output_grid_dims[dd]) break;
                idx[dd] = 0;
                pt[dd] = output_origin[dd];
                ++dd;
            }
        }
    };


    // Avoid pre-computing LUTs if not using fast evaluation
    if (!use_fast_evaluation_ || num_dims_ > 6) {
        fallback();
        return;
    }


    // Precompute LUTs (only when needed)
    std::vector<axis_lut_t> lut(num_dims_);

    for (int d = 0; d < num_dims_; ++d) {
        const int Nd = output_grid_dims[d];

        lut[d].o.resize(Nd);
        lut[d].b0.resize(Nd);
        lut[d].b1.resize(Nd);
        lut[d].b2.resize(Nd);
        lut[d].b3.resize(Nd);

        REAL x = output_origin[d];
        REAL B[4];

        for (int i = 0; i < Nd; ++i, x += output_step_size[d]) {
            const int span = find_knot_interval(x, d);
            lut[d].o[i] = span - SPLINE_DEGREE;

            evaluate_fast_cubic_basis(x, span, d, B);
            lut[d].b0[i] = B[0];
            lut[d].b1[i] = B[1];
            lut[d].b2[i] = B[2];
            lut[d].b3[i] = B[3];
        }
    }

    // Dispatch
    switch (num_dims_) {
        case 1: evaluate_regular_grid_1d_lut(*this, lut, output_grid_dims, output_values); break;
        case 2: evaluate_regular_grid_2d_lut(*this, lut, output_grid_dims, output_values); break;
        case 3: evaluate_regular_grid_3d_lut(*this, lut, output_grid_dims, output_values); break;
        case 4: evaluate_regular_grid_4d_lut(*this, lut, output_grid_dims, output_values); break;
        case 5: evaluate_regular_grid_5d_lut(*this, lut, output_grid_dims, output_values); break;
        case 6: evaluate_regular_grid_6d_lut(*this, lut, output_grid_dims, output_values); break;
        default: fallback(); break;
    }
}


static inline void evaluate_regular_grid_1d_lut(
    const Natural_BSpline_ND& spline,
    const std::vector<axis_lut_t>& lut,
    const int* out_dims,
    REAL* out)
{
    const auto& strides = spline.control_point_strides();
    const int s0 = strides[0];
    const REAL* coeff = spline.coefficients().data();

    const std::size_t n0 = (std::size_t)out_dims[0];

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ii0 = 0; ii0 < (int)n0; ++ii0) {
        const std::size_t i0 = (std::size_t)ii0;

        const int o0 = lut[0].o[i0];
        const REAL B0_0 = lut[0].b0[i0], B0_1 = lut[0].b1[i0], B0_2 = lut[0].b2[i0], B0_3 = lut[0].b3[i0];

        const int base0 = o0 * s0;

        REAL result = 0;
        result += B0_0 * coeff[base0 + 0 * s0];
        result += B0_1 * coeff[base0 + 1 * s0];
        result += B0_2 * coeff[base0 + 2 * s0];
        result += B0_3 * coeff[base0 + 3 * s0];

        out[i0] = result;
    }
}


static inline void evaluate_regular_grid_2d_lut(
    const Natural_BSpline_ND& spline,
    const std::vector<axis_lut_t>& lut,
    const int* out_dims,
    REAL* out)
{
    const auto& strides = spline.control_point_strides();
    const int s0 = strides[0], s1 = strides[1];
    const REAL* coeff = spline.coefficients().data();

    const std::size_t n0 = (std::size_t)out_dims[0];
    const std::size_t n1 = (std::size_t)out_dims[1];

    const std::size_t stride1 = n0;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ii1 = 0; ii1 < (int)n1; ++ii1) {
        const std::size_t i1 = (std::size_t)ii1;

        const int o1 = lut[1].o[i1];
        const REAL B1_0 = lut[1].b0[i1], B1_1 = lut[1].b1[i1], B1_2 = lut[1].b2[i1], B1_3 = lut[1].b3[i1];

        for (std::size_t i0 = 0; i0 < n0; ++i0) {
            const int o0 = lut[0].o[i0];
            const REAL B0_0 = lut[0].b0[i0], B0_1 = lut[0].b1[i0], B0_2 = lut[0].b2[i0], B0_3 = lut[0].b3[i0];

            const std::size_t lin = i0 + i1 * stride1;
            const int base0 = o0 * s0 + o1 * s1;

            REAL result = 0;

            // a,b stencil (16 terms); unroll b (dim1) inside a-loop
            for (int a = 0; a < 4; ++a) {
                const REAL Wa = (a == 0 ? B0_0 : (a == 1 ? B0_1 : (a == 2 ? B0_2 : B0_3)));
                const int ia = base0 + a * s0;

                result += Wa * B1_0 * coeff[ia + 0 * s1];
                result += Wa * B1_1 * coeff[ia + 1 * s1];
                result += Wa * B1_2 * coeff[ia + 2 * s1];
                result += Wa * B1_3 * coeff[ia + 3 * s1];
            }

            out[lin] = result;
        }
    }
}


static inline void evaluate_regular_grid_3d_lut(
    const Natural_BSpline_ND& spline,
    const std::vector<axis_lut_t>& lut,
    const int* out_dims,
    REAL* out)
{
    const auto& strides = spline.control_point_strides();
    const int s0 = strides[0], s1 = strides[1], s2 = strides[2];
    const REAL* coeff = spline.coefficients().data();

    const std::size_t n0 = (std::size_t)out_dims[0];
    const std::size_t n1 = (std::size_t)out_dims[1];
    const std::size_t n2 = (std::size_t)out_dims[2];

    const std::size_t stride0 = 1;
    const std::size_t stride1 = n0;
    const std::size_t stride2 = n0 * n1;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ii2 = 0; ii2 < (int)n2; ++ii2) {
        const std::size_t i2 = (std::size_t)ii2;

        const int o2 = lut[2].o[i2];
        const REAL B2_0 = lut[2].b0[i2], B2_1 = lut[2].b1[i2], B2_2 = lut[2].b2[i2], B2_3 = lut[2].b3[i2];

        for (std::size_t i1 = 0; i1 < n1; ++i1) {
            const int o1 = lut[1].o[i1];
            const REAL B1_0 = lut[1].b0[i1], B1_1 = lut[1].b1[i1], B1_2 = lut[1].b2[i1], B1_3 = lut[1].b3[i1];

            for (std::size_t i0 = 0; i0 < n0; ++i0) {
                const int o0 = lut[0].o[i0];
                const REAL B0_0 = lut[0].b0[i0], B0_1 = lut[0].b1[i0], B0_2 = lut[0].b2[i0], B0_3 = lut[0].b3[i0];

                const std::size_t lin = i0 * stride0 + i1 * stride1 + i2 * stride2;
                const int base0 = o0 * s0 + o1 * s1 + o2 * s2;

                REAL result = 0;

                for (int a = 0; a < 4; ++a) {
                    const REAL Wa = (a == 0 ? B0_0 : (a == 1 ? B0_1 : (a == 2 ? B0_2 : B0_3)));
                    const int ia = base0 + a * s0;

                    for (int b = 0; b < 4; ++b) {
                        const REAL Wab = Wa * (b == 0 ? B1_0 : (b == 1 ? B1_1 : (b == 2 ? B1_2 : B1_3)));
                        const int iab = ia + b * s1;

                        result += Wab * B2_0 * coeff[iab + 0 * s2];
                        result += Wab * B2_1 * coeff[iab + 1 * s2];
                        result += Wab * B2_2 * coeff[iab + 2 * s2];
                        result += Wab * B2_3 * coeff[iab + 3 * s2];
                    }
                }

                out[lin] = result;
            }
        }
    }
}


static inline void evaluate_regular_grid_4d_lut(
    const Natural_BSpline_ND& spline,
    const std::vector<axis_lut_t>& lut,
    const int* out_dims,
    REAL* out)
{
    const auto& strides = spline.control_point_strides();
    const int s0 = strides[0], s1 = strides[1], s2 = strides[2], s3 = strides[3];
    const REAL* coeff = spline.coefficients().data();

    const std::size_t n0 = (std::size_t)out_dims[0];
    const std::size_t n1 = (std::size_t)out_dims[1];
    const std::size_t n2 = (std::size_t)out_dims[2];
    const std::size_t n3 = (std::size_t)out_dims[3];

    const std::size_t stride0 = 1;
    const std::size_t stride1 = n0;
    const std::size_t stride2 = n0 * n1;
    const std::size_t stride3 = n0 * n1 * n2;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ii3 = 0; ii3 < (int)n3; ++ii3) {
        const std::size_t i3 = (std::size_t)ii3;

        const int o3 = lut[3].o[i3];
        const REAL B3_0 = lut[3].b0[i3], B3_1 = lut[3].b1[i3], B3_2 = lut[3].b2[i3], B3_3 = lut[3].b3[i3];

        for (std::size_t i2 = 0; i2 < n2; ++i2) {
            const int o2 = lut[2].o[i2];
            const REAL B2_0 = lut[2].b0[i2], B2_1 = lut[2].b1[i2], B2_2 = lut[2].b2[i2], B2_3 = lut[2].b3[i2];

            for (std::size_t i1 = 0; i1 < n1; ++i1) {
                const int o1 = lut[1].o[i1];
                const REAL B1_0 = lut[1].b0[i1], B1_1 = lut[1].b1[i1], B1_2 = lut[1].b2[i1], B1_3 = lut[1].b3[i1];

                for (std::size_t i0 = 0; i0 < n0; ++i0) {
                    const int o0 = lut[0].o[i0];
                    const REAL B0_0 = lut[0].b0[i0], B0_1 = lut[0].b1[i0], B0_2 = lut[0].b2[i0], B0_3 = lut[0].b3[i0];

                    const std::size_t lin =
                        i0 * stride0 + i1 * stride1 + i2 * stride2 + i3 * stride3;

                    const int base0 = o0 * s0 + o1 * s1 + o2 * s2 + o3 * s3;

                    REAL result = 0;

                    for (int a = 0; a < 4; ++a) {
                        const REAL Wa = (a == 0 ? B0_0 : (a == 1 ? B0_1 : (a == 2 ? B0_2 : B0_3)));
                        const int ia = base0 + a * s0;

                        for (int b = 0; b < 4; ++b) {
                            const REAL Wab = Wa * (b == 0 ? B1_0 : (b == 1 ? B1_1 : (b == 2 ? B1_2 : B1_3)));
                            const int iab = ia + b * s1;

                            for (int c = 0; c < 4; ++c) {
                                const REAL Wabc = Wab * (c == 0 ? B2_0 : (c == 1 ? B2_1 : (c == 2 ? B2_2 : B2_3)));
                                const int iabc = iab + c * s2;

                                result += Wabc * B3_0 * coeff[iabc + 0 * s3];
                                result += Wabc * B3_1 * coeff[iabc + 1 * s3];
                                result += Wabc * B3_2 * coeff[iabc + 2 * s3];
                                result += Wabc * B3_3 * coeff[iabc + 3 * s3];
                            }
                        }
                    }

                    out[lin] = result;
                }
            }
        }
    }
}


static inline void evaluate_regular_grid_5d_lut(
    const Natural_BSpline_ND& spline,
    const std::vector<axis_lut_t>& lut,
    const int* out_dims,
    REAL* out)
{
    const auto& strides = spline.control_point_strides();
    const int s0 = strides[0], s1 = strides[1], s2 = strides[2], s3 = strides[3], s4 = strides[4];

    const REAL* coeff = spline.coefficients().data();

    const std::size_t n0 = (std::size_t)out_dims[0];
    const std::size_t n1 = (std::size_t)out_dims[1];
    const std::size_t n2 = (std::size_t)out_dims[2];
    const std::size_t n3 = (std::size_t)out_dims[3];
    const std::size_t n4 = (std::size_t)out_dims[4];

    const std::size_t stride0 = 1;
    const std::size_t stride1 = n0;
    const std::size_t stride2 = n0 * n1;
    const std::size_t stride3 = n0 * n1 * n2;
    const std::size_t stride4 = n0 * n1 * n2 * n3;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ii4 = 0; ii4 < (int)n4; ++ii4) {
        const std::size_t i4 = (std::size_t)ii4;

        const int o4 = lut[4].o[i4];
        const REAL B4_0 = lut[4].b0[i4], B4_1 = lut[4].b1[i4], B4_2 = lut[4].b2[i4], B4_3 = lut[4].b3[i4];

        for (std::size_t i3 = 0; i3 < n3; ++i3) {
            const int o3 = lut[3].o[i3];
            const REAL B3_0 = lut[3].b0[i3], B3_1 = lut[3].b1[i3], B3_2 = lut[3].b2[i3], B3_3 = lut[3].b3[i3];

            for (std::size_t i2 = 0; i2 < n2; ++i2) {
                const int o2 = lut[2].o[i2];
                const REAL B2_0 = lut[2].b0[i2], B2_1 = lut[2].b1[i2], B2_2 = lut[2].b2[i2], B2_3 = lut[2].b3[i2];

                for (std::size_t i1 = 0; i1 < n1; ++i1) {
                    const int o1 = lut[1].o[i1];
                    const REAL B1_0 = lut[1].b0[i1], B1_1 = lut[1].b1[i1], B1_2 = lut[1].b2[i1], B1_3 = lut[1].b3[i1];

                    for (std::size_t i0 = 0; i0 < n0; ++i0) {
                        const int o0 = lut[0].o[i0];
                        const REAL B0_0 = lut[0].b0[i0], B0_1 = lut[0].b1[i0], B0_2 = lut[0].b2[i0], B0_3 = lut[0].b3[i0];

                        const std::size_t lin =
                            i0 * stride0 + i1 * stride1 + i2 * stride2 + i3 * stride3 + i4 * stride4;

                        const int base0 = o0 * s0 + o1 * s1 + o2 * s2 + o3 * s3 + o4 * s4;

                        REAL result = 0;

                        for (int a = 0; a < 4; ++a) {
                            const REAL Wa = (a == 0 ? B0_0 : (a == 1 ? B0_1 : (a == 2 ? B0_2 : B0_3)));
                            const int ia = base0 + a * s0;

                            for (int b = 0; b < 4; ++b) {
                                const REAL Wab = Wa * (b == 0 ? B1_0 : (b == 1 ? B1_1 : (b == 2 ? B1_2 : B1_3)));
                                const int iab = ia + b * s1;

                                for (int c = 0; c < 4; ++c) {
                                    const REAL Wabc = Wab * (c == 0 ? B2_0 : (c == 1 ? B2_1 : (c == 2 ? B2_2 : B2_3)));
                                    const int iabc = iab + c * s2;

                                    for (int d = 0; d < 4; ++d) {
                                        const REAL Wabcd = Wabc * (d == 0 ? B3_0 : (d == 1 ? B3_1 : (d == 2 ? B3_2 : B3_3)));
                                        const int iabcd = iabc + d * s3;

                                        result += Wabcd * B4_0 * coeff[iabcd + 0 * s4];
                                        result += Wabcd * B4_1 * coeff[iabcd + 1 * s4];
                                        result += Wabcd * B4_2 * coeff[iabcd + 2 * s4];
                                        result += Wabcd * B4_3 * coeff[iabcd + 3 * s4];
                                    }
                                }
                            }
                        }

                        out[lin] = result;
                    }
                }
            }
        }
    }
}


static inline void evaluate_regular_grid_6d_lut(
    const Natural_BSpline_ND& spline,
    const std::vector<axis_lut_t>& lut,
    const int* out_dims,
    REAL* out)
{
    const auto& strides = spline.control_point_strides();
    const int s0 = strides[0], s1 = strides[1], s2 = strides[2], s3 = strides[3], s4 = strides[4], s5 = strides[5];
    const REAL* coeff = spline.coefficients().data();

    const std::size_t n0 = (std::size_t)out_dims[0];
    const std::size_t n1 = (std::size_t)out_dims[1];
    const std::size_t n2 = (std::size_t)out_dims[2];
    const std::size_t n3 = (std::size_t)out_dims[3];
    const std::size_t n4 = (std::size_t)out_dims[4];
    const std::size_t n5 = (std::size_t)out_dims[5];

    const std::size_t stride0 = 1;
    const std::size_t stride1 = n0;
    const std::size_t stride2 = n0 * n1;
    const std::size_t stride3 = n0 * n1 * n2;
    const std::size_t stride4 = n0 * n1 * n2 * n3;
    const std::size_t stride5 = n0 * n1 * n2 * n3 * n4;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ii5 = 0; ii5 < (int)n5; ++ii5) {
        const std::size_t i5 = (std::size_t)ii5;

        const int o5 = lut[5].o[i5];
        const REAL B5_0 = lut[5].b0[i5], B5_1 = lut[5].b1[i5], B5_2 = lut[5].b2[i5], B5_3 = lut[5].b3[i5];

        for (std::size_t i4 = 0; i4 < n4; ++i4) {
            const int o4 = lut[4].o[i4];
            const REAL B4_0 = lut[4].b0[i4], B4_1 = lut[4].b1[i4], B4_2 = lut[4].b2[i4], B4_3 = lut[4].b3[i4];

            for (std::size_t i3 = 0; i3 < n3; ++i3) {
                const int o3 = lut[3].o[i3];
                const REAL B3_0 = lut[3].b0[i3], B3_1 = lut[3].b1[i3], B3_2 = lut[3].b2[i3], B3_3 = lut[3].b3[i3];

                for (std::size_t i2 = 0; i2 < n2; ++i2) {
                    const int o2 = lut[2].o[i2];
                    const REAL B2_0 = lut[2].b0[i2], B2_1 = lut[2].b1[i2], B2_2 = lut[2].b2[i2], B2_3 = lut[2].b3[i2];

                    for (std::size_t i1 = 0; i1 < n1; ++i1) {
                        const int o1 = lut[1].o[i1];
                        const REAL B1_0 = lut[1].b0[i1], B1_1 = lut[1].b1[i1], B1_2 = lut[1].b2[i1], B1_3 = lut[1].b3[i1];

                        for (std::size_t i0 = 0; i0 < n0; ++i0) {
                            const int o0 = lut[0].o[i0];
                            const REAL B0_0 = lut[0].b0[i0], B0_1 = lut[0].b1[i0], B0_2 = lut[0].b2[i0], B0_3 = lut[0].b3[i0];

                            const std::size_t lin =
                                i0 * stride0 + i1 * stride1 + i2 * stride2 + i3 * stride3 + i4 * stride4 + i5 * stride5;

                            const int base0 = o0 * s0 + o1 * s1 + o2 * s2 + o3 * s3 + o4 * s4 + o5 * s5;

                            REAL result = 0;

                            for (int a = 0; a < 4; ++a) {
                                const REAL Wa = (a == 0 ? B0_0 : (a == 1 ? B0_1 : (a == 2 ? B0_2 : B0_3)));
                                const int ia = base0 + a * s0;

                                for (int b = 0; b < 4; ++b) {
                                    const REAL Wab = Wa * (b == 0 ? B1_0 : (b == 1 ? B1_1 : (b == 2 ? B1_2 : B1_3)));
                                    const int iab = ia + b * s1;

                                    for (int c = 0; c < 4; ++c) {
                                        const REAL Wabc = Wab * (c == 0 ? B2_0 : (c == 1 ? B2_1 : (c == 2 ? B2_2 : B2_3)));
                                        const int iabc = iab + c * s2;

                                        for (int d = 0; d < 4; ++d) {
                                            const REAL Wabcd = Wabc * (d == 0 ? B3_0 : (d == 1 ? B3_1 : (d == 2 ? B3_2 : B3_3)));
                                            const int iabcd = iabc + d * s3;

                                            for (int e = 0; e < 4; ++e) {
                                                const REAL Wabcde = Wabcd * (e == 0 ? B4_0 : (e == 1 ? B4_1 : (e == 2 ? B4_2 : B4_3)));
                                                const int iabcde = iabcd + e * s4;

                                                result += Wabcde * B5_0 * coeff[iabcde + 0 * s5];
                                                result += Wabcde * B5_1 * coeff[iabcde + 1 * s5];
                                                result += Wabcde * B5_2 * coeff[iabcde + 2 * s5];
                                                result += Wabcde * B5_3 * coeff[iabcde + 3 * s5];
                                            }
                                        }
                                    }
                                }
                            }

                            out[lin] = result;
                        }
                    }
                }
            }
        }
    }
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
