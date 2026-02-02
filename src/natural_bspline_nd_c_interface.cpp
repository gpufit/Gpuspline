#include "spline.h"
#include "natural_bspline_nd.h"
#include <algorithm>

int calculate_coefficients_natural_bspline_nd(
    int num_dims,
    const int* dims,
    REAL* data,
    REAL* knots_out,
    REAL* coeff_strides_out,
    REAL* coefficients_out)
{
    // Convert dims to std::vector<int>
    std::vector<int> shape(num_dims);
    for (int i = 0; i < num_dims; ++i)
        shape[i] = dims[i];

    // Create spline and compute coefficients
    Natural_BSpline_ND spline(shape, data);

    // Copy coefficients
    const std::vector<REAL>& coeffs = spline.coefficients();
    std::copy(coeffs.begin(), coeffs.end(), coefficients_out);

    // ---- coeff strides (axis0 fastest) ----
    if (coeff_strides_out)
    {
        const std::vector<int>& strides = spline.control_point_strides();
        for (int d = 0; d < num_dims; ++d)
            coeff_strides_out[d] = static_cast<REAL>(strides[d]); // stored as REAL for ABI simplicity
    }

    // ---- knots (concatenate axis0, axis1, ...) ----
    if (knots_out)
    {
        const auto& kv = spline.knot_vectors();

        std::size_t off = 0;
        for (int d = 0; d < num_dims; ++d)
        {
            const std::vector<REAL>& k = kv[d];
            std::copy(k.begin(), k.end(), knots_out + off);
            off += k.size();
        }
    }

    return 0;
}


int calculate_coefficients_natural_bspline_nd_portable(int argc, void* argv[])
{
    return calculate_coefficients_natural_bspline_nd(
        *((int*) argv[0]),
        (int*) argv[1],
        (REAL*) argv[2],
        (REAL*) argv[3],
        (REAL*) argv[4],
        (REAL*) argv[5]);
}


int calculate_values_natural_bspline_nd(
    int num_dims,
    const int* data_dims,
    const REAL* coefficients,
    int n_input_coords,
    const REAL* input_coords,
    int flag_fast_evaluate, 
    REAL* output_values)
{
    using ND = Natural_BSpline_ND;

    // Convert data_dims to vector<int>
    std::vector<int> shape(num_dims);
    for (int i = 0; i < num_dims; ++i)
        shape[i] = data_dims[i];

    // Reconstruct coefficients (total = product(data_dims[d] + 2))
    std::size_t total_coeffs = 1;
    for (int i = 0; i < num_dims; ++i)
        total_coeffs *= static_cast<std::size_t>(data_dims[i] + 2);
    std::vector<REAL> coeff_vec(coefficients, coefficients + total_coeffs);

    // Create spline
    ND spline(shape, coeff_vec);

    // Set fast evaluation
    spline.set_fast_evaluation(flag_fast_evaluate != 0);

    // Evaluate
    spline.evaluate_batch(n_input_coords, input_coords, output_values);

    return 0;
}


int calculate_values_natural_bspline_nd_portable(int argc, void* argv[])
{
    return calculate_values_natural_bspline_nd(
        *((int*) argv[0]),
        (int*) argv[1],
        (REAL*) argv[2],
        *((int*) argv[3]),
        (REAL*) argv[4],
        *((int*)argv[5]),
        (REAL*) argv[6]);
}


int calculate_sizes_natural_bspline_nd(
    int num_dims,
    const int* dims,
    size_t* knots_len_out,
    size_t* coeff_strides_len_out,
    size_t* coefficients_len_out)
{
    if (num_dims <= 0 || !dims)
        return 1;

    size_t knots_len = 0;
    size_t coeff_len = 1;

    for (int d = 0; d < num_dims; ++d)
    {
        int const N = dims[d];
        if (N <= 0)
            return 2;

        knots_len += (size_t)(N + 6);

        // control points per axis = N + 2
        coeff_len *= (size_t)(N + 2);

    }

    if (knots_len_out) *knots_len_out = (size_t)knots_len;
    if (coeff_strides_len_out) *coeff_strides_len_out = num_dims;
    if (coefficients_len_out) *coefficients_len_out = (size_t)coeff_len;

    return 0;
}


int calculate_sizes_natural_bspline_nd_portable(int argc, void* argv[])
{
    return calculate_sizes_natural_bspline_nd(
        *((int*)argv[0]),
        (int*)argv[1],
        (size_t*)argv[2],
        (size_t*)argv[3],
        (size_t*)argv[4]);
}