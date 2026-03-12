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
    const int* spline_domain_dims,
    const REAL* spline_coefficients,
    int n_input_coords,
    const REAL* input_coords,       // coords array dimensions: [num_dims, n_input_coords]
    int flag_fast_evaluate, 
    REAL* output_values)
{
    using ND = Natural_BSpline_ND;

    // Convert spline_domain_dims to vector<int>
    std::vector<int> shape(num_dims);
    for (int i = 0; i < num_dims; ++i)
        shape[i] = spline_domain_dims[i];

    // Reconstruct spline_coefficients (total = product(spline_domain_dims[d] + 2))
    std::size_t total_coeffs = 1;
    for (int i = 0; i < num_dims; ++i)
        total_coeffs *= static_cast<std::size_t>(spline_domain_dims[i] + 2);
    std::vector<REAL> coeff_vec(spline_coefficients, spline_coefficients + total_coeffs);

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


int calculate_values_regular_grid_natural_bspline_nd(
    int num_dims,
    const int* spline_domain_dims,      // defines coefficient shape: (dims[d] + 2)
    const REAL* spline_coefficients,
    const int* output_grid_dims,        // number of output samples per dim
    const REAL* output_origin,          // length num_dims
    const REAL* output_step_size,       // length num_dims
    int flag_fast_evaluate,
    REAL* output_values)                // length = product(output_grid_dims[d]), x fastest
{
    using ND = Natural_BSpline_ND;

    if (num_dims <= 0 ||
        !spline_domain_dims || !spline_coefficients ||
        !output_grid_dims || !output_origin || !output_step_size ||
        !output_values)
        return 1;

    // spline shape
    std::vector<int> shape(num_dims);
    for (int d = 0; d < num_dims; ++d)
        shape[d] = spline_domain_dims[d];

    // coeff count = product(spline_domain_dims[d] + 2)
    std::size_t total_coeffs = 1;
    for (int d = 0; d < num_dims; ++d)
        total_coeffs *= static_cast<std::size_t>(spline_domain_dims[d] + 2);

    std::vector<REAL> coeff_vec(spline_coefficients, spline_coefficients + total_coeffs);

    ND spline(shape, coeff_vec);
    spline.set_fast_evaluation(flag_fast_evaluate != 0);

    spline.evaluate_regular_grid(output_grid_dims, output_origin, output_step_size, output_values);

    return 0;
}


int calculate_values_regular_grid_natural_bspline_nd_portable(int argc, void* argv[])
{
    return calculate_values_regular_grid_natural_bspline_nd(
        *((int*)argv[0]),
        (int*)argv[1],
        (REAL*)argv[2],
        (int*)argv[3],
        (REAL*)argv[4],
        (REAL*)argv[5],
        *((int*)argv[6]),
        (REAL*)argv[7]);
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