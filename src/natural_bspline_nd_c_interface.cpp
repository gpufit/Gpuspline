#include "natural_bspline_nd.h"

int calculate_coefficients_bspline_nd(
    int num_dims,
    const int* dims,
    REAL* data,
    int* multi_indices_out,  // flattened, shape = (4^D) × D
    REAL* coefficients_out
)
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

    // Copy multi-indices
    const auto& multi_indices = spline.multi_indices();  
    for (int i = 0; i < multi_indices.size(); ++i) {
        for (int d = 0; d < num_dims; ++d) {
            multi_indices_out[i * num_dims + d] = multi_indices[i][d];
        }
    }

    return 0;
}


int calculate_coefficients_bspline_nd_portable(int argc, void* argv[])
{
    return calculate_coefficients_bspline_nd(
        *((int*) argv[0]),
        (int*) argv[1],
        (REAL*) argv[2],
        (int*) argv[3],
        (REAL*) argv[4]);
}


int calculate_values_bspline_nd(
    int num_dims,
    const int* data_dims,
    const REAL* coefficients,
    const int* multi_indices,
    int n_input_coords,
    const REAL* input_coords,
    REAL* output_values)
{
    using ND = Natural_BSpline_ND;

    // Convert data_dims to vector<int>
    std::vector<int> shape(num_dims);
    for (int i = 0; i < num_dims; ++i)
        shape[i] = data_dims[i];

    // Reconstruct multi_indices
    int num_basis = static_cast<int>(std::pow(ND::SPLINE_DEGREE + 1, num_dims));
    std::vector<std::array<int, ND::MAX_NDIMS>> multi_vec(num_basis);
    for (int i = 0; i < num_basis; ++i)
        for (int d = 0; d < num_dims; ++d)
            multi_vec[i][d] = multi_indices[i * num_dims + d];

    // Reconstruct coefficients (total = product(data_dims[d] + 2))
    std::size_t total_coeffs = 1;
    for (int i = 0; i < num_dims; ++i)
        total_coeffs *= static_cast<std::size_t>(data_dims[i] + 2);
    std::vector<REAL> coeff_vec(coefficients, coefficients + total_coeffs);

    // Create spline
    ND spline(shape, coeff_vec, multi_vec);

    // Prepare input points
    std::vector<std::vector<REAL>> points(n_input_coords, std::vector<REAL>(num_dims));
    for (int i = 0; i < n_input_coords; ++i)
        for (int d = 0; d < num_dims; ++d)
            points[i][d] = input_coords[i * num_dims + d];

    // Evaluate
    std::vector<REAL> result;
    spline.evaluate_batch(points, result);

    // Copy to output
    std::copy(result.begin(), result.end(), output_values);

    return 0;
}


int calculate_values_bspline_nd_portable(int argc, void* argv[])
{
    return calculate_values_bspline_nd(
        *((int*) argv[0]),
        (int*) argv[1],
        (REAL*) argv[2],
        (int*) argv[3],
        *((int*) argv[4]),
        (REAL*) argv[5],
        (REAL*) argv[6]);
}
