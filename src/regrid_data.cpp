#include "spline_classes.h"
#include <algorithm>


// Gaussian weighting function
REAL gaussian_weight(REAL w, REAL v, REAL w_i, REAL v_i, REAL sigma) {
    REAL dist2 = (w - w_i) * (w - w_i) + (v - v_i) * (v - v_i);
    return std::exp(-dist2 / (2.0 * sigma * sigma));
}

// Function to compute the median distance to the K nearest neighbors
REAL compute_adaptive_sigma(const REAL* w_samples, const REAL* v_samples, int num_samples, REAL w_target, REAL v_target, int K) {
    std::vector<REAL> distances;
    for (int i = 0; i < num_samples; i++) {
        REAL dist2 = (w_samples[i] - w_target) * (w_samples[i] - w_target) +
            (v_samples[i] - v_target) * (v_samples[i] - v_target);
        distances.push_back(std::sqrt(dist2));
    }

    std::sort(distances.begin(), distances.end());
    return distances[K / 2];
}


// Function to interpolate 3D volumes onto a regular (w, v) grid
void regrid_3d_volumes(
    // Inputs
    int x_size,
    int y_size,
    int z_size,
    int num_samples,
    const REAL* volumes,
    const REAL* w_samples,
    const REAL* v_samples,
    int w_grid_size,
    int v_grid_size,
    int K,
    REAL c,
    const REAL* w_grid,
    const REAL* v_grid,
    // Outputs
    REAL* output_volumes,
    int* output_vol_n_avg,
    int* output_vol_avg_indices,
    REAL* output_vol_avg_weights)
{
    for (int wi = 0; wi < w_grid_size; wi++) {
        for (int vi = 0; vi < v_grid_size; vi++) {
            REAL w_target = w_grid[wi];
            REAL v_target = v_grid[vi];

            // Compute initial adaptive sigma
            REAL sigma = c * compute_adaptive_sigma(w_samples, v_samples, num_samples, w_target, v_target, K);

            // Initialize weight sums, output volume, and count of used volumes
            REAL weight_sum = 0.0;
            int count = 0;
            int index_offset = (wi * v_grid_size + vi) * num_samples;

            for (int z = 0; z < z_size; z++) {
                for (int y = 0; y < y_size; y++) {
                    for (int x = 0; x < x_size; x++) {
                        output_volumes[x + y * x_size + z * x_size * y_size + wi * x_size * y_size * z_size + vi * x_size * y_size * z_size * w_grid_size] = 0.0;
                    }
                }
            }

            // Perform Gaussian-weighted averaging with adaptive sigma widening
            int iteration = 0;
            int max_iterations = 5; // Limit sigma widening attempts
            while (weight_sum == 0.0 && iteration < max_iterations) {
                for (int i = 0; i < num_samples; i++) {
                    REAL weight = gaussian_weight(w_target, v_target, w_samples[i], v_samples[i], sigma);
                    if (weight > 0.0) {
                        output_vol_avg_indices[index_offset + count] = i;
                        output_vol_avg_weights[index_offset + count] = weight;
                        count++;
                    }

                    for (int z = 0; z < z_size; z++) {
                        for (int y = 0; y < y_size; y++) {
                            for (int x = 0; x < x_size; x++) {
                                int input_idx = x + y * x_size + z * x_size * y_size + i * x_size * y_size * z_size;
                                int output_idx = x + y * x_size + z * x_size * y_size + wi * x_size * y_size * z_size + vi * x_size * y_size * z_size * w_grid_size;

                                output_volumes[output_idx] += weight * volumes[input_idx];
                            }
                        }
                    }
                    weight_sum += weight;
                }

                // If no weight contribution, increase sigma and retry
                if (weight_sum == 0.0) {
                    sigma *= 1.5; // Increase sigma by 50%
                }
                iteration++;
            }

            // Normalize output volume
            if (weight_sum > 0.0) {
                for (int z = 0; z < z_size; z++) {
                    for (int y = 0; y < y_size; y++) {
                        for (int x = 0; x < x_size; x++) {
                            int output_idx = x + y * x_size + z * x_size * y_size + wi * x_size * y_size * z_size + vi * x_size * y_size * z_size * w_grid_size;
                            output_volumes[output_idx] /= weight_sum;
                        }
                    }
                }
            }

            // Store the number of volumes used for this grid point
            output_vol_n_avg[wi * v_grid_size + vi] = count;
        }
    }
}



int regrid_3d_volumes_portable(int argc, void *argv[])
{
    regrid_3d_volumes(
        *((int *)argv[0]),
        *((int *)argv[1]),
        *((int *)argv[2]),
        *((int *)argv[3]),
        (REAL *)argv[4],
        (REAL *)argv[5],
        (REAL *)argv[6],
        *((int *)argv[7]),
        *((int *)argv[8]),
        *((int *)argv[9]),
        *((REAL *)argv[10]),
        (REAL *)argv[11],
        (REAL *)argv[12],
        (REAL *)argv[13], 
        (int *)argv[14],
        (int *)argv[15],
        (REAL *)argv[16]);

    return 0;
}