#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "natural_bspline_1d.h"
#include "natural_bspline_nd.h"

constexpr int N_1D = 11;
constexpr int NX_2D = 11;
constexpr int NY_2D = 9;

// === 1D Utilities ===
void generate_1d_data(std::vector<REAL>& data) {
    data.resize(N_1D);
    for (int i = 0; i < N_1D; ++i)
        data[i] = std::sin(2.0 * M_PI * i / (N_1D - 1));
}

void export_1d_ground_truth(const std::vector<REAL>& data, const std::string& filename) {
    std::ofstream out(filename);
    out << "x,f_true\n";
    for (int i = 0; i < N_1D; ++i)
        out << i << "," << data[i] << "\n";
}

void export_1d_spline_eval(Natural_BSpline_1D& spline, const std::string& filename) {
    std::ofstream out(filename);
    out << "x,f_slow,df_slow,f_fast,df_fast\n";

    for (REAL x = 0.0; x <= N_1D - 1; x += 0.1f) {
        spline.set_fast_evaluation(false);
        REAL f_slow = spline.evaluate(x);
        REAL df_slow = spline.evaluate_derivative(x);

        spline.set_fast_evaluation(true);
        REAL f_fast = spline.evaluate(x);
        REAL df_fast = spline.evaluate_derivative(x);

        out << x << "," << f_slow << "," << df_slow << "," << f_fast << "," << df_fast << "\n";
    }
}

void export_nd_vs_1d(Natural_BSpline_1D& spline_1d, const Natural_BSpline_ND& spline_nd, const std::string& filename) {
    std::ofstream out(filename);
    out << "x,f_1d,df_1d,f_nd,df_nd\n";
    for (REAL x = 0.0; x <= N_1D - 1; x += 0.1f) {
        spline_1d.set_fast_evaluation(false);
        REAL f1 = spline_1d.evaluate(x);
        REAL d1 = spline_1d.evaluate_derivative(x);

        std::vector<REAL> x_nd = { x };
        REAL f2 = spline_nd.evaluate(x_nd.data());
        REAL d2 = spline_nd.evaluate_derivative(x_nd.data(), 0);

        out << x << "," << f1 << "," << d1 << "," << f2 << "," << d2 << "\n";
    }
}

// === 2D Utilities ===
void generate_2d_data(std::vector<REAL>& data) {
    data.resize(NX_2D * NY_2D);
    for (int j = 0; j < NY_2D; ++j) {
        for (int i = 0; i < NX_2D; ++i) {
            data[j * NX_2D + i] = std::sin(2.0 * M_PI * i / (NX_2D - 1)) * std::cos(2.0 * M_PI * j / (NY_2D - 1));
        }
    }
}

void export_2d_slice_x_axis(const Natural_BSpline_ND& spline, const std::string& filename) {
    std::ofstream out(filename);
    out << "x,y_fixed,f_val,df_dx,df_dy\n";
    REAL y = 4.0; // mid-slice along y

    for (REAL x = 0.0; x <= NX_2D - 1; x += 0.1f) {
        std::vector<REAL> pt = { x, y };
        REAL f = spline.evaluate(pt.data());
        REAL dfdx = spline.evaluate_derivative(pt.data(), 0);
        REAL dfdy = spline.evaluate_derivative(pt.data(), 1);
        out << x << "," << y << "," << f << "," << dfdx << "," << dfdy << "\n";
    }
}

void export_2d_slice_y_axis(const Natural_BSpline_ND& spline, const std::string& filename) {
    std::ofstream out(filename);
    out << "y,x_fixed,f_val,df_dx,df_dy\n";
    REAL x = 5.0; // mid-slice along x

    for (REAL y = 0.0; y <= NY_2D - 1; y += 0.1f) {
        std::vector<REAL> pt = { x, y };
        REAL f = spline.evaluate(pt.data());
        REAL dfdx = spline.evaluate_derivative(pt.data(), 0);
        REAL dfdy = spline.evaluate_derivative(pt.data(), 1);
        out << y << "," << x << "," << f << "," << dfdx << "," << dfdy << "\n";
    }
}

void run_small_2d_sanity_test(bool fast_eval) {
    std::cout << "\n=== Small 3x3 Sanity Check (" << (fast_eval ? "fast" : "slow") << " evaluation) ===\n";

    int Nx = 3, Ny = 3;
    std::vector<int> shape2d = { Nx, Ny };
    std::vector<REAL> data(Nx * Ny);

    // Fill grid: f(x, y) = x + 10*y
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            data[iy * Nx + ix] = REAL(ix + 10 * iy);
        }
    }

    // Build spline
    Natural_BSpline_ND spline(shape2d, data.data());

    // Configure 1D evaluation mode for internal solves
    Natural_BSpline_1D dummy(Nx);  // reuse as setter
    dummy.set_fast_evaluation(fast_eval);

    // Query points
    std::vector<REAL> x_vals = { 0.0, 1.0, 2.0 };
    std::vector<REAL> y_vals = { 0.0, 1.0, 2.0 };

    for (REAL y : y_vals) {
        for (REAL x : x_vals) {
            std::vector<REAL> pt = { x, y };
            REAL val = spline.evaluate(pt.data());
            std::cout << "spline(" << x << ", " << y << ") = " << val << "\n";
        }
    }
}


void run_rectangular_2d_sanity_test() {
    std::cout << "\n=== Rectangular 2D Sanity Test (left-to-right gradient) ===\n";

    int Nx = 5, Ny = 3;
    std::vector<int> shape2d = { Nx, Ny };
    std::vector<REAL> data(Nx * Ny);

    // Fill with horizontal gradient: f(x, y) = x
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
            data[j * Nx + i] = static_cast<REAL>(i);  // x-coordinate only

    // Create spline
    Natural_BSpline_ND spline(shape2d, data.data());

    // Evaluate at exact grid points
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            std::vector<REAL> pt = { static_cast<REAL>(i), static_cast<REAL>(j) };
            REAL val = spline.evaluate(pt.data());
            std::cout << "f(" << i << ", " << j << ") = " << val << "\n";
        }
    }
}


void test_reconstruction_from_coefficients() {
    std::cout << "\n=== Testing ND Spline Reconstruction ===\n";

    std::vector<int> shape = { NX_2D, NY_2D };
    std::vector<REAL> data;
    generate_2d_data(data);

    // Original spline
    Natural_BSpline_ND spline_orig(shape, data.data());

    const std::vector<REAL>& coeffs = spline_orig.coefficients();
    const auto& multi = spline_orig.multi_indices();

    // Reconstructed spline
    Natural_BSpline_ND spline_recon(shape, coeffs, multi);

    // Compare evaluation at several points
    for (REAL y = 0.0; y <= NY_2D - 1; y += 0.5f) {
        for (REAL x = 0.0; x <= NX_2D - 1; x += 0.5f) {
            std::vector<REAL> pt = { x, y };
            REAL f1 = spline_orig.evaluate(pt.data());
            REAL f2 = spline_recon.evaluate(pt.data());
            if (std::abs(f1 - f2) > 1e-5f)
                std::cerr << "Mismatch at (" << x << "," << y << "): " << f1 << " vs " << f2 << "\n";
        }
    }

    std::cout << "Reconstruction test complete.\n";
}



// === Main ===
int main() {
    std::cout << "=== Natural B-Spline Test ===\n";

    // 1D Test
    std::vector<REAL> data_1d;
    generate_1d_data(data_1d);
    Natural_BSpline_1D spline_1d(N_1D, data_1d.data());

    std::vector<int> shape_1d = { N_1D };
    Natural_BSpline_ND spline_1d_nd(shape_1d, data_1d.data());

    export_1d_ground_truth(data_1d, "cpp_ground_truth.csv");
    export_1d_spline_eval(spline_1d, "cpp_natural_spline.csv");
    export_nd_vs_1d(spline_1d, spline_1d_nd, "cpp_nd_vs_1d_comparison.csv");

    std::cout << "1D spline tests written.\n";

    // 2D Test
    std::vector<REAL> data_2d;
    generate_2d_data(data_2d);
    std::vector<int> shape_2d = { NX_2D, NY_2D };
    Natural_BSpline_ND spline_2d(shape_2d, data_2d.data());
    spline_2d.set_fast_evaluation(true);

    export_2d_slice_x_axis(spline_2d, "cpp_2d_slice_x.csv");
    export_2d_slice_y_axis(spline_2d, "cpp_2d_slice_y.csv");

    std::cout << "2D spline slice tests written.\n";


    run_small_2d_sanity_test(false);  // Slow basis evaluation
    run_small_2d_sanity_test(true);   // Fast basis evaluation


    test_reconstruction_from_coefficients();
    run_rectangular_2d_sanity_test();

    return 0;
}

 