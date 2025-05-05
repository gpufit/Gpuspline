#include <iostream>
#include <fstream>
#include <iomanip>
#include "spline_classes.h"


void print_solver_debug_info(const Natural_BSpline_1D& spline) {
    const auto& A = spline.get_collocation_matrix();
    const auto& b = spline.get_rhs_vector();
    int M = spline.num_control_points();

    std::cout << "\n=== Collocation Matrix A (" << M << " x " << M << ") ===\n";
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            std::cout << std::fixed << std::setw(10) << std::setprecision(6) << A(i, j) << " ";
        }
        std::cout << '\n';
    }

    std::cout << "\n=== RHS Vector b (" << M << ") ===\n";
    for (int i = 0; i < M; ++i) {
        std::cout << "  b[" << std::setw(2) << i << "] = " << std::fixed << std::setprecision(6) << b(i) << "\n";
    }

    std::cout << "\n=== Spline Coefficients ===\n";
    const REAL* coeffs = spline.coefficients();
    for (int i = 0; i < M; ++i) {
        std::cout << "  coeff[" << std::setw(2) << i << "] = " << std::fixed << std::setprecision(6) << coeffs[i] << "\n";
    }
    std::cout << std::endl;
}


void print_basis_curve_to_csv(const Natural_BSpline_1D& spline, int input_span, REAL x0, const std::string& filename) {
    std::ofstream out(filename);
    out << "x,B0,B1,B2,B3\n";

    constexpr int NUM_POINTS = 100;
    constexpr int p = 3;
    int span = input_span;

    for (int i = 0; i <= NUM_POINTS; ++i) {
        REAL x = x0 + static_cast<REAL>(i) / NUM_POINTS;  // x in [0, 1]
        REAL basis[p + 1] = {};
        spline.evaluate_basis(x, span, 0, basis);

        out << x - x0;
        for (int j = 0; j <= p; ++j) {
            out << "," << basis[j];
        }
        out << "\n";
    }

    out.close();
    std::cout << "Basis curve data saved to " << filename << "\n";
}


int main() {

    // Create the test data
    const int N = 11;
    REAL data[N];
    for (int i = 0; i < N; ++i) {
        data[i] = std::sin(2.0f * M_PI * i / REAL((N - 1)));
    }

    // Construct spline
    Natural_BSpline_1D spline_1d(N, data);

    // Set up ND spline (1D case)
    std::vector<int> shape = { N };
    std::vector<REAL> data_vec(data, data + N);
    Natural_BSpline_ND spline_nd(shape, data_vec.data());

    // Output spline interpolation results
    std::ofstream spline_out("cpp_natural_spline.csv");
    spline_out << "x,f_slow,df_slow,f_fast,df_fast\n";

    // Output raw data points (ground truth)
    std::ofstream truth_out("cpp_ground_truth.csv");
    truth_out << "x,f_true\n";
    for (int i = 0; i < N; ++i) {
        truth_out << i << "," << data[i] << "\n";
    }

    // Sample spline values at fine-grained x positions
    for (REAL x = 0.0; x <= REAL(N - 1); x += 0.1f) {
        spline_1d.set_fast_evaluation(false);
        REAL f_slow = spline_1d.evaluate(x);
        REAL df_slow = spline_1d.evaluate_derivative(x);

        spline_1d.set_fast_evaluation(true);
        REAL f_fast = spline_1d.evaluate(x);
        REAL df_fast = spline_1d.evaluate_derivative(x);

        spline_out << x << "," << f_slow << "," << df_slow << ","
            << f_fast << "," << df_fast << "\n";
    }

    spline_out.close();
    truth_out.close();

    //print_solver_debug_info(spline);

    std::cout << "=== C++ Knot Vector ===\n";
    for (int i = 0; i < spline_1d.num_knots(); ++i) {
        std::cout << "t[" << i << "] = " << spline_1d.knot(i) << "\n";
    }

    std::cout << "=== C++ Coefficients ===\n";
    for (int i = 0; i < spline_1d.num_control_points(); ++i) {
        std::cout << "c[" << i << "] = " << spline_1d.coefficients()[i] << "\n";
    }

    std::cout << "Output written to cpp_natural_spline.csv and cpp_ground_truth.csv\n";

    //std::cout << "\n=== Comparing slow vs. fast basis function values ===\n";
    //for (REAL x = 0.0; x <= static_cast<REAL>(spline.num_data_points() - 1); x += 0.25) {
    //    std::cout << "x = " << x << "\n";
    //    spline.debug_compare_basis(x);
    //}

    std::ofstream nd_out("cpp_nd_vs_1d_comparison.csv");
    nd_out << "x,f_1d,df_1d,f_nd,df_nd\n";

    std::cout << "cpp_nd_vs_1d_comparison";
    std::cout << "x,f_1d,df_1d,f_nd,df_nd\n";

    for (REAL x = 0.0; x <= REAL(N - 1); x += 0.1f) {
        spline_1d.set_fast_evaluation(false);
        REAL f_1d = spline_1d.evaluate(x);
        REAL df_1d = spline_1d.evaluate_derivative(x);

        std::vector<REAL> input = { x };
        REAL f_nd = spline_nd.evaluate(input);
        REAL df_nd = spline_nd.evaluate_derivative(input, 0);  // derivative in x

        nd_out << x << "," << f_1d << "," << df_1d << "," << f_nd << "," << df_nd << "\n";
        std::cout << x << "," << f_1d << "," << df_1d << "," << f_nd << "," << df_nd << "\n";
    }

    nd_out.close();
    std::cout << "ND vs 1D comparison written to cpp_nd_vs_1d_comparison.csv\n";


    //print_basis_curve_to_csv(spline_1d, 3, 0.0, "basis_span3.csv");
    //print_basis_curve_to_csv(spline_1d, 4, 1.0, "basis_span4.csv");
    //print_basis_curve_to_csv(spline_1d, 10, 7.0, "basis_span10.csv");
    //print_basis_curve_to_csv(spline_1d, 11, 8.0, "basis_span11.csv");

// ===== 2D Natural Spline Test =====
    std::cout << "\nStarting 2D spline test...\n";

    // Grid size
    int Nx = 11;
    int Ny = 7;
    std::vector<int> shape2d = { Nx, Ny };

    // Create 2D data: f(x,y) = sin(2*pi*x/Nx) * cos(2*pi*y/Ny)
    std::vector<REAL> data2d(Nx * Ny);
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            data2d[iy * Nx + ix] = std::sin(2.0f * M_PI * ix / REAL(Nx - 1)) * std::cos(2.0f * M_PI * iy / REAL(Ny - 1));
        }
    }

    // Build 2D spline
    Natural_BSpline_ND spline_2d(shape2d, data2d.data());

    // Open CSVs
    std::ofstream ground_truth_out("cpp_2d_ground_truth.csv");
    std::ofstream spline2d_out("cpp_2d_spline_output.csv");
    ground_truth_out << "x,y,f_true\n";
    spline2d_out << "x,y,f_interp\n";

    // Fine sampling resolution
    const REAL dx = 0.1f;
    const REAL dy = 0.1f;

    // Sample over [0, Nx-1] x [0, Ny-1]
    for (REAL y = 0.0f; y <= REAL(Ny - 1); y += dy) {
        for (REAL x = 0.0f; x <= REAL(Nx - 1); x += dx) {
            std::vector<REAL> point = { x, y };

            // True function value
            REAL f_true = std::sin(2.0f * M_PI * x / REAL(Nx - 1)) * std::cos(2.0f * M_PI * y / REAL(Ny - 1));
            ground_truth_out << x << "," << y << "," << f_true << "\n";

            // Spline interpolation
            REAL f_interp = spline_2d.evaluate(point);
            spline2d_out << x << "," << y << "," << f_interp << "\n";
        }
    }

    ground_truth_out.close();
    spline2d_out.close();
    std::cout << "2D ground truth written to cpp_2d_ground_truth.csv\n";
    std::cout << "2D spline output written to cpp_2d_spline_output.csv\n";


    // Set up small test
    Nx = 3;
    Ny = 3;
    shape2d = { Nx, Ny };
    std::vector<REAL> new_data2d(Nx * Ny);

    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            new_data2d[iy * Nx + ix] = REAL(ix + 10 * iy);
        }
    }

    Natural_BSpline_ND new_spline_2d(shape2d, new_data2d.data());

    // Evaluate at some points
    std::vector<REAL> x_query = { 0.0, 1.0, 2.0 };
    std::vector<REAL> y_query = { 0.0, 1.0, 2.0 };

    for (REAL y : y_query) {
        for (REAL x : x_query) {
            std::vector<REAL> point = { x, y };
            REAL val = new_spline_2d.evaluate(point);
            std::cout << "spline(" << x << ", " << y << ") = " << val << "\n";
        }
    }



    return 0;
}