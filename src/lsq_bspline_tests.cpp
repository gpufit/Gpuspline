#include <iostream>
#include <fstream>
#include "spline_classes.h"

int main() {

    //const int N = 10;
    //REAL data[N] = {0.0, 0.1, 0.4, 0.9, 1.6, 2.5, 3.6, 4.9, 6.4, 8.1};

    // Define N and test data
    const int N = 10;
    REAL data[N];
    for (int i = 0; i < N; ++i) {
        data[i] = std::sin(2.0f * M_PI * i / REAL(N));
    }

    // Construct least-squares spline
    LSQ_BSpline_1D spline(N, data);

    //std::cout << "x\tf(x)\tf'(x)\n";
    //for (REAL x = 0.0; x <= 9.0; x += 0.5f) {
    //    REAL fx = spline.evaluate(x);
    //    REAL dfx = spline.evaluate_derivative(x);
    //    std::cout << x << "\t" << fx << "\t" << dfx << "\n";
    //}

    //std::cout << "\nC++ spline coefficients:\n";
    //for (int i = 0; i < spline.num_control_points(); ++i) {
    //    std::cout << i << ": " << spline.coefficient_at(i) << "\n";
    //}

    // Spline output at fine resolution
    std::ofstream spline_out("cpp_spline_output.csv");
    spline_out << "x,f_cpp,df_cpp\n";
    REAL max_x = REAL(N-1);
    for (REAL x = 0.0; x <= max_x; x += 0.1f) {
        REAL fx = spline.evaluate(x);
        REAL dfx = spline.evaluate_derivative(x);
        spline_out << x << "," << fx << "," << dfx << "\n";
    }
    spline_out.close();

    // Ground truth at original N points
    std::ofstream truth_out("cpp_ground_truth.csv");
    truth_out << "x,f_true\n";
    for (int i = 0; i < N; ++i) {
        REAL x = static_cast<REAL>(i);
        REAL f_true = data[i];
        truth_out << x << "," << f_true << "\n";
    }
    truth_out.close();

    return 0;
}