#include <iostream>
#include "spline_classes.h"

int main() {

    //const int N = 10;
    //REAL data[N] = {0.0, 0.1, 0.4, 0.9, 1.6, 2.5, 3.6, 4.9, 6.4, 8.1};

    // Define N and test data
    const int N = 10;
    float data[N];
    for (int i = 0; i < N; ++i) {
        data[i] = std::sin(2.0f * M_PI * i / N);
    }

    // Construct spline
    BSpline1D spline(N, data);

    // Evaluate spline and its derivative at regular intervals
    std::cout << "x\tf(x)\tf'(x)\n";
    for (REAL x = 0.0; x <= 9.0; x += 0.5) {
        REAL fx = spline.evaluate(x);
        REAL dfx = spline.evaluate_derivative(x);
        std::cout << x << "\t" << fx << "\t" << dfx << "\n";
    }


    std::cout << "C++ spline coefficients:\n";
    for (int i = 0; i < spline.num_control_points(); ++i) {
        std::cout << i << ": " << spline.coefficient_at(i) << "\n";
    }


    return 0;
}