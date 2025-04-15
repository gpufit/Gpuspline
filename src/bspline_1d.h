#ifndef BSPLINE_1D_H_INCLUDED
#define BSPLINE_1D_H_INCLUDED

class BSpline1D {
public:
    static const int SPLINE_DEGREE;
    int num_control_points_;
    REAL* coefficients_;

    // Constructors
    BSpline1D(int num_control_points);
    BSpline1D(int num_control_points, const REAL* values);

    // Destructor
    ~BSpline1D();

    // Evaluation
    REAL evaluate(REAL x) const;
    REAL evaluate_derivative(REAL x) const;
    void calculate_values(const REAL* x_values, int num_points, REAL* output_values) const;
    void calculate_derivatives(const REAL* x_values, int num_points, REAL* output_derivatives) const;


    // Coefficient handling
    void calculate_coefficients(const REAL* values);
    void init_zero();

private:
    std::vector<REAL> knots_;

    // Internal methods
    void init_not_a_knot();
    int find_knot_interval(REAL x) const;
    void evaluate_basis(REAL x, int interval, REAL* basis) const;
    void evaluate_basis_derivative(REAL x, int interval, REAL* dbasis) const;
};

#endif