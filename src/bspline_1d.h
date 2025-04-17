#ifndef BSPLINE_1D_H_INCLUDED
#define BSPLINE_1D_H_INCLUDED


#ifdef _WIN32
#  ifdef SPLINES_BUILD_DLL
#    define SPLINE_API __declspec(dllexport)
#  else
#    define SPLINE_API __declspec(dllimport)
#  endif
#else
#  define SPLINE_API
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class SPLINE_API BSpline1D {
public:
    static constexpr int SPLINE_DEGREE = 3;
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

    // Get the number of control points
    int num_control_points() const { return num_control_points_; }

    // Get the i-th coefficient
    REAL coefficient_at(int i) const {
        if (i >= 0 && i < num_control_points_) {
            return coefficients_[i];
        }
        else {
            return 0.0;  // Or throw if you prefer strict bounds
        }
    }

private:
    std::vector<REAL> knots_;

    // Internal methods
    void init_not_a_knot();
    void init_uniform_clamped();
    int find_knot_interval(REAL x) const;
    void evaluate_basis(REAL x, int interval, REAL* basis) const;
    void evaluate_basis_derivative(REAL x, int interval, REAL* dbasis) const;
};

#endif