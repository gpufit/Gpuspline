#ifndef LSQ_BSPLINE_1D_H_INCLUDED
#define LSQ_BSPLINE_1D_H_INCLUDED


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

class SPLINE_API LSQ_BSpline_1D {
public:

    static constexpr int SPLINE_DEGREE = 3;

    // Constructors
    LSQ_BSpline_1D();  // Default constructor (empty)
    LSQ_BSpline_1D(int num_data_points);  // Zero-initialised coefficients
    LSQ_BSpline_1D(int num_data_points, const REAL* values);  // Fit to data via least-squares

    // Evaluation
    REAL evaluate(REAL x) const;
    REAL evaluate_derivative(REAL x) const;

    void calculate_values(const REAL* x_values, int num_points, REAL* output_values) const;
    void calculate_derivatives(const REAL* x_values, int num_points, REAL* output_derivatives) const;

    // Coefficient calculation (least-squares fit)
    void calculate_coefficients(const REAL* values);  // Fit from num_data_points
    void init_zero();  // Zero-initialise coefficients

    int num_control_points() const { return num_control_points_; }
    REAL coefficient_at(int i) const { return coefficients_[i]; }

    // Destructor
    ~LSQ_BSpline_1D();

private:
    void init_uniform_clamped_knots();  // Clamped knot vector allocation
    void evaluate_fast_basis(REAL t, REAL* basis) const;
    void evaluate_fast_basis_derivative(REAL t, REAL* dbasis) const;
    void evaluate_basis_for_fitting(REAL x, int span, REAL* basis, const REAL* knots);

    // Members
    int num_data_points_;      // Number of input samples
    int num_control_points_;   // = num_data_points + degree + 3
    int num_knots_;            // = num_control_points + degree + 1
    REAL* coefficients_;       // Control points
    REAL* knots_;              // Clamped uniform knot vector
};

#endif