#ifndef NATURAL_BSPLINE_1D_H_INCLUDED
#define NATURAL_BSPLINE_1D_H_INCLUDED


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


class SPLINE_API Natural_BSpline_1D {
public:
    static constexpr int SPLINE_DEGREE = 3;

    // Constructors
    Natural_BSpline_1D();
    Natural_BSpline_1D(int num_data_points);
    Natural_BSpline_1D(int num_data_points, const REAL* values);

    // Destructor
    ~Natural_BSpline_1D();

    // Coefficient calculation
    void calculate_coefficients(const REAL* values);

    // Evaluation
    REAL evaluate(REAL x) const;
    REAL evaluate_derivative(REAL x) const;

    void calculate_values(const REAL* x_values, int num_points, REAL* output_values) const;
    void calculate_derivatives(const REAL* x_values, int num_points, REAL* output_derivatives) const;

    // Runtime toggle for fast basis evaluation
    bool get_fast_evaluation();
    void set_fast_evaluation(bool flag);

    // Metadata
    int num_knots() const { return num_knots_; }
    REAL knot(int i) const { return knots_[i]; }
    int num_control_points() const { return num_control_points_; }
    int num_data_points() const { return num_data_points_; }
    bool coefficients_ready() const { return coefficients_ready_; }
    const REAL* coefficients() const { return coefficients_; }
    const Eigen::MatrixXd& get_collocation_matrix() const { return A_; }
    const Eigen::VectorXd& get_rhs_vector() const { return b_; }

    // Debug
    void debug_compare_basis(REAL x);

    void evaluate_basis(const REAL x, const int span, const int derivative_order, REAL* result) const;

private:
    int num_data_points_;
    int num_control_points_;
    int num_knots_;

    bool knots_ready_;
    bool coefficients_ready_;

    REAL* knots_;
    REAL* coefficients_;

    bool use_fast_evaluation_;

    // Stored as double precision for compatibility with Eigen::MatrixXd
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;

    void init_knot_vector();
    int find_knot_interval(REAL x) const;

    void evaluate_fast_cubic_basis(const REAL x, const int span, REAL* basis) const;
    void evaluate_fast_cubic_basis_derivative(const REAL x, const int span, REAL* dbasis) const;
};



#endif