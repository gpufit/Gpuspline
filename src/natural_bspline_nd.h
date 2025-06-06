#ifndef NATURAL_BSPLINE_ND_H_INCLUDED
#define NATURAL_BSPLINE_ND_H_INCLUDED


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


class SPLINE_API Natural_BSpline_ND {
public:
    static constexpr int MAX_NDIMS = 8;
    static constexpr int SPLINE_DEGREE = 3;

    // Constructors
    Natural_BSpline_ND(const std::vector<int>& grid_shape, const REAL* values);

    Natural_BSpline_ND(const std::vector<int>& data_dims,
        const std::vector<REAL>& coefficients,
        const std::vector<std::array<int, MAX_NDIMS>>& multi_indices);

    // Destructor
    ~Natural_BSpline_ND();

    // Runtime toggle for fast basis evaluation
    bool get_fast_evaluation();
    void set_fast_evaluation(bool flag);

    // Evaluate the spline at a given N-D point (x)
    REAL evaluate(const std::vector<REAL>& x) const;

    // Evaluate partial derivative with respect to axis `dim`
    REAL evaluate_derivative(const std::vector<REAL>& x, int dim) const;

    void evaluate_batch(const std::vector<std::vector<REAL>>& points, std::vector<REAL>& output) const;

    void evaluate_batch_derivatives(const std::vector<std::vector<REAL>>& points, int dimension, std::vector<REAL>& output) const;

    // Accessors
    int dimension() const { return num_dims_; }
    const std::vector<int>& shape() const { return data_dims_; }
    const std::vector<REAL>& coefficients() const { return coefficients_; }

private:
    int num_dims_;
    int total_data_points_;
    int total_control_points_;

    bool use_fast_evaluation_;

    std::vector<int> data_dims_;
    //std::vector<int> data_strides_;
    std::vector<int> control_point_dims_;
    //std::vector<int> control_point_strides_;

    std::vector<std::array<int, MAX_NDIMS>> multi_indices_; // precomputed multi-index combinations

    std::vector<REAL> coefficients_;               // Flattened tensor of spline coeffs
    std::vector<std::vector<REAL>> knot_vectors_;  // one knot vector per dimension

    // Internal utilities
    void init_knot_vectors();
    void calculate_coefficients(const REAL* values);
    int find_knot_interval(REAL x, int dimension) const;

    // Basis function evaluation
    static void evaluate_basis_fitpack(REAL x, const std::vector<REAL>& knots, int span, int derivative_order, REAL* basis);

    void evaluate_fast_cubic_basis(const REAL x, const int span, const int axis, REAL* basis) const;
    void evaluate_fast_cubic_basis_derivative(const REAL x, const int span, const int axis, REAL* dbasis) const;

    // Flatten multi-index (i0, i1, ..., id-1) to linear index
    void build_multi_indices();
    //int flatten_multi_index(int i, int slice_idx) const;
    //int flatten_multi_index_control(int i, int slice_idx) const;
};

#endif