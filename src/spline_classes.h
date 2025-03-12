#ifndef SPLINE_CLASSES_H_INCLUDED
#define SPLINE_CLASSES_H_INCLUDED

#include "definitions.h"
#include "math_utils.h"
#include "equation_system.h"
#include <vector>
#include <cmath>
#include <cstring>  // For memset
#include <cstddef>
#include <limits>
#include <algorithm>


// TODO maybe we also want to calculate derivatives

// TODO SplineXD can all derive from a common spline base class

// implementation of 1-3D cubic spline interpolation as partly described in http://mathworld.wolfram.com/CubicSpline.html

class Spline1D
{ // 1D cubic interpolating spline calculations
public:
    // end point condition is that 1st and 2nd derivative of spline is zero at end points (?)
    
    // data_size must be larger 1
    
    Spline1D(std::size_t const data_size);
    
    Spline1D(
        std::size_t const n_intervals,
        REAL * coefficients);
    
    ~Spline1D();

    // must be called before all below, data must be equally spaced
    void initialize(REAL const* data);

    // must be called before all below, spline interval equals data spacing
    void calculate_coefficients(REAL * coefficients);

    // given a list of x values calls the function below iteratively
    void calculate_values(
        REAL * spline_values,
        REAL const * x_values,
        std::size_t size_x);

    // calculates a spline value ( outside of range (x < 0 | x > number data points - 1) gives 0 + offset)
    REAL calculate_value(REAL const x);

    void interpolate(
        REAL * interpolated_data,
        REAL const * data,
        REAL const * x_values,
        std::size_t const size_x);

    // change the ordering of a coefficients array
    void Spline1D::convert_csaps_coefficients(
        REAL * csaps_coefficients,
        std::size_t const n_spline_intervals,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

public:
    // TODO rename to coefficients_per_interval
    static std::size_t const n_coefficients_per_point = 4;

private:
    // holds the spline coefficients (number of data points - 1) * 4 coefficients
    // coefficients are stored by order first (increasing), then by interval
    REAL * coefficients_;
    REAL const * data_;
    std::size_t const data_size_;
    // number of spline points (intervals) is data_size - 1
    std::size_t const n_intervals_;
    bool coefficients_calculated_;
    bool data_initialized_;
};

// TODO this assumes quadratic areas in 2D, make it so that rectangular is allowed (dimension vector)
class Spline2D
{ // 2D cubic interpolating spline
public:
    // end conditions ? (as in https://github.com/ZhuangLab/storm-analysis/blob/master/storm_analysis/spliner/spline2D.py)

    // TODO instead of giving size_x, size_y give std::vector<std::size_t> dimensions and check for correct length inside

    // data_sizes must be larger 1

    Spline2D(
        std::size_t const data_size_x,
        std::size_t const data_size_y);
    
    Spline2D(
        std::size_t const n_intervals_x,
        std::size_t const n_intervals_y,
        REAL * coefficients);
    
    ~Spline2D();

    // set 2D data (first dimension (x) first ordering)
    void initialize(REAL const* data);

    void calculate_coefficients(REAL * coefficients);

    // x is fastest changing dimension in spline_values
    void calculate_values(
        REAL * spline_values,
        REAL const * x_values,
        REAL const * y_values,
        std::size_t const size_x, // TODO replace size_x, size_y with vector
        std::size_t const size_y);

    REAL calculate_value(
        REAL const x,
        REAL const y);

    void interpolate(
        REAL * interpolated_data,
        REAL const * data,
        REAL const * x_values,
        REAL const * y_values,
        std::size_t const size_x,
        std::size_t const size_y);

    // change the ordering of a coefficients array
    void Spline2D::convert_csaps_coefficients(
        REAL * csaps_coefficients,
        std::size_t const n_spline_intervals_x,
        std::size_t const n_spline_intervals_y,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

public:
    static std::size_t const n_coefficients_per_point = 4 * 4;
    // TODO total number of spline intervals = product of (dimensions - 1)
    std::size_t const n_intervals_;

private:
    // (sqrt(number of data points) - 1) ^2 * 16 coefficients
    // coefficients are sorted first by increasing order (last dimension changes fastest) then by interval (last dimension changes fastest)
    REAL * coefficients_;
    REAL const * data_;
    std::size_t const data_size_x_;
    std::size_t const data_size_y_;
    // TODO intervals (data_sizex - 1)
    std::size_t const n_intervals_x_;
    // TODO intervals (data_sizey - 1)
    std::size_t const n_intervals_y_;
    bool coefficients_calculated_;
    bool data_initialized_;
};

// TODO make this on rectangular grids (dimensions vector)
class Spline3D
{ // 3D cubic interpolating spline
public:

    // end point conditions ? (as in https://github.com/ZhuangLab/storm-analysis/blob/master/storm_analysis/spliner/spline3D.py)

    // TODO use pointer array instead "std::size_t *" - or vector? - or as it is?

    // data_sizes must be larger 1

    Spline3D(
        std::size_t const data_size_x,
        std::size_t const data_size_y,
        std::size_t const data_size_z);

    Spline3D(
        std::size_t const n_intervals_x,
        std::size_t const n_intervals_y,
        std::size_t const n_intervals_z,
        REAL * coefficients);

    ~Spline3D();

    void initialize(REAL const* data);

    void calculate_coefficients(REAL * coefficients);

    void calculate_values(
        REAL * spline_values,
        REAL const * x_values,
        REAL const * y_values,
        REAL const * z_values,
        std::size_t const size_x,
        std::size_t const size_y,
        std::size_t const size_z);

    REAL calculate_value(
        REAL const x,
        REAL const y,
        REAL const z);

    void interpolate(
        REAL * interpolated_data,
        REAL const * data,
        REAL const * x_values,
        REAL const * y_values,
        REAL const * z_values,
        std::size_t const size_x,
        std::size_t const size_y,
        std::size_t const size_z);

    void Spline3D::convert_csaps_coefficients(
        REAL * csaps_coefficients,
        std::size_t const n_spline_intervals_x,
        std::size_t const n_spline_intervals_y,
        std::size_t const n_spline_intervals_z,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

public:
    static std::size_t const n_coefficients_per_point = 4 * 4 * 4;
    std::size_t const n_intervals_;

private:
    // Prod(dimensions - 1) * 64 coefficients, example 
    // order is coefficients in interval first, then
    REAL * coefficients_;
    REAL const * data_;
    std::size_t const data_size_x_;
    std::size_t const data_size_y_;
    std::size_t const data_size_z_;
    std::size_t const n_intervals_x_;
    std::size_t const n_intervals_y_;
    std::size_t const n_intervals_z_;
    bool coefficients_calculated_;
    bool data_initialized_;
};

// TODO make this on rectangular grids (dimensions vector)
class Spline4D
{ // 4D cubic interpolating spline
public:

    // TODO use pointer array instead "std::size_t *" - or vector? - or as it is?

    // data_sizes must be larger 1

    Spline4D(
        std::size_t const data_size_x,
        std::size_t const data_size_y,
        std::size_t const data_size_z,
        std::size_t const data_size_w);

    Spline4D(
        std::size_t const n_intervals_x,
        std::size_t const n_intervals_y,
        std::size_t const n_intervals_z,
        std::size_t const n_intervals_w,
        REAL * coefficients);

    ~Spline4D();

    void initialize(REAL const* data);

    void calculate_coefficients(REAL * coefficients);

    void calculate_values(
        REAL * spline_values,
        REAL const * x_values,
        REAL const * y_values,
        REAL const * z_values,
        REAL const * w_values,
        std::size_t const size_x,
        std::size_t const size_y,
        std::size_t const size_z,
        std::size_t const size_w);

    REAL calculate_value(
        REAL const x,
        REAL const y,
        REAL const z,
        REAL const w);

    void interpolate(
        REAL * interpolated_data,
        REAL const * data,
        REAL const * x_values,
        REAL const * y_values,
        REAL const * z_values,
        REAL const * w_values,
        std::size_t const size_x,
        std::size_t const size_y,
        std::size_t const size_z,
        std::size_t const size_w);

    void Spline4D::convert_csaps_coefficients(
        REAL * csaps_coefficients,
        std::size_t const n_spline_intervals_x,
        std::size_t const n_spline_intervals_y,
        std::size_t const n_spline_intervals_z,
        std::size_t const n_spline_intervals_w,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

public:
    static std::size_t const n_coefficients_per_point = 4 * 4 * 4 * 4;
    std::size_t const n_intervals_;

private:
    // Prod(dimensions - 1) * 64 coefficients, example 
    // order is coefficients in interval first, then
    REAL * coefficients_;
    REAL const * data_;
    std::size_t const data_size_x_;
    std::size_t const data_size_y_;
    std::size_t const data_size_z_;
    std::size_t const data_size_w_;
    std::size_t const n_intervals_x_;
    std::size_t const n_intervals_y_;
    std::size_t const n_intervals_z_;
    std::size_t const n_intervals_w_;
    bool coefficients_calculated_;
    bool data_initialized_;
};


class Spline5D
{ // 5D cubic interpolating spline
public:

    // Constructor with data sizes
    Spline5D(
        std::size_t const data_size_x,
        std::size_t const data_size_y,
        std::size_t const data_size_z,
        std::size_t const data_size_w,
        std::size_t const data_size_v);

    // Constructor with precomputed coefficients
    Spline5D(
        std::size_t const n_intervals_x,
        std::size_t const n_intervals_y,
        std::size_t const n_intervals_z,
        std::size_t const n_intervals_w,
        std::size_t const n_intervals_v,
        REAL * coefficients);

    ~Spline5D();

    // Initialize spline with raw data
    void initialize(REAL const* data);

    // Compute spline coefficients
    void calculate_coefficients(REAL * coefficients);

    // Evaluate the spline at multiple points
    void calculate_values(
        REAL * spline_values,
        REAL const * x_values,
        REAL const * y_values,
        REAL const * z_values,
        REAL const * w_values,
        REAL const * v_values,
        std::size_t const size_x,
        std::size_t const size_y,
        std::size_t const size_z,
        std::size_t const size_w,
        std::size_t const size_v);

    // Evaluate the spline at a single point
    REAL calculate_value(
        REAL const x,
        REAL const y,
        REAL const z,
        REAL const w,
        REAL const v);

    // Perform full interpolation of a 5D dataset
    void interpolate(
        REAL * interpolated_data,
        REAL const * data,
        REAL const * x_values,
        REAL const * y_values,
        REAL const * z_values,
        REAL const * w_values,
        REAL const * v_values,
        std::size_t const size_x,
        std::size_t const size_y,
        std::size_t const size_z,
        std::size_t const size_w,
        std::size_t const size_v);

    // Convert external CSAPS coefficients to internal format
    void convert_csaps_coefficients(
        REAL * csaps_coefficients,
        std::size_t const n_spline_intervals_x,
        std::size_t const n_spline_intervals_y,
        std::size_t const n_spline_intervals_z,
        std::size_t const n_spline_intervals_w,
        std::size_t const n_spline_intervals_v,
        REAL * grid_spacing_array,
        REAL * reordered_coefficients);

public:
    // Number of coefficients per grid point (4^5 for 5D cubic splines)
    static std::size_t const n_coefficients_per_point = 4 * 4 * 4 * 4 * 4;
    std::size_t const n_intervals_;

private:
    // Coefficients storage: (prod(dimensions - 1)) * (4^5) coefficients
    REAL * coefficients_;
    REAL const * data_;
    std::size_t const data_size_x_;
    std::size_t const data_size_y_;
    std::size_t const data_size_z_;
    std::size_t const data_size_w_;
    std::size_t const data_size_v_;
    std::size_t const n_intervals_x_;
    std::size_t const n_intervals_y_;
    std::size_t const n_intervals_z_;
    std::size_t const n_intervals_w_;
    std::size_t const n_intervals_v_;
    bool coefficients_calculated_;
    bool data_initialized_;
};



class BSpline1D {
public:
    // Constructors and Destructor
    BSpline1D(int Nx);
    BSpline1D(int Nx, REAL* coefficients);
    ~BSpline1D();

    // Initialize coefficients from raw data
    void initialize(const REAL* data);

    // Compute spline coefficients
    void calculate_coefficients(REAL* coefficients);

    // given a list of x values calls the function below iteratively
    void calculate_values(
        REAL * bspline_values,
        REAL const * x_values,
        std::size_t size_x);

    // Evaluate the B-spline function
    REAL evaluate(REAL x) const;

    // Evaluate the B-spline function and its first derivative
    void evaluate_with_derivative(REAL x, REAL& function_value, REAL& dx) const;

private:
    // Grid size
    int Nx_;

    // Store total control points
    int N_control_;
    
    // data storage
    REAL const * data_;
    bool data_initialized_;

    // Coefficient storage
    REAL* coefficients_;
    bool coefficients_calculated_;

    // Knot vector storage
    REAL* knots_;

    // Utility functions
    void allocate_knots();
    void compute_basis_functions(REAL* basis_values, REAL x) const;
    void compute_derivative_basis_functions(REAL* basis_derivatives, REAL x) const;
};


class BSpline5D {
public:
    // Constructors and Destructor
    BSpline5D(int Nx, int Ny, int Nz, int Nw, int Nv);
    BSpline5D(int Nx, int Ny, int Nz, int Nw, int Nv, REAL* coefficients);
    ~BSpline5D();

    // Initialize coefficients from raw data
    void initialize(const REAL* data);

    // Compute spline coefficients
    void calculate_coefficients(REAL* coefficients, const REAL* data);

    // Evaluate the B-spline function
    REAL evaluate(REAL x, REAL y, REAL z, REAL w, REAL v) const;

    // Evaluate the B-spline function and its first derivatives
    void evaluate_with_derivatives(REAL x, REAL y, REAL z, REAL w, REAL v,
        REAL& function_value,
        REAL& dx, REAL& dy, REAL& dz, REAL& dw, REAL& dv) const;

private:
    // Grid sizes
    int Nx_, Ny_, Nz_, Nw_, Nv_;

    // Coefficient storage
    REAL* coefficients_;
    bool coefficients_calculated_;
    bool data_initialized_;

    // Utility functions
    void compute_basis_functions(REAL* basis_values, REAL x) const;
    void compute_derivative_basis_functions(REAL* basis_derivatives, REAL x) const;
};


#endif
