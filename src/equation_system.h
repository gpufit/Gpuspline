#ifndef EQUATION_SYSTEM_INCLUDED
#define EQUATION_SYSTEM_INCLUDED

#include <vector>

class Spline1D;

class EquationSystem
{
public:
    EquationSystem(std::size_t const N);

    ~EquationSystem();

    void prepare_matrix();

    void decompose_matrix();

    void solve(REAL * solution);
    void solve();

    void resize(std::size_t const N);

    std::size_t N_;
    std::vector<REAL> matrix_;
    std::vector<REAL> vector_;
    std::vector<REAL> solution_;

    virtual void set_matrix() = 0;
};

class EquationSystem_1D : public EquationSystem
{
public:
    EquationSystem_1D(std::size_t const N);
    EquationSystem_1D();

    void set_vector(REAL const * data);

    void set_matrix();
};

class EquationSystem_2D : public EquationSystem
{
public:
    EquationSystem_2D();

    void set_vector(
        std::vector<Spline1D> & splines,
        std::size_t const i,
        std::size_t const j);

    void set_matrix();
};

class EquationSystem_3D : public EquationSystem
{
public:
    EquationSystem_3D();

    void set_vector(
        std::vector<Spline1D> & splines,
        std::size_t const interpolated_size_x,
        std::size_t const interpolated_size_y,
        std::size_t const i,
        std::size_t const j,
        std::size_t const k);

    void set_matrix();
};

class EquationSystem_4D : public EquationSystem
{
public:
    EquationSystem_4D();

    void set_vector(
        std::vector<Spline1D> & splines,
        std::size_t const interpolated_size_x,
        std::size_t const interpolated_size_y,
        std::size_t const interpolated_size_z,
        std::size_t const i,
        std::size_t const j,
        std::size_t const k, 
        std::size_t const l);

    void set_matrix();
};

extern EquationSystem_1D equation_system_1d;
extern EquationSystem_2D equation_system_2d;
extern EquationSystem_3D equation_system_3d;
extern EquationSystem_4D equation_system_4d;

#endif // !EQUATION_SYSTEM_INCLUDED