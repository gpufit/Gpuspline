#include "spline_classes.h"
#include "equation_system.h"

EquationSystem::EquationSystem(std::size_t const N) :
    N_(N),
    matrix_(N * N),
    vector_(N),
    solution_(N)
{}

EquationSystem::~EquationSystem() {}

void EquationSystem::prepare_matrix()
{
    set_matrix();
    decompose_matrix();
}

void EquationSystem::decompose_matrix()
{
    for (std::size_t i = 0; i < N_; i++)
    {
        for (std::size_t j = i + 1; j < N_; j++)
        {
            matrix_[j * N_ + i] /= matrix_[i* N_ + i];

            for (std::size_t k = i + 1; k < N_; k++)
            {
                matrix_[j * N_ + k] -= matrix_[j * N_ + i] * matrix_[i* N_ + k];
            }
        }
    }
}

void EquationSystem::solve(REAL * solution)
{
    for (std::size_t i = 0; i < N_; i++)
    {
        solution[i] = vector_[i];

        for (std::size_t k = 0; k < i; k++)
        {
            solution[i] -= matrix_[i* N_ + k] * solution[k];
        }
    }

    for (std::size_t i = N_ - 1; i >= 0 && i < N_; i--)
    {
        for (std::size_t k = i + 1; k < N_; k++)
        {
            solution[i] -= matrix_[i* N_ + k] * solution[k];
        }

        solution[i] = solution[i] / matrix_[i* N_ + i];
    }
}

void EquationSystem::solve()
{
    for (std::size_t i = 0; i < N_; i++)
    {
        solution_[i] = vector_[i];

        for (std::size_t k = 0; k < i; k++)
        {
            solution_[i] -= matrix_[i* N_ + k] * solution_[k];
        }
    }

    for (std::size_t i = N_ - 1; i >= 0 && i < N_; i--)
    {
        for (std::size_t k = i + 1; k < N_; k++)
        {
            solution_[i] -= matrix_[i* N_ + k] * solution_[k];
        }

        solution_[i] = solution_[i] / matrix_[i* N_ + i];
    }
}

void EquationSystem::resize(std::size_t N)
{
    if (N != N_)
    {
        N_ = N;

        matrix_.resize(N_ * N_);
        vector_.resize(N_);
        solution_.resize(N_);

        std::fill(matrix_.begin(), matrix_.end(), 0.f);

        prepare_matrix();
    }
}

EquationSystem_1D::EquationSystem_1D(std::size_t const N) : EquationSystem(N)
{
    if (N > 0)
        prepare_matrix();
}

EquationSystem_1D::EquationSystem_1D() : EquationSystem_1D(0)
{}

// see Wolfram page, creates a matrix (just containing constant numbers)
void EquationSystem_1D::set_matrix()
{
    matrix_[0 * N_ + 0] = 2;
    matrix_[(N_ - 1) * N_ + (N_ - 1)] = 2;
    matrix_[0 * N_ + 1] = 1;
    matrix_[(N_ - 1) * N_ + (N_ - 2)] = 1;

    for (std::size_t i = 1; i < N_ - 1; i++)
    {
        matrix_[i * N_ + (i - 1)] = 1;
        matrix_[i * N_ + i] = 4;
        matrix_[i * N_ + (i + 1)] = 1;
    }
}

// see Wolfram page, creates vector b from data
void EquationSystem_1D::set_vector(REAL const * data)
{
    vector_[0] = 3 * (data[1] - data[0]);
    vector_[N_ - 1] = 3 * (data[N_ - 1] - data[N_ - 2]);

    for (std::size_t i = 1; i < N_ - 1; i++)
    {
        vector_[i] = 3 * (data[i + 1] - data[i - 1]);
    }
}

EquationSystem_2D::EquationSystem_2D() : EquationSystem(Spline2D::n_coefficients_per_point)
{
    prepare_matrix();
}

// compute matrix holding x^ay^b for the 16x16 points in the interval
void EquationSystem_2D::set_matrix()
{
    for (std::size_t i = 0; i < 4; i++)
    {
        REAL dx = static_cast<REAL>(i) / 3;
        for (std::size_t j = 0; j < 4; j++)
        {
            REAL dy = static_cast<REAL>(j) / 3;
            for (std::size_t k = 0; k < 4; k++)
            {
                for (std::size_t l = 0; l < 4; l++)
                {
                    // k and l fastest changing, i and j slowest
                    matrix_[(i * 4 + j) * 16 + (k * 4 + l)]
                        = static_cast<REAL>(std::pow(dx, k) * std::pow(dy, l));
                }
            }
        }
    }
}

// vector holding the 16 data values in the interval
void EquationSystem_2D::set_vector(
    std::vector<Spline1D> & splines,
    std::size_t const i,
    std::size_t const j)
{
    // i = spline_interval along x-axis
    // j = spline_interval along y-xis
    for (std::size_t k = 0; k < 4; k++)
    {
        Spline1D & spline = splines[3 * i + k];
        for (std::size_t l = 0; l < 4; l++)
        {
            REAL y = static_cast<REAL>(j) + static_cast<REAL>(l) / 3;
            // order is y-axis first, then x-axis
            vector_[k * 4 + l] = spline.calculate_value(y);
        }
    }
}

EquationSystem_3D::EquationSystem_3D() : EquationSystem(Spline3D::n_coefficients_per_point)
{
    prepare_matrix();
}

// compute matrix holding x^ay^bz^c for the 64x64 data values in the interval
void EquationSystem_3D::set_matrix()
{
    for (std::size_t i = 0; i < 4; i++)
    {
        REAL const dx = static_cast<REAL>(i) / 3;
        for (std::size_t j = 0; j < 4; j++)
        {
            REAL const dy = static_cast<REAL>(j) / 3;
            for (std::size_t k = 0; k < 4; k++)
            {
                REAL const dz = static_cast<REAL>(k) / 3;
                for (std::size_t l = 0; l < 4; l++)
                {
                    for (std::size_t m = 0; m < 4; m++)
                    {
                        for (std::size_t n = 0; n < 4; n++)
                        {
                            matrix_[(i * 16 + j * 4 + k) * 64 + (l * 16 + m * 4 + n)]
                                = static_cast<REAL>
                                (std::pow(dx, l) * std::pow(dy, m) * std::pow(dz, n));
                        }
                    }
                }
            }
        }
    }
}

// vector holding the 64 data values in the interval
void EquationSystem_3D::set_vector(
    std::vector<Spline1D> & splines,
    std::size_t const interpolated_size_x,
    std::size_t const interpolated_size_y,
    std::size_t const i,
    std::size_t const j,
    std::size_t const k)
{
    // (i,j,k) spline interval index in 3D
    for (std::size_t m = 0; m < 4; m++)
    {
        for (std::size_t n = 0; n < 4; n++)
        {
            Spline1D & spline
                = splines[i * 3 * interpolated_size_y + j * 3 + m * interpolated_size_y + n];

            for (std::size_t o = 0; o < 4; o++)
            {
                REAL z = static_cast<REAL>(k) + static_cast<REAL>(o) / 3;
                vector_[m * 16 + n * 4 + o] = spline.calculate_value(z);
            }
        }
    }
}