#include "spline_classes.h"
#include "equation_system.h"

// TODO better documentation of what each function does and requires and in which order these functions should be called

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

// compute matrix holding x^a y^b for the 16x16 points in the interval
void EquationSystem_2D::set_matrix()
{
    for (std::size_t m1 = 0; m1 < 4; m1++)
    {
        REAL dx = static_cast<REAL>(m1) / 3.0;

        for (std::size_t n1 = 0; n1 < 4; n1++)
        {
            REAL dy = static_cast<REAL>(n1) / 3.0;

            for (std::size_t m2 = 0; m2 < 4; m2++)
            {
                for (std::size_t n2 = 0; n2 < 4; n2++)
                {
                    // m2 and n2 fastest changing, m1 and n1 slowest
                    matrix_[(m1 * 4 + n1) * 16 + (m2 * 4 + n2)]
                        = static_cast<REAL>(std::pow(dx, m2) * std::pow(dy, n2));
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
    // j = spline_interval along y-axis
    for (std::size_t k = 0; k < 4; k++)
    {
        Spline1D & spline = splines[3 * i + k];
        for (std::size_t m = 0; m < 4; m++)
        {
            REAL y = static_cast<REAL>(j) + static_cast<REAL>(m) / 3.0;
            // order is y-axis first, then x-axis
            vector_[k * 4 + m] = spline.calculate_value(y);
        }
    }
}

EquationSystem_3D::EquationSystem_3D() : EquationSystem(Spline3D::n_coefficients_per_point)
{
    prepare_matrix();
}

// compute matrix holding x^a y^b z^c for the 64x64 data values in the interval
void EquationSystem_3D::set_matrix()
{
    for (std::size_t m1 = 0; m1 < 4; m1++)
    {
        REAL const dx = static_cast<REAL>(m1) / 3.0;

        for (std::size_t n1 = 0; n1 < 4; n1++)
        {
            REAL const dy = static_cast<REAL>(n1) / 3.0;

            for (std::size_t o1 = 0; o1 < 4; o1++)
            {
                REAL const dz = static_cast<REAL>(o1) / 3.0;

                for (std::size_t m2 = 0; m2 < 4; m2++)
                {
                    for (std::size_t n2 = 0; n2 < 4; n2++)
                    {
                        for (std::size_t o2 = 0; o2 < 4; o2++)
                        {

                            // m2, n2, o2 fastest changing, m1, n1, o1 slowest
                            matrix_[(m1 * 16 + n1 * 4 + o1) * 64 + (m2 * 16 + n2 * 4 + o2)]
                                = static_cast<REAL>(std::pow(dx, m2) * std::pow(dy, n2) * std::pow(dz, o2));

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
                REAL z = static_cast<REAL>(k) + static_cast<REAL>(o) / 3.0;
                vector_[m * 16 + n * 4 + o] = spline.calculate_value(z);
            }
        }
    }
}


EquationSystem_4D::EquationSystem_4D() : EquationSystem(Spline4D::n_coefficients_per_point)
{
    prepare_matrix();
}

// compute matrix holding x^a y^b z^c u^d for the 256x256 data values in the interval
void EquationSystem_4D::set_matrix()
{

    for (std::size_t m1 = 0; m1 < 4; m1++)
    {
        REAL const dx = static_cast<REAL>(m1) / 3.0;

        for (std::size_t n1 = 0; n1 < 4; n1++)
        {
            REAL const dy = static_cast<REAL>(n1) / 3.0;

            for (std::size_t o1 = 0; o1 < 4; o1++)
            {
                REAL const dz = static_cast<REAL>(o1) / 3.0;

                for (std::size_t p1 = 0; p1 < 4; p1++)
                {
                    REAL const du = static_cast<REAL>(p1) / 3.0;


                    for (std::size_t m2 = 0; m2 < 4; m2++)
                    {
                        for (std::size_t n2 = 0; n2 < 4; n2++)
                        {
                            for (std::size_t o2 = 0; o2 < 4; o2++)
                            {
                                for (std::size_t p2 = 0; p2 < 4; p2++)
                                {

                                    // m2, n2, o2, p2 fastest changing, m1, n1, o1, p1 slowest
                                    matrix_[(m1 * 64 + n1 * 16 + o1 * 4 + p1) * 256 + (m2 * 64 + n2 * 16 + o2 * 4 + p2)]
                                        = static_cast<REAL>(std::pow(dx, m2) * std::pow(dy, n2) * std::pow(dz, o2) * std::pow(du, p2));

                                }
                            }
                        }
                    }


                }
            }
        }
    }

}

// vector holding the 256 data values in the interval
void EquationSystem_4D::set_vector(
    std::vector<Spline1D> & splines,
    std::size_t const interpolated_size_x,
    std::size_t const interpolated_size_y,
    std::size_t const interpolated_size_z,
    std::size_t const i,
    std::size_t const j,
    std::size_t const k, 
    std::size_t const m)
{


    // (i,j,k,m) spline interval index in 4D
    for (std::size_t m = 0; m < 4; m++)
    {
        for (std::size_t n = 0; n < 4; n++)
        {
            for (std::size_t o = 0; o < 4; o++)
            {

                Spline1D & spline
                    = splines[i * 3 * interpolated_size_y * interpolated_size_z + j * 3 * interpolated_size_z + k * 3 + m * interpolated_size_y * interpolated_size_z + n * interpolated_size_z + o];

                for (std::size_t p = 0; p < 4; p++)
                {
                    REAL w = static_cast<REAL>(m)+static_cast<REAL>(p) / 3.0;
                    vector_[m * 64 + n * 16 + o * 4 + p] = spline.calculate_value(w);
                }
            }
        }
    }


}
