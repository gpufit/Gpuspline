#include "spline_classes.h"
#include "spline.h"
#include <iostream>
#include <vector>

// internal test/example

bool equal(REAL a, REAL b)
{
    return std::abs(a - b) < std::numeric_limits<REAL>::epsilon() * 10;
}

int test_spline_2d()
{
    // data size
    int const size_x = 10;
    int const size_y = 12;

    // interpolated values size
    REAL const delta = 1.f;
    int const spline_size_x = static_cast<std::size_t>(size_x / delta);
    int const spline_size_y = static_cast<std::size_t>(size_y / delta);

    std::vector<REAL> psf(size_x * size_y);
    std::vector<REAL> spline(spline_size_x * spline_size_y);
    std::vector<REAL> x_values(spline_size_x);
    std::vector<REAL> y_values(spline_size_y);
    std::vector<REAL> coefficients((size_x - 1) * (size_y - 1) * 16);

    // gaussian peak parameters
    REAL const a = 1;
    REAL const x0 = (REAL(size_x) - 1) / 2;
    REAL const y0 = (REAL(size_y) - 3) / 2;
    REAL const s = 1.5f;
    REAL const b = 0;

    for (int point_index_y = 0; point_index_y < size_y; point_index_y++)
    {
        for (int point_index_x = 0; point_index_x < size_x; point_index_x++)
        {
            int const point_index = point_index_y * size_x + point_index_x;
            REAL const argx = ((point_index_x - x0)*(point_index_x - x0)) / (2 * s * s);
            REAL const argy = ((point_index_y - y0)*(point_index_y - y0)) / (2 * s * s);
            REAL const ex = exp(-argx) * exp(-argy);
            psf[point_index] = a * ex + b;
        }
    }

    // interpolation x values
    for (int point_index_x = 0; point_index_x < spline_size_x; point_index_x++)
        x_values[point_index_x] = point_index_x * delta;

    // interpolation y values
    for (int point_index_y = 0; point_index_y < spline_size_y; point_index_y++)
        y_values[point_index_y] = point_index_y * delta;

    // calculate spline coefficients
    calculate_coefficients_2d(psf.data(), size_x, size_y, coefficients.data());
    
    // calculate spline values
    calculate_values_2d(coefficients.data(), size_x-1, size_y-1, size_x, size_y, x_values.data(), y_values.data(), spline.data());

    // test if data is the same for delta == 1
    bool test_passed = true;
    if (equal(delta, 1))
    {
        for (int point_index_y = 0; point_index_y < size_y; point_index_y++)
        {
            for (int point_index_x = 0; point_index_x < size_x; point_index_x++)
            {
                int const point_index = point_index_y * size_x + point_index_x;
                if (!equal(psf[point_index], spline[point_index]))
                {
                    std::cout << "index y=" << point_index_y << " x=" << point_index_x << " psf=" << psf[point_index] << " not equal to spline=" << spline[point_index] << std::endl;
                    test_passed = false;
                }
            }
        }
    }
    if (test_passed)
    {
        std::cout << std::endl << "Test passed!" << std::endl << std::endl;
        return 1;
    }
    else
    {
        std::cout << std::endl << "Test failed!" << std::endl << std::endl;
        return -1;
    }

    return 0;
}

int test_spline_3d()
{
    int const size_x = 10;
    int const size_y = 12;
    int const size_z = 17;
    int const size_xy = size_x * size_y;

    REAL const delta = .5f;
    int const spline_size_x = static_cast<std::size_t>(size_x / delta);
    int const spline_size_y = static_cast<std::size_t>(size_y / delta);
    int const spline_size_z = static_cast<std::size_t>(size_z / delta);

    std::vector<REAL> psf(size_x * size_y * size_z);
    std::vector<REAL> spline(spline_size_x * spline_size_y * spline_size_z);
    std::vector<REAL> x_values(spline_size_x);
    std::vector<REAL> y_values(spline_size_y);
    std::vector<REAL> z_values(spline_size_z);
    std::vector<REAL> coefficients((size_x - 1) * (size_y - 1) * (size_z - 1) * 64);

    // gaussian peak parameters
    REAL const a = 1;
    REAL const x0 = (REAL(size_x) - 1) / 2;
    REAL const y0 = (REAL(size_y) - 1) / 2;
    REAL const z0 = (REAL(size_z) - 1) / 2;
    REAL const s = 1.5f;
    REAL const b = 0;

    for (int point_index_z = 0; point_index_z < size_z; point_index_z++)
    {
        for (int point_index_y = 0; point_index_y < size_y; point_index_y++)
        {
            for (int point_index_x = 0; point_index_x < size_x; point_index_x++)
            {
                int const point_index
                    = point_index_z * size_xy + point_index_y * size_x + point_index_x;

                REAL const argx = ((point_index_x - x0)*(point_index_x - x0)) / (2 * s * s);
                REAL const argy = ((point_index_y - y0)*(point_index_y - y0)) / (2 * s * s);
                REAL const argz = ((point_index_z - z0)*(point_index_z - z0)) / (2 * s * s);
                REAL const ex = exp(-argx) * exp(-argy) * exp(-argz);
                psf[point_index] = a * ex + b;
            }
        }
    }

    // interpolation x values
    for (int point_index_x = 0; point_index_x < spline_size_x; point_index_x++)
        x_values[point_index_x] = point_index_x * delta;

    // interpolation y values
    for (int point_index_y = 0; point_index_y < spline_size_y; point_index_y++)
        y_values[point_index_y] = point_index_y * delta;

    // interpolation z values
    for (int point_index_z = 0; point_index_z < spline_size_z; point_index_z++)
        z_values[point_index_z] = point_index_z * delta;

    // calculate spline coefficients
    calculate_coefficients_3d(
        psf.data(),
        size_x, size_y, size_z,
        coefficients.data());

    // calculate spline values
    calculate_values_3d(
        coefficients.data(),
        size_x-1, size_y-1, size_z-1,
        size_x, size_y, size_z,
        x_values.data(), y_values.data(), z_values.data(),
        spline.data());

    return 1;
}

int main()
{
    test_spline_2d();
    test_spline_3d();
}