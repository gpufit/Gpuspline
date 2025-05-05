
#include <iostream>
#include <cinttypes>
#include <tuple>
#include <vector>
#include <string>
#include <limits>

namespace fitpack {

/*
 * B-spline evaluation routine.
 */
 
void _deBoor_D(const double *t, double x, int k, int ell, int m, double *result);

}
