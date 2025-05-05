
#pragma once

#include "definitions.h"  // for REAL

inline void evaluate_fast_cubic_basis_first_span(const REAL t, REAL* basis) {
    REAL t2 = t * t;
    REAL t3 = t2 * t;

    basis[0] = 1.0 - 3.0 * t + 3.0 * t2 - t3;         // B0(x)
    basis[1] = 3.0 * t - 4.5 * t2 + 1.75 * t3;        // B1(x)
    basis[2] = 1.5 * t2 - (11.0 / 12.0) * t3;         // B2(x)
    basis[3] = t3 / 6.0;                              // B3(x)
}

inline void evaluate_fast_cubic_basis_second_span(const REAL t, REAL* basis) {
    REAL t2 = t * t;
    REAL t3 = t2 * t;

    basis[0] = 0.25 - 0.75 * t + 0.75 * t2 - 0.25 * t3;                    // B0(x)
    basis[1] = (7.0 / 12.0) + 0.25 * t - 1.25 * t2 + (7.0 / 12.0) * t3;    // B1(x)
    basis[2] = (1.0 / 6.0) + 0.5 * t + 0.5 * t2 - 0.5 * t3;                // B2(x)
    basis[3] = t3 / 6.0;                                                   // B3(x)
}

inline void evaluate_fast_cubic_basis_interior_span(const REAL t, REAL* basis) {
    REAL t2 = t * t;
    REAL t3 = t2 * t;
    REAL omt = 1.0 - t;

    basis[0] = (1.0 / 6.0) * omt * omt * omt;                               // B0(x)
    basis[1] = (1.0 / 6.0) * (3.0 * t3 - 6.0 * t2 + 4.0);                   // B1(x)
    basis[2] = (1.0 / 6.0) * (-3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0);        // B2(x)
    basis[3] = (1.0 / 6.0) * t3;                                            // B3(x)
}

inline void evaluate_fast_cubic_basis_second_last_span(const REAL t, REAL* basis) {
    REAL omt = 1.0 - t;
    REAL omt2 = omt * omt;
    REAL omt3 = omt2 * omt;

    basis[0] = omt3 / 6.0;                                                       // B0(x)
    basis[1] = (1.0 / 6.0) + 0.5 * omt + 0.5 * omt2 - 0.5 * omt3;                // B1(x)
    basis[2] = (7.0 / 12.0) + 0.25 * omt - 1.25 * omt2 + (7.0 / 12.0) * omt3;    // B2(x)
    basis[3] = 0.25 - 0.75 * omt + 0.75 * omt2 - 0.25 * omt3;                    // B3(x)
}

inline void evaluate_fast_cubic_basis_last_span(REAL t, REAL* basis) {
    REAL omt = 1.0 - t;
    REAL omt2 = omt * omt;
    REAL omt3 = omt2 * omt;

    basis[0] = omt3 / 6.0;                                  // B0(x)
    basis[1] = 1.5 * omt2 - (11.0 / 12.0) * omt3;           // B1(x)
    basis[2] = 3.0 * omt - 4.5 * omt2 + 1.75 * omt3;        // B2(x)
    basis[3] = 1.0 - 3.0 * omt + 3.0 * omt2 - omt3;         // B3(x)
}

inline void evaluate_fast_cubic_basis_derivative_first_span(const REAL t, REAL* dbasis) {
    REAL t2 = t * t;

    dbasis[0] = -3.0 + 6.0 * t - 3.0 * t2;
    dbasis[1] = 3.0 - 9.0 * t + 5.25 * t2;
    dbasis[2] = 3.0 * t - 2.75 * t2;
    dbasis[3] = 0.5 * t2;
}

inline void evaluate_fast_cubic_basis_derivative_second_span(const REAL t, REAL* dbasis) {
    REAL t2 = t * t;

    dbasis[0] = -0.75 + 1.5 * t - 0.75 * t2;
    dbasis[1] = 0.25 - 2.5 * t + 1.75 * t2;
    dbasis[2] = 0.5 + t - 1.5 * t2;
    dbasis[3] = 0.5 * t2;
}

inline void evaluate_fast_cubic_basis_derivative_interior_span(const REAL t, REAL* dbasis) {
    REAL t2 = t * t;
    REAL omt = 1.0 - t;

    dbasis[0] = -0.5 * omt * omt;
    dbasis[1] = 1.5 * t2 - 2.0 * t;
    dbasis[2] = -1.5 * t2 + t + 0.5;
    dbasis[3] = 0.5 * t2;
}

inline void evaluate_fast_cubic_basis_derivative_second_last_span(const REAL t, REAL* dbasis) {
    REAL omt = 1.0 - t;
    REAL omt2 = omt * omt;

    dbasis[0] = -0.5 * omt2;
    dbasis[1] = -0.5 - omt + 1.5 * omt2;
    dbasis[2] = -0.25 + 2.5 * omt - 1.75 * omt2;
    dbasis[3] = 0.75 - 1.5 * omt + 0.75 * omt2;
}

inline void evaluate_fast_cubic_basis_derivative_last_span(const REAL t, REAL* dbasis) {
    REAL omt = 1.0 - t;
    REAL omt2 = omt * omt;

    dbasis[0] = -0.5 * omt2;
    dbasis[1] = -3.0 * omt + 2.75 * omt2;
    dbasis[2] = -3.0 + 9.0 * omt - 5.25 * omt2;
    dbasis[3] = 3.0 - 6.0 * omt + 3.0 * omt2;
}