#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
#include "SpiceUsr.h"
}

using namespace std;

// ============================================================================
// Types and constants
// ============================================================================
constexpr int NX = 6;
constexpr int NSTM = NX * NX;
constexpr int NTOTAL = NX + NSTM;

using Vec3  = array<double, 3>;
using Mat3  = array<array<double, 3>, 3>;
using Mat6  = array<array<double, 6>, 6>;
using State = array<double, NTOTAL>;
using Matrix = vector<vector<double>>;

struct PlanetParams {
    double GM_E;       // Earth GM [m^3/s^2]
    double GM_M;       // Moon GM  [m^3/s^2]
    double R_E;        // Earth radius [m]
    double R_M;        // Moon radius [m]
    int n_maxM;        // Moon SH max degree/order
    int normalized;    // 1 normalized, 0 unnormalized
    double GM_S;       // Sun GM [m^3/s^2]
    int n_maxE;        // Earth SH max degree/order
    double R_S;        // Sun radius [m]
    double mass;       // Spacecraft mass [kg]
    double area;       // Spacecraft effective area [m^2]
    double eta;        // SRP scale factor [-]
};

struct SimulationConfig {
    int n_maxM = 10;                      // Moon SH max degree/order
    int n_maxE = 0;                       // Earth SH max degree/order
    string trajectory_file = "output_trajectory.txt";
    string output_file     = "measurements_output.txt";
};

struct GravityOutput {
    double U{};
    Vec3 dU{};
    Mat3 ddU{};
};

struct SRPOutput {
    Vec3 aSRP{};
    Mat3 daSRP_dr{};
    Vec3 daSRP_dEta{};
};

struct GSOutput {
    double g1{}, g2{}, g3{};
    double s1{}, s2{}, s3{}, s4{}, s5{}, s6{};
};

// ============================================================================
// Basic linear algebra utilities
// ============================================================================
double norm3(const Vec3& a) {
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

Vec3 add3(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec3 sub3(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 scale3(double s, const Vec3& a) {
    return {s*a[0], s*a[1], s*a[2]};
}

Vec3 mat3_vec(const Mat3& A, const Vec3& x) {
    Vec3 y{};
    for (int i = 0; i < 3; ++i) {
        y[i] = 0.0;
        for (int j = 0; j < 3; ++j) y[i] += A[i][j] * x[j];
    }
    return y;
}

Mat3 zeros3() {
    Mat3 A{};
    for (auto& row : A) row.fill(0.0);
    return A;
}

Mat3 eye3() {
    Mat3 I = zeros3();
    for (int i = 0; i < 3; ++i) I[i][i] = 1.0;
    return I;
}

Mat3 transpose3(const Mat3& A) {
    Mat3 AT = zeros3();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            AT[i][j] = A[j][i];
    return AT;
}

Mat3 mat3_mul(const Mat3& A, const Mat3& B) {
    Mat3 C = zeros3();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

Mat3 mat3_add(const Mat3& A, const Mat3& B) {
    Mat3 C = zeros3();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

Mat3 mat3_scale(double s, const Mat3& A) {
    Mat3 B = zeros3();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            B[i][j] = s * A[i][j];
    return B;
}

// ============================================================================
// SPICE helpers
// ============================================================================
void load_kernels() {
    furnsh_c("./kernels/kernels_LRO.tm");
}

Vec3 get_body_position_m(const string& target,
                         double et,
                         const string& frame,
                         const string& abcorr,
                         const string& observer) {
    SpiceDouble state[6];
    SpiceDouble lt;

    spkezr_c(target.c_str(), et, frame.c_str(), abcorr.c_str(), observer.c_str(), state, &lt);

    // SPICE returns km and km/s. Convert position to meters.
    return {state[0] * 1.0e3, state[1] * 1.0e3, state[2] * 1.0e3};
}

Mat3 get_rotation(const string& frame_from, const string& frame_to, double et) {
    SpiceDouble Rspice[3][3];
    pxform_c(frame_from.c_str(), frame_to.c_str(), et, Rspice);

    Mat3 R = zeros3();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            R[i][j] = Rspice[i][j];
    return R;
}

// ============================================================================
// Gravity-model helper functions: B and G/S coefficients
// ============================================================================
void getB_unnormalized(int n_max, double x, double y, double z, double r_n, double Re,
                       Matrix& b_real, Matrix& b_imag) {
    int N = n_max + 3;
    b_real.assign(N, vector<double>(N, 0.0));
    b_imag.assign(N, vector<double>(N, 0.0));

    for (int m = 0; m <= n_max + 2; ++m) {
        for (int n = m; n <= n_max + 2; ++n) {
            if (m == n) {
                if (m == 0) {
                    b_real[n][m] = Re / r_n;
                    b_imag[n][m] = 0.0;
                } else {
                    b_real[n][n] = (2.0*n - 1.0) * Re / r_n *
                                   ((x/r_n) * b_real[n-1][n-1] - (y/r_n) * b_imag[n-1][n-1]);
                    b_imag[n][n] = (2.0*n - 1.0) * Re / r_n *
                                   ((y/r_n) * b_real[n-1][n-1] + (x/r_n) * b_imag[n-1][n-1]);
                }
            } else {
                if (n >= 2) {
                    b_real[n][m] = (2.0*n - 1.0) / (n - m) * (Re*z) / (r_n*r_n) * b_real[n-1][m]
                                 - (n + m - 1.0) / (n - m) * pow(Re/r_n, 2) * b_real[n-2][m];
                    b_imag[n][m] = (2.0*n - 1.0) / (n - m) * (Re*z) / (r_n*r_n) * b_imag[n-1][m]
                                 - (n + m - 1.0) / (n - m) * pow(Re/r_n, 2) * b_imag[n-2][m];
                } else {
                    b_real[n][m] = (2.0*n - 1.0) / (n - m) * (Re*z) / (r_n*r_n) * b_real[n-1][m];
                    b_imag[n][m] = (2.0*n - 1.0) / (n - m) * (Re*z) / (r_n*r_n) * b_imag[n-1][m];
                }
            }
        }
    }
}

void getB_normalized(int n_max, double x, double y, double z, double r_n, double Re,
                     Matrix& b_real, Matrix& b_imag) {
    int N = n_max + 3;
    b_real.assign(N, vector<double>(N, 0.0));
    b_imag.assign(N, vector<double>(N, 0.0));

    for (int m = 0; m <= n_max + 2; ++m) {
        for (int n = m; n <= n_max + 2; ++n) {
            if (m == n) {
                if (m == 0) {
                    b_real[n][m] = Re / r_n;
                    b_imag[n][m] = 0.0;
                } else {
                    double delta_1_n = (n == 1) ? 1.0 : 0.0;
                    double fac = sqrt((1.0 + delta_1_n) * (2.0*n + 1.0) / (2.0*n)) * Re / r_n;
                    b_real[n][n] = fac * ((x/r_n) * b_real[n-1][n-1] - (y/r_n) * b_imag[n-1][n-1]);
                    b_imag[n][n] = fac * ((y/r_n) * b_real[n-1][n-1] + (x/r_n) * b_imag[n-1][n-1]);
                }
            } else {
                double fac1 = sqrt((4.0*n*n - 1.0) / (n*n - m*m)) * (Re*z) / (r_n*r_n);
                if (n >= 2) {
                    double fac2 = sqrt((2.0*n + 1.0) * ((n - 1.0)*(n - 1.0) - m*m) /
                                      ((2.0*n - 3.0) * (n*n - m*m))) * pow(Re/r_n, 2);
                    b_real[n][m] = fac1 * b_real[n-1][m] - fac2 * b_real[n-2][m];
                    b_imag[n][m] = fac1 * b_imag[n-1][m] - fac2 * b_imag[n-2][m];
                } else {
                    b_real[n][m] = fac1 * b_real[n-1][m];
                    b_imag[n][m] = fac1 * b_imag[n-1][m];
                }
            }
        }
    }
}

GSOutput getGS_unnormalized(int n, int m) {
    GSOutput out;
    double delta = (m == 0) ? 1.0 : 0.0;

    out.g1 = 0.5 * (1.0 + delta);
    out.g2 = 0.5 * (n - m + 2.0) * (n - m + 1.0);
    out.g3 = -(n - m + 1.0);

    out.s1 = 1.0;
    if (m == 0) out.s2 = (n + 2.0) * (n + 1.0);
    else if (m == 1) out.s2 = (n + 1.0) * n;
    else out.s2 = (n - m + 2.0) * (n - m + 1.0);

    out.s3 = (n - m + 4.0) * (n - m + 3.0) * (n - m + 2.0) * (n - m + 1.0);
    out.s4 = (m == 0) ? (n + 1.0) : (n - m + 1.0);
    out.s5 = (n - m + 3.0) * (n - m + 2.0) * (n - m + 1.0);
    out.s6 = (n - m + 2.0) * (n - m + 1.0);

    return out;
}

GSOutput getGS_normalized(int n, int m) {
    GSOutput out;
    double d0 = (m == 0) ? 1.0 : 0.0;
    double d1 = (m == 1) ? 1.0 : 0.0;
    double d2 = (m == 2) ? 1.0 : 0.0;

    out.g1 = 0.5 * (1.0 + d0) * sqrt((2.0 - d0) * (2.0*n + 1.0) *
                                     (n + m + 2.0) * (n + m + 1.0) /
                                     (2.0 * (2.0*n + 3.0)));
    out.g2 = 0.5 * sqrt((2.0 - d0) * (2.0*n + 1.0) *
                        (n - m + 2.0) * (n - m + 1.0) /
                        ((2.0 - d1) * (2.0*n + 3.0)));
    out.g3 = -sqrt((2.0*n + 1.0) * (n + m + 1.0) *
                   (n - m + 1.0) / (2.0*n + 3.0));

    out.s1 = sqrt(0.5 * (2.0 - d0) * (2.0*n + 1.0) *
                  (n + m + 4.0) * (n + m + 3.0) *
                  (n + m + 2.0) * (n + m + 1.0) / (2.0*n + 5.0));
    out.s2 = sqrt((2.0*n + 1.0) * (n + m + 2.0) * (n + m + 1.0) *
                  (n - m + 2.0) * (n - m + 1.0) / (2.0*n + 5.0));
    out.s3 = sqrt((2.0 - d0) * (2.0*n + 1.0) *
                  (n - m + 4.0) * (n - m + 3.0) *
                  (n - m + 2.0) * (n - m + 1.0) /
                  ((2.0 - d2) * (2.0*n + 5.0)));
    out.s4 = sqrt(0.5 * (2.0 - d0) * (2.0*n + 1.0) *
                  (n + m + 3.0) * (n + m + 2.0) *
                  (n + m + 1.0) * (n - m + 1.0) / (2.0*n + 5.0));
    out.s5 = sqrt((2.0 - d0) * (2.0*n + 1.0) *
                  (n + m + 1.0) * (n - m + 3.0) *
                  (n - m + 2.0) * (n - m + 1.0) /
                  ((2.0 - d1) * (2.0*n + 5.0)));
    out.s6 = sqrt((2.0*n + 1.0) * (n + m + 2.0) * (n + m + 1.0) *
                  (n - m + 2.0) * (n - m + 1.0) / (2.0*n + 5.0));

    return out;
}

// ============================================================================
// Force model functions translated from MATLAB
// ============================================================================
GravityOutput potentialGradient_nm(const Matrix& C_mat,
                                    const Matrix& S_mat,
                                    int n_max,
                                    const Vec3& r,
                                    double Re,
                                    double GM,
                                    int normalized) {
    GravityOutput out;
    out.ddU = zeros3();

    double r_n = norm3(r);
    double x = r[0], y = r[1], z = r[2];

    Matrix b_real, b_imag;
    if (normalized == 1) getB_normalized(n_max, x, y, z, r_n, Re, b_real, b_imag);
    else getB_unnormalized(n_max, x, y, z, r_n, Re, b_real, b_imag);

    double ddU_ddx = 0.0, ddU_ddz = 0.0;
    double ddU_dxdy = 0.0, ddU_dxdz = 0.0, ddU_dydz = 0.0;

    for (int n = 0; n <= n_max; ++n) {
        for (int m = 0; m <= n; ++m) {
            double C = C_mat[n][m];
            double S = S_mat[n][m];

            out.U += GM / Re * (b_real[n][m] * C + b_imag[n][m] * S);

            GSOutput gs = (normalized == 1) ? getGS_normalized(n, m) : getGS_unnormalized(n, m);

            if (m == 0) {
                out.dU[0] += GM / (Re*Re) * gs.g1 * (-C * b_real[n+1][m+1] - S * b_imag[n+1][m+1]);
                out.dU[1] += GM / (Re*Re) * gs.g1 * ( S * b_real[n+1][m+1] - C * b_imag[n+1][m+1]);
            } else {
                out.dU[0] += GM / (Re*Re) *
                    (gs.g1 * (-C*b_real[n+1][m+1] - S*b_imag[n+1][m+1]) +
                     gs.g2 * ( C*b_real[n+1][m-1] + S*b_imag[n+1][m-1]));
                out.dU[1] += GM / (Re*Re) *
                    (gs.g1 * ( S*b_real[n+1][m+1] - C*b_imag[n+1][m+1]) +
                     gs.g2 * ( S*b_real[n+1][m-1] - C*b_imag[n+1][m-1]));
            }

            out.dU[2] += GM / (Re*Re) * gs.g3 * (C*b_real[n+1][m] + S*b_imag[n+1][m]);

            // Hessian terms
            if (m == 0) {
                ddU_ddx += GM / pow(Re, 3) * 0.5 * C_mat[n][0] *
                           (gs.s1*b_real[n+2][2] - gs.s2*b_real[n+2][0]);
                ddU_dxdy += GM / pow(Re, 3) * 0.5 * C_mat[n][0] * gs.s1 * b_imag[n+2][2];
            } else if (m == 1) {
                ddU_ddx += GM / pow(Re, 3) * 0.25 *
                           (C_mat[n][1] * (gs.s1*b_real[n+2][3] - 3.0*gs.s2*b_real[n+2][1]) +
                            S_mat[n][1] * (gs.s1*b_imag[n+2][3] - gs.s2*b_imag[n+2][1]));
                ddU_dxdy += GM / pow(Re, 3) * 0.25 *
                            (-S_mat[n][1] * (gs.s1*b_real[n+2][3] + gs.s2*b_real[n+2][1]) +
                              C_mat[n][1] * (gs.s1*b_imag[n+2][3] - gs.s2*b_imag[n+2][1]));
            } else {
                ddU_ddx += GM / pow(Re, 3) * 0.25 *
                           (C * (gs.s1*b_real[n+2][m+2] - 2.0*gs.s2*b_real[n+2][m] + gs.s3*b_real[n+2][m-2]) +
                            S * (gs.s1*b_imag[n+2][m+2] - 2.0*gs.s2*b_imag[n+2][m] + gs.s3*b_imag[n+2][m-2]));
                ddU_dxdy += GM / pow(Re, 3) * 0.25 *
                            (-S * (gs.s1*b_real[n+2][m+2] - gs.s3*b_real[n+2][m-2]) +
                              C * (gs.s1*b_imag[n+2][m+2] - gs.s3*b_imag[n+2][m-2]));
            }

            if (m == 0) {
                ddU_dxdz += GM / pow(Re, 3) * gs.s4 * C_mat[n][0] * b_real[n+2][1];
                ddU_dydz += GM / pow(Re, 3) * gs.s4 * C_mat[n][0] * b_imag[n+2][1];
            } else {
                ddU_dxdz += GM / pow(Re, 3) * 0.5 *
                            (C * (gs.s4*b_real[n+2][m+1] - gs.s5*b_real[n+2][m-1]) +
                             S * (gs.s4*b_imag[n+2][m+1] - gs.s5*b_imag[n+2][m-1]));
                ddU_dydz += GM / pow(Re, 3) * 0.5 *
                            (-S * (gs.s4*b_real[n+2][m+1] + gs.s5*b_real[n+2][m-1]) +
                              C * (gs.s4*b_imag[n+2][m+1] + gs.s5*b_imag[n+2][m-1]));
            }

            ddU_ddz += GM / pow(Re, 3) * gs.s6 * (C*b_real[n+2][m] + S*b_imag[n+2][m]);
        }
    }

    double ddU_ddy = -ddU_ddx - ddU_ddz;
    out.ddU = {{{ddU_ddx,  ddU_dxdy, ddU_dxdz},
                {ddU_dxdy, ddU_ddy,  ddU_dydz},
                {ddU_dxdz, ddU_dydz, ddU_ddz}}};

    return out;
}

// Parallel version intended for the Moon spherical-harmonic gravity call only.
GravityOutput potentialGradient_nm_parallel_moon(const Matrix& C_mat,
                                    const Matrix& S_mat,
                                    int n_max,
                                    const Vec3& r,
                                    double Re,
                                    double GM,
                                    int normalized) {
    GravityOutput out;
    out.ddU = zeros3();

    double r_n = norm3(r);
    double x = r[0], y = r[1], z = r[2];

    Matrix b_real, b_imag;
    if (normalized == 1) getB_normalized(n_max, x, y, z, r_n, Re, b_real, b_imag);
    else getB_unnormalized(n_max, x, y, z, r_n, Re, b_real, b_imag);

    double U_acc = 0.0;
    double dUx = 0.0, dUy = 0.0, dUz = 0.0;
    double ddU_ddx = 0.0, ddU_ddz = 0.0;
    double ddU_dxdy = 0.0, ddU_dxdz = 0.0, ddU_dydz = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(none) schedule(static) \
    shared(C_mat, S_mat, b_real, b_imag, n_max, GM, Re, normalized) \
    reduction(+:U_acc,dUx,dUy,dUz,ddU_ddx,ddU_ddz,ddU_dxdy,ddU_dxdz,ddU_dydz)
#endif
    for (int n = 0; n <= n_max; ++n) {
        for (int m = 0; m <= n; ++m) {
            double C = C_mat[n][m];
            double S = S_mat[n][m];

            U_acc += GM / Re * (b_real[n][m] * C + b_imag[n][m] * S);

            GSOutput gs = (normalized == 1) ? getGS_normalized(n, m) : getGS_unnormalized(n, m);

            if (m == 0) {
                dUx += GM / (Re*Re) * gs.g1 * (-C * b_real[n+1][m+1] - S * b_imag[n+1][m+1]);
                dUy += GM / (Re*Re) * gs.g1 * ( S * b_real[n+1][m+1] - C * b_imag[n+1][m+1]);
            } else {
                dUx += GM / (Re*Re) *
                    (gs.g1 * (-C*b_real[n+1][m+1] - S*b_imag[n+1][m+1]) +
                     gs.g2 * ( C*b_real[n+1][m-1] + S*b_imag[n+1][m-1]));
                dUy += GM / (Re*Re) *
                    (gs.g1 * ( S*b_real[n+1][m+1] - C*b_imag[n+1][m+1]) +
                     gs.g2 * ( S*b_real[n+1][m-1] - C*b_imag[n+1][m-1]));
            }

            dUz += GM / (Re*Re) * gs.g3 * (C*b_real[n+1][m] + S*b_imag[n+1][m]);

            // Hessian terms
            if (m == 0) {
                ddU_ddx += GM / pow(Re, 3) * 0.5 * C_mat[n][0] *
                           (gs.s1*b_real[n+2][2] - gs.s2*b_real[n+2][0]);
                ddU_dxdy += GM / pow(Re, 3) * 0.5 * C_mat[n][0] * gs.s1 * b_imag[n+2][2];
            } else if (m == 1) {
                ddU_ddx += GM / pow(Re, 3) * 0.25 *
                           (C_mat[n][1] * (gs.s1*b_real[n+2][3] - 3.0*gs.s2*b_real[n+2][1]) +
                            S_mat[n][1] * (gs.s1*b_imag[n+2][3] - gs.s2*b_imag[n+2][1]));
                ddU_dxdy += GM / pow(Re, 3) * 0.25 *
                            (-S_mat[n][1] * (gs.s1*b_real[n+2][3] + gs.s2*b_real[n+2][1]) +
                              C_mat[n][1] * (gs.s1*b_imag[n+2][3] - gs.s2*b_imag[n+2][1]));
            } else {
                ddU_ddx += GM / pow(Re, 3) * 0.25 *
                           (C * (gs.s1*b_real[n+2][m+2] - 2.0*gs.s2*b_real[n+2][m] + gs.s3*b_real[n+2][m-2]) +
                            S * (gs.s1*b_imag[n+2][m+2] - 2.0*gs.s2*b_imag[n+2][m] + gs.s3*b_imag[n+2][m-2]));
                ddU_dxdy += GM / pow(Re, 3) * 0.25 *
                            (-S * (gs.s1*b_real[n+2][m+2] - gs.s3*b_real[n+2][m-2]) +
                              C * (gs.s1*b_imag[n+2][m+2] - gs.s3*b_imag[n+2][m-2]));
            }

            if (m == 0) {
                ddU_dxdz += GM / pow(Re, 3) * gs.s4 * C_mat[n][0] * b_real[n+2][1];
                ddU_dydz += GM / pow(Re, 3) * gs.s4 * C_mat[n][0] * b_imag[n+2][1];
            } else {
                ddU_dxdz += GM / pow(Re, 3) * 0.5 *
                            (C * (gs.s4*b_real[n+2][m+1] - gs.s5*b_real[n+2][m-1]) +
                             S * (gs.s4*b_imag[n+2][m+1] - gs.s5*b_imag[n+2][m-1]));
                ddU_dydz += GM / pow(Re, 3) * 0.5 *
                            (-S * (gs.s4*b_real[n+2][m+1] + gs.s5*b_real[n+2][m-1]) +
                              C * (gs.s4*b_imag[n+2][m+1] + gs.s5*b_imag[n+2][m-1]));
            }

            ddU_ddz += GM / pow(Re, 3) * gs.s6 * (C*b_real[n+2][m] + S*b_imag[n+2][m]);
        }
    }

    out.U = U_acc;
    out.dU = {dUx, dUy, dUz};

    double ddU_ddy = -ddU_ddx - ddU_ddz;
    out.ddU = {{{ddU_ddx,  ddU_dxdy, ddU_dxdz},
                {ddU_dxdy, ddU_ddy,  ddU_dydz},
                {ddU_dxdz, ddU_dydz, ddU_ddz}}};

    return out;
}

State state_add_scaled(const State& x, const State& k, double scale) {
    State y{};
    for (int i = 0; i < NTOTAL; ++i) y[i] = x[i] + scale*k[i];
    return y;
}

State state_linear_combination(const State& x,
                               double h,
                               const vector<pair<double, const State*>>& terms) {
    State y = x;
    for (int i = 0; i < NTOTAL; ++i) {
        for (const auto& term : terms) {
            y[i] += h * term.first * (*(term.second))[i];
        }
    }
    return y;
}

vector<double> linspace(double tmin, double tmax, int N) {
    vector<double> t(N);
    if (N == 1) {
        t[0] = tmin;
        return t;
    }
    double dt = (tmax - tmin) / static_cast<double>(N - 1);
    for (int k = 0; k < N; ++k) t[k] = tmin + k*dt;
    return t;
}

State initialize_state_with_stm(const Vec3& r0, const Vec3& v0) {
    State X{};
    X[0] = r0[0]; X[1] = r0[1]; X[2] = r0[2];
    X[3] = v0[0]; X[4] = v0[1]; X[5] = v0[2];

    // PHI0 = reshape(eye(6), [36, 1]) in MATLAB column-major ordering.
    for (int col = 0; col < NX; ++col)
        for (int row = 0; row < NX; ++row)
            X[NX + row + NX*col] = (row == col) ? 1.0 : 0.0;
    return X;
}

// ============================================================================
// Example coefficient and parameter initialization
// Replace this with your own SH coefficient loading routine.
// ============================================================================
Matrix make_zero_coeffs(int nmax) {
    return Matrix(nmax + 1, vector<double>(nmax + 1, 0.0));
}

void set_point_mass_coefficients(Matrix& C) {
    if (!C.empty() && !C[0].empty()) C[0][0] = 1.0;
}

Matrix read_csv_matrix(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open coefficient file: " + filename);
    }

    Matrix A;
    string line;

    while (getline(file, line)) {
        if (line.empty()) continue;

        vector<double> row;
        string value;
        stringstream ss(line);

        while (getline(ss, value, ',')) {
            if (!value.empty()) row.push_back(stod(value));
        }

        if (!row.empty()) A.push_back(row);
    }

    if (A.empty()) throw runtime_error("Coefficient file is empty: " + filename);

    size_t ncols = A[0].size();
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i].size() != ncols) {
            throw runtime_error("Coefficient file has inconsistent row lengths: " + filename);
        }
    }

    return A;
}

int get_available_max_degree(const Matrix& A) {
    return static_cast<int>(min(A.size(), A[0].size())) - 1;
}

Matrix truncate_coeff_matrix(const Matrix& A_full, int nmax_requested) {
    int nmax_available = get_available_max_degree(A_full);

    if (nmax_requested < 0) {
        throw runtime_error("Requested spherical-harmonic degree must be nonnegative.");
    }
    if (nmax_requested > nmax_available) {
        throw runtime_error("Requested spherical-harmonic degree is larger than coefficient file supports.");
    }

    Matrix A = make_zero_coeffs(nmax_requested);
    for (int n = 0; n <= nmax_requested; ++n) {
        for (int m = 0; m <= n; ++m) {
            A[n][m] = A_full[n][m];
        }
    }
    return A;
}

void load_moon_sh_coefficients(const string& C_file,
                               const string& S_file,
                               int nmax_requested,
                               Matrix& C_M,
                               Matrix& S_M) {
    Matrix C_full = read_csv_matrix(C_file);
    Matrix S_full = read_csv_matrix(S_file);

    int nmax_available = min(get_available_max_degree(C_full), get_available_max_degree(S_full));
    if (nmax_requested > nmax_available) {
        throw runtime_error("Requested Moon SH degree exceeds available C/S coefficient data.");
    }

    C_M = truncate_coeff_matrix(C_full, nmax_requested);
    S_M = truncate_coeff_matrix(S_full, nmax_requested);

    // The provided files already have C00 = 1 and S00 = 0, but enforce this.
    C_M[0][0] = 1.0;
    S_M[0][0] = 0.0;
}

string trim_copy(const string& s) {
    size_t first = 0;
    while (first < s.size() && isspace(static_cast<unsigned char>(s[first]))) ++first;

    size_t last = s.size();
    while (last > first && isspace(static_cast<unsigned char>(s[last - 1]))) --last;

    return s.substr(first, last - first);
}

string strip_inline_comment(const string& line) {
    size_t pos = line.find('#');
    if (pos == string::npos) return line;
    return line.substr(0, pos);
}

unordered_map<string, string> read_key_value_file(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open config file: " + filename);
    }

    unordered_map<string, string> kv;
    string line;
    int lineNumber = 0;

    while (getline(file, line)) {
        ++lineNumber;
        line = trim_copy(strip_inline_comment(line));
        if (line.empty()) continue;

        size_t eq = line.find('=');
        if (eq == string::npos) {
            throw runtime_error("Invalid config line " + to_string(lineNumber) +
                                ": expected key = value");
        }

        string key = trim_copy(line.substr(0, eq));
        string val = trim_copy(line.substr(eq + 1));

        if (key.empty() || val.empty()) {
            throw runtime_error("Invalid config line " + to_string(lineNumber) +
                                ": empty key or value");
        }

        kv[key] = val;
    }

    return kv;
}

string get_required_string(const unordered_map<string, string>& kv, const string& key) {
    auto it = kv.find(key);
    if (it == kv.end()) throw runtime_error("Missing required config key: " + key);
    return it->second;
}

double get_required_double(const unordered_map<string, string>& kv, const string& key) {
    return stod(get_required_string(kv, key));
}

int get_required_int(const unordered_map<string, string>& kv, const string& key) {
    return stoi(get_required_string(kv, key));
}

SimulationConfig read_simulation_config(const string& filename) {
    auto kv = read_key_value_file(filename);

    SimulationConfig cfg;
    cfg.n_maxM = get_required_int(kv, "n_maxM");
    cfg.n_maxE = get_required_int(kv, "n_maxE");
    cfg.trajectory_file = get_required_string(kv, "trajectory_file");
    cfg.output_file = get_required_string(kv, "output_file");

    if (cfg.n_maxM < 0) throw runtime_error("Config n_maxM must be nonnegative.");
    if (cfg.n_maxE < 0) throw runtime_error("Config n_maxE must be nonnegative.");
    if (cfg.output_file.empty()) throw runtime_error("Config output_file cannot be empty.");

    return cfg;
}

PlanetParams make_planet_params(const SimulationConfig& cfg) {
    PlanetParams p{};
    p.GM_E = 398600435507023;
    p.GM_M = 4.9028001224453001E12;
    p.R_E  = 6378140;
    p.R_M  = 1738000;
    p.n_maxM = cfg.n_maxM;
    p.normalized = 1;
    p.GM_S = 1.32712440041279e+20;
    p.n_maxE = cfg.n_maxE;
    p.R_S = 696000000;
    return p;
}

// ============================================================================
// Trajectory input and measurement simulation
// ============================================================================
struct TrajectoryData {
    vector<double> time;
    vector<Vec3> r;
    vector<Vec3> v;
};

TrajectoryData read_trajectory_txt(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open trajectory file: " + filename);
    }

    TrajectoryData traj;
    string line;

    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        stringstream ss(line);
        double t, x, y, z, vx, vy, vz;
        if (ss >> t >> x >> y >> z >> vx >> vy >> vz) {
            traj.time.push_back(t);
            traj.r.push_back({x, y, z});
            traj.v.push_back({vx, vy, vz});
        }
    }

    if (traj.time.empty()) {
        throw runtime_error("Trajectory file is empty or has the wrong format: " + filename);
    }

    return traj;
}

void save_true_measurements_txt(const string& filename,
                                const vector<double>& time,
                                const vector<array<double, 6>>& Y_mE) {
    ofstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open output measurement file: " + filename);
    }

    file << scientific << setprecision(numeric_limits<double>::max_digits10);
    file << "# time[s_past_J2000_TDB] Txx[mE] Txy[mE] Txz[mE] Tyy[mE] Tyz[mE] Tzz[mE]\n";

    for (size_t k = 0; k < time.size(); ++k) {
        file << time[k] << " "
             << Y_mE[k][0] << " "
             << Y_mE[k][1] << " "
             << Y_mE[k][2] << " "
             << Y_mE[k][3] << " "
             << Y_mE[k][4] << " "
             << Y_mE[k][5] << "\n";
    }
}

Mat3 gravity_tensor_in_j2000(const Matrix& C,
                             const Matrix& S,
                             int nmax,
                             const Vec3& r_j2000,
                             const Mat3& J2000_BODY,
                             double R,
                             double GM,
                             int normalized,
                             bool use_parallel_moon) {
    Mat3 BODY_J2000 = transpose3(J2000_BODY);
    Vec3 r_body = mat3_vec(BODY_J2000, r_j2000);

    GravityOutput g;
    if (use_parallel_moon) {
        g = potentialGradient_nm_parallel_moon(C, S, nmax, r_body, R, GM, normalized);
    } else {
        g = potentialGradient_nm(C, S, nmax, r_body, R, GM, normalized);
    }

    return mat3_mul(mat3_mul(J2000_BODY, g.ddU), BODY_J2000);
}

array<double, 6> pack_upper_trace_mE(const Mat3& T) {
    constexpr double SI_TO_mE = 1.0e12; // 1 mE = 1e-12 s^-2
    return {T[0][0] * SI_TO_mE,
            T[0][1] * SI_TO_mE,
            T[0][2] * SI_TO_mE,
            T[1][1] * SI_TO_mE,
            T[1][2] * SI_TO_mE,
            T[2][2] * SI_TO_mE};
}

vector<array<double, 6>> compute_measurements_true_from_trajectory(
    const TrajectoryData& traj,
    const PlanetParams& p,
    const vector<Matrix>& C_mat,
    const vector<Matrix>& S_mat) {

    const size_t Nt = traj.time.size();
    vector<array<double, 6>> Y(Nt);

    const Matrix& C_E = C_mat[0];
    const Matrix& S_E = S_mat[0];
    const Matrix& C_M = C_mat[1];
    const Matrix& S_M = S_mat[1];

    for (size_t k = 0; k < Nt; ++k) {
        const double et = traj.time[k];
        const Vec3& r_sc_moon_j2000 = traj.r[k];

        Mat3 J2000_EARTH = get_rotation("IAU_EARTH", "J2000", et);
        Mat3 J2000_MOON  = get_rotation("MOON_PA",   "J2000", et);

        Vec3 r_EM = get_body_position_m("EARTH", et, "J2000", "NONE", "MOON");
        Vec3 r_ES = sub3(r_sc_moon_j2000, r_EM);

        Mat3 T_M = gravity_tensor_in_j2000(C_M, S_M, p.n_maxM,
                                           r_sc_moon_j2000,
                                           J2000_MOON,
                                           p.R_M, p.GM_M,
                                           p.normalized,
                                           true);

        Mat3 T_E = gravity_tensor_in_j2000(C_E, S_E, p.n_maxE,
                                           r_ES,
                                           J2000_EARTH,
                                           p.R_E, p.GM_E,
                                           p.normalized,
                                           false);

        Mat3 T = mat3_add(T_M, T_E);
        Y[k] = pack_upper_trace_mE(T);
    }

    return Y;
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char* argv[]) {
    try {
        string configFile = (argc > 1) ? argv[1] : "simulation_config.txt";
        SimulationConfig cfg = read_simulation_config(configFile);

        //const string trajectoryFile = (argc > 1) ? argv[1] : "trajectory_output.txt";
        //const string outputFile     = (argc > 2) ? argv[2] : "true_measurements_output.txt";

        load_kernels();

#ifdef _OPENMP
        cout << "OpenMP enabled with up to " << omp_get_max_threads()
             << " threads for the Moon SH measurement summation." << endl;
#else
        cout << "OpenMP not enabled. Compile with OpenMP flags to parallelize the Moon SH summation." << endl;
#endif

        PlanetParams p = make_planet_params(cfg);

        vector<Matrix> C_mat(2), S_mat(2);
        C_mat[0] = make_zero_coeffs(p.n_maxE);
        S_mat[0] = make_zero_coeffs(p.n_maxE);
        set_point_mass_coefficients(C_mat[0]);

        load_moon_sh_coefficients("Cnm.txt", "Snm.txt", p.n_maxM, C_mat[1], S_mat[1]);
        cout << "Loaded Moon SH coefficients up to degree/order " << p.n_maxM << endl;

        cout << "Reading trajectory from: " << cfg.trajectory_file << endl;
        TrajectoryData traj = read_trajectory_txt(cfg.trajectory_file);
        cout << "Loaded " << traj.time.size() << " trajectory samples." << endl;

        cout << "Computing true gravity-gradient measurements..." << endl;
        vector<array<double, 6>> Y = compute_measurements_true_from_trajectory(traj, p, C_mat, S_mat);
        cout << "  DONE ..." << endl;

        save_true_measurements_txt(cfg.output_file, traj.time, Y);
        cout << "Saved measurements to: " << cfg.output_file << endl;

        cout << fixed << setprecision(9);
        cout << "First measurement [mE] = "
             << Y.front()[0] << " " << Y.front()[1] << " " << Y.front()[2] << " "
             << Y.front()[3] << " " << Y.front()[4] << " " << Y.front()[5] << endl;

        cout << "Last measurement [mE]  = "
             << Y.back()[0] << " " << Y.back()[1] << " " << Y.back()[2] << " "
             << Y.back()[3] << " " << Y.back()[4] << " " << Y.back()[5] << endl;

        kclear_c();
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        kclear_c();
        return 1;
    }

    return 0;
}
