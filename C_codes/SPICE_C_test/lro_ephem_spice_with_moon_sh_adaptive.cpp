#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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
    // Edit these paths/names for your local kernel set.
    //furnsh_c("naif0012.tls");
    //furnsh_c("de421.bsp");
    furnsh_c("/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm");
    // Required if you use MOON_PA. Uncomment and provide the actual files.
    // furnsh_c("pck00010.tpc");
    // furnsh_c("moon_080317.tf");
    // furnsh_c("moon_pa_de421_1900-2050.bpc");
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

Vec3 get_body_velocity_m(const string& target,
                         double et,
                         const string& frame,
                         const string& abcorr,
                         const string& observer) {
    SpiceDouble state[6];
    SpiceDouble lt;

    spkezr_c(target.c_str(), et, frame.c_str(), abcorr.c_str(), observer.c_str(), state, &lt);

    // SPICE returns km and km/s. Convert position to meters.
    return {state[3] * 1.0e3, state[4] * 1.0e3, state[5] * 1.0e3};
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

SRPOutput SRP(const Vec3& rs, double eta, double m, double A) {
    SRPOutput out;

    double Gamma = 5.67e-8;       // [kg/(s^3 K^4)] Stefan-Boltzmann constant in SI units
    double c     = 2.99792458e8;  // [m/s]
    double Rs    = 6.96e8;        // [m]
    double Ts    = 5778.0;        // [K]
    double Cr    = 1.0;

    double rsn = norm3(rs);
    double P = (Gamma * Rs*Rs * pow(Ts, 4)) / (c * m);
    double factor = eta * P * Cr * A / pow(rsn, 3);

    out.aSRP = scale3(factor, rs);
    out.daSRP_dEta = scale3(P * Cr * A / pow(rsn, 3), rs);

    double pf = eta * Gamma * pow(Ts, 4) * Rs*Rs * Cr * A / (m * c * pow(rsn, 3));
    Mat3 I = eye3();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out.daSRP_dr[i][j] = pf * (I[i][j] - 3.0 * rs[i] * rs[j] / (rsn*rsn));
        }
    }

    return out;
}

double shadow_function(double R_sun, double R_planet, const Vec3& r_sun, const Vec3& r_sc) {
    double d_sc = norm3(r_sc);                  // planet -> SC
    double d_sun = norm3(sub3(r_sun, r_sc));    // SC -> Sun

    double theta_p = asin(R_planet / d_sc);
    double theta_s = asin(R_sun / d_sun);

    Vec3 sc_to_sun = sub3(r_sun, r_sc);
    double B = -(r_sc[0]*sc_to_sun[0] + r_sc[1]*sc_to_sun[1] + r_sc[2]*sc_to_sun[2]) / (d_sc*d_sun);
    B = max(-1.0, min(1.0, B));
    double theta_ps = acos(B);

    double S = 0.0;
    double J1 = fabs(theta_p - theta_s);
    double J2 = theta_p + theta_s;

    if (theta_p + theta_s <= theta_ps) {
        S = 0.0;
    } else if (theta_p - theta_s > theta_ps) {
        S = M_PI * theta_s * theta_s;
    } else if (theta_s - theta_p > theta_ps) {
        S = M_PI * theta_p * theta_p;
    } else if (J1 < theta_ps && theta_ps < J2) {
        double arg1 = (theta_s*theta_s + theta_ps*theta_ps - theta_p*theta_p) / (2.0*theta_s*theta_ps);
        double arg2 = (theta_p*theta_p + theta_ps*theta_ps - theta_s*theta_s) / (2.0*theta_p*theta_ps);
        arg1 = max(-1.0, min(1.0, arg1));
        arg2 = max(-1.0, min(1.0, arg2));

        double CAF = acos(arg1);
        double CBD = acos(arg2);

        double SAFC = 0.5 * CAF * theta_s * theta_s;
        double SAEC = 0.5 * (theta_s * sin(CAF)) * (theta_s * cos(CAF));
        double SBDC = 0.5 * CBD * theta_p * theta_p;
        double SBEC = 0.5 * (theta_p * sin(CBD)) * (theta_p * cos(CBD));

        S = 2.0*(SAFC - SAEC) + 2.0*(SBDC - SBEC);
    } else {
        cerr << "Warning: no shadow condition fulfilled. Setting F = 1.\n";
        return 1.0;
    }

    return 1.0 - S / (M_PI * theta_s * theta_s);
}

Mat6 compute_jacobian(const Mat3& T) {
    Mat6 J{};
    for (auto& row : J) row.fill(0.0);

    for (int i = 0; i < 3; ++i) {
        J[i][i + 3] = 1.0;
        for (int j = 0; j < 3; ++j) J[i + 3][j] = T[i][j];
    }
    return J;
}

// ============================================================================
// EOM and STM propagation
// ============================================================================
State EOM_LRO_EPHEM(double t,
                    const State& x,
                    const PlanetParams& p,
                    const vector<Matrix>& C_mat,
                    const vector<Matrix>& S_mat) {
    State dx{};

    Vec3 r_sc = {x[0], x[1], x[2]};
    Vec3 v_sc = {x[3], x[4], x[5]};

    const string ref = "J2000";
    const string abcorr = "NONE";
    const string observer = "MOON";

    Vec3 Epos = get_body_position_m("EARTH", t, ref, abcorr, observer);
    Vec3 Mpos = get_body_position_m("MOON",  t, ref, abcorr, observer);
    Vec3 Spos = get_body_position_m("SUN",   t, ref, abcorr, observer);

    Vec3 r1 = sub3(r_sc, Epos);   // SC - Earth
    Vec3 r2 = sub3(r_sc, Mpos);   // SC - Moon
    Vec3 r3 = sub3(r_sc, Spos);   // SC - Sun

    Vec3 r_ME = sub3(Mpos, Epos); // Moon - Earth
    Vec3 r_MS = sub3(Mpos, Spos); // Moon - Sun

    Mat3 J2000_EARTH = get_rotation("IAU_EARTH", "J2000", t);
    Mat3 J2000_MOON  = get_rotation("MOON_PA",   "J2000", t);

    Mat3 EARTH_J2000 = transpose3(J2000_EARTH);
    Mat3 MOON_J2000  = transpose3(J2000_MOON);

    double F = shadow_function(p.R_S, p.R_M, Spos, r2);
    SRPOutput srp = SRP(r3, p.eta, p.mass, p.area);

    const Matrix& C_E = C_mat[0];
    const Matrix& S_E = S_mat[0];
    const Matrix& C_M = C_mat[1];
    const Matrix& S_M = S_mat[1];

    GravityOutput Earth = potentialGradient_nm(C_E, S_E, p.n_maxE,
                                               mat3_vec(EARTH_J2000, r1),
                                               p.R_E, p.GM_E, p.normalized);

    GravityOutput Earth_tide = potentialGradient_nm(C_E, S_E, p.n_maxE,
                                                    mat3_vec(EARTH_J2000, r_ME),
                                                    p.R_E, p.GM_E, p.normalized);

    GravityOutput Moon = potentialGradient_nm(C_M, S_M, p.n_maxM,
                                              mat3_vec(MOON_J2000, r2),
                                              p.R_M, p.GM_M, p.normalized);

    // MATLAB code used Cmat_M/Smat_M and n=0 for the Sun point-mass-like term.
    GravityOutput Sun = potentialGradient_nm(C_M, S_M, 0, r3, p.R_S, p.GM_S, p.normalized);

    Vec3 dU1 = mat3_vec(J2000_EARTH, Earth.dU);
    Vec3 dU2 = mat3_vec(J2000_MOON,  Moon.dU);
    Vec3 dU3 = Sun.dU;

    Mat3 ddU1 = mat3_mul(mat3_mul(J2000_EARTH, Earth.ddU), EARTH_J2000);
    Mat3 ddU2 = mat3_mul(mat3_mul(J2000_MOON,  Moon.ddU),  MOON_J2000);
    Mat3 ddU3 = Sun.ddU;

    Vec3 a_tidal_E = scale3(-1.0, mat3_vec(J2000_EARTH, Earth_tide.dU));
    Vec3 a_tidal_S = scale3(p.GM_S / pow(norm3(r_MS), 3), r_MS);
    Vec3 aSRP = scale3(F, srp.aSRP);

    Vec3 acc = add3(add3(add3(dU2, dU1), dU3), add3(add3(a_tidal_E, a_tidal_S), aSRP));

    Mat3 T = mat3_add(mat3_add(mat3_add(ddU2, ddU1), ddU3), mat3_scale(F, srp.daSRP_dr));
    Mat6 J = compute_jacobian(T);

    // Unpack Phi using MATLAB column-major ordering.
    double Phi[NX][NX]{};
    for (int col = 0; col < NX; ++col)
        for (int row = 0; row < NX; ++row)
            Phi[row][col] = x[NX + row + NX*col];

    double Phi_dot[NX][NX]{};
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NX; ++j)
            for (int k = 0; k < NX; ++k)
                Phi_dot[i][j] += J[i][k] * Phi[k][j];

    dx[0] = v_sc[0];
    dx[1] = v_sc[1];
    dx[2] = v_sc[2];
    dx[3] = acc[0];
    dx[4] = acc[1];
    dx[5] = acc[2];

    for (int col = 0; col < NX; ++col)
        for (int row = 0; row < NX; ++row)
            dx[NX + row + NX*col] = Phi_dot[row][col];

    return dx;
}

// ============================================================================
// Adaptive integrator utilities.
//
// This section replaces the previous fixed-step RK4 with an adaptive
// Dormand-Prince RK4(5) method. It is similar in spirit to MATLAB ode45:
// the solution is advanced with a fifth-order formula, while an embedded
// fourth-order formula provides a local truncation-error estimate used to
// accept/reject steps and change the internal step size.
//
// The output grid is still controlled by the vector `time` in main(). The
// adaptive integrator simply takes as many internal steps as needed between
// time[k] and time[k+1].
// ============================================================================
struct IntegratorOptions {
    double relTol = 1.0e-13;
    double absTol = 1.0e-13;
    double initialStep = 0.1;   // Initial internal step-size guess [s]
    double minStep = 1.0e-12;   // Minimum allowed internal step [s]
    double maxStep = 60.0;      // Maximum allowed internal step [s]
    int maxAttempts = 100000;   // Prevents infinite loops if tolerance is too strict
};

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

struct RK45StepResult {
    State x5{};
    State x4{};
    double errNorm = 0.0;
};

RK45StepResult rk45_dormand_prince_step(double t,
                                         const State& x,
                                         double h,
                                         const PlanetParams& p,
                                         const vector<Matrix>& C_mat,
                                         const vector<Matrix>& S_mat,
                                         const IntegratorOptions& opts) {
    // Dormand-Prince coefficients, same family used by MATLAB ode45.
    State k1 = EOM_LRO_EPHEM(t, x, p, C_mat, S_mat);

    State x2 = state_linear_combination(x, h, {{1.0/5.0, &k1}});
    State k2 = EOM_LRO_EPHEM(t + h*(1.0/5.0), x2, p, C_mat, S_mat);

    State x3 = state_linear_combination(x, h, {{3.0/40.0, &k1},
                                               {9.0/40.0, &k2}});
    State k3 = EOM_LRO_EPHEM(t + h*(3.0/10.0), x3, p, C_mat, S_mat);

    State x4s = state_linear_combination(x, h, {{44.0/45.0, &k1},
                                                {-56.0/15.0, &k2},
                                                {32.0/9.0, &k3}});
    State k4 = EOM_LRO_EPHEM(t + h*(4.0/5.0), x4s, p, C_mat, S_mat);

    State x5s = state_linear_combination(x, h, {{19372.0/6561.0, &k1},
                                                {-25360.0/2187.0, &k2},
                                                {64448.0/6561.0, &k3},
                                                {-212.0/729.0, &k4}});
    State k5 = EOM_LRO_EPHEM(t + h*(8.0/9.0), x5s, p, C_mat, S_mat);

    State x6 = state_linear_combination(x, h, {{9017.0/3168.0, &k1},
                                               {-355.0/33.0, &k2},
                                               {46732.0/5247.0, &k3},
                                               {49.0/176.0, &k4},
                                               {-5103.0/18656.0, &k5}});
    State k6 = EOM_LRO_EPHEM(t + h, x6, p, C_mat, S_mat);

    State x7 = state_linear_combination(x, h, {{35.0/384.0, &k1},
                                               {500.0/1113.0, &k3},
                                               {125.0/192.0, &k4},
                                               {-2187.0/6784.0, &k5},
                                               {11.0/84.0, &k6}});
    State k7 = EOM_LRO_EPHEM(t + h, x7, p, C_mat, S_mat);

    // Fifth-order solution.
    State x_5 = x7;

    // Embedded fourth-order solution.
    State x_4 = state_linear_combination(x, h, {{5179.0/57600.0, &k1},
                                                {7571.0/16695.0, &k3},
                                                {393.0/640.0, &k4},
                                                {-92097.0/339200.0, &k5},
                                                {187.0/2100.0, &k6},
                                                {1.0/40.0, &k7}});

    // Weighted RMS error norm, similar to standard adaptive ODE solvers.
    double sumSq = 0.0;
    for (int i = 0; i < NTOTAL; ++i) {
        double scale = opts.absTol + opts.relTol * max(fabs(x[i]), fabs(x_5[i]));
        double e = (x_5[i] - x_4[i]) / scale;
        sumSq += e * e;
    }

    RK45StepResult result;
    result.x5 = x_5;
    result.x4 = x_4;
    result.errNorm = sqrt(sumSq / static_cast<double>(NTOTAL));
    return result;
}

State integrate_adaptive_rk45(double t0,
                              const State& X0,
                              double tf,
                              const PlanetParams& p,
                              const vector<Matrix>& C_mat,
                              const vector<Matrix>& S_mat,
                              const IntegratorOptions& opts) {
    State X = X0;
    double t = t0;

    if (tf == t0) return X;

    const double direction = (tf > t0) ? 1.0 : -1.0;
    double h = direction * min(opts.maxStep, max(opts.minStep, fabs(opts.initialStep)));

    int attempts = 0;
    while (direction * (tf - t) > 0.0) {
        if (++attempts > opts.maxAttempts) {
            throw runtime_error("Adaptive RK45 exceeded maxAttempts. Try relaxing tolerances or increasing maxAttempts.");
        }

        // Do not step beyond the requested output time.
        if (direction * (t + h - tf) > 0.0) {
            h = tf - t;
        }

        RK45StepResult step = rk45_dormand_prince_step(t, X, h, p, C_mat, S_mat, opts);

        // Accept if normalized error is <= 1. Also accept exactly-zero error.
        if (step.errNorm <= 1.0 || step.errNorm == 0.0) {
            X = step.x5;
            t += h;
        }

        // Step-size update. For a 4/5 embedded method, use exponent 1/5.
        double safety = 0.9;
        double factor;
        if (step.errNorm == 0.0) {
            factor = 5.0;
        } else {
            factor = safety * pow(1.0 / step.errNorm, 0.2);
            factor = min(5.0, max(0.1, factor));
        }

        double hNew = h * factor;
        double absHNew = fabs(hNew);

        if (absHNew < opts.minStep) {
            if (step.errNorm > 1.0) {
                throw runtime_error("Adaptive RK45 step size fell below minStep before satisfying tolerance.");
            }
            absHNew = opts.minStep;
        }

        absHNew = min(absHNew, opts.maxStep);
        h = direction * absHNew;
    }

    return X;
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

PlanetParams make_planet_params() {
    PlanetParams p{};
    p.GM_E = 398600435507023;
    p.GM_M = 4.9028001224453001E12;
    p.R_E  = 6378140;
    p.R_M  = 1738000;
    p.n_maxM = 100;  // Change this to truncate the Moon SH field.
    p.normalized = 1;
    p.GM_S = 1.32712440041279e+20;
    p.n_maxE = 0;
    p.R_S = 696000000;
    p.mass = 220.0;
    p.area = 2.8;
    p.eta = 0;
    return p;
}

// ============================================================================
// Main
// ============================================================================
int main() {
    try {
        load_kernels();

        PlanetParams p = make_planet_params();

        vector<Matrix> C_mat(2), S_mat(2);

        // Earth is kept as point mass here. Replace with Earth SH loading if needed.
        C_mat[0] = make_zero_coeffs(p.n_maxE);
        S_mat[0] = make_zero_coeffs(p.n_maxE);
        set_point_mass_coefficients(C_mat[0]);

        // Moon spherical-harmonic coefficients.
        // Cnm.txt and Snm.txt contain coefficients up to degree/order 10.
        // To use a lower degree, only change p.n_maxM in make_planet_params().
        load_moon_sh_coefficients("Cnm.txt", "Snm.txt", p.n_maxM, C_mat[1], S_mat[1]);
        cout << "Loaded Moon SH coefficients up to degree/order " << p.n_maxM << endl;

        SpiceDouble tmin, tmax;
        str2et_c("2015 MAR 20 00:00:00 UTC", &tmin);
        str2et_c("2015 MAR 20 06:00:00 UTC", &tmax);

        double fs = 1.0;  // instrumentParams_GG(1,5) [Hz]
        int N = static_cast<int>((tmax - tmin) * fs) + 1;
        vector<double> time = linspace(tmin, tmax, N);

        Vec3 r0 = get_body_position_m("LRO", tmin, "J2000", "NONE", "MOON");
        Vec3 v0 = get_body_velocity_m("LRO", tmin, "J2000", "NONE", "MOON");
        State X = initialize_state_with_stm(r0, v0);

        vector<State> state;
        state.reserve(N);

        // Adaptive RK4(5) settings. These are intentionally similar to your MATLAB
        // tolerances. If the integrator takes extremely small steps, relax to 1e-11
        // or 1e-12 first, then tighten after validating the force model.
        IntegratorOptions odeOpts;
        odeOpts.relTol = 1.0e-10;
        odeOpts.absTol = 1.0e-10;
        odeOpts.initialStep = 0.1;
        odeOpts.minStep = 1.0e-10;
        odeOpts.maxStep = 30.0;

        for (int k = 0; k < N; ++k) {
            state.push_back(X);

            if (k < N - 1) {
                X = integrate_adaptive_rk45(time[k], X, time[k+1],
                                            p, C_mat, S_mat, odeOpts);
            }
        }

        cout << "  DONE ..." << endl;
        cout << "Simulated time [s] = " << (tmax - tmin) << endl;
        cout << fixed << setprecision(9);
        cout << "Initial r [m]   = " << state[0][0] << "  "
                                  << state[0][1] << "  "
                                  << state[0][2] << endl;
        cout << "Initial v [m/s] = " << state[0][3] << "  "
                                  << state[0][4] << "  "
                                  << state[0][5] << endl;
        cout << "Final r [m]   = " << state.back()[0] << "  "
                                  << state.back()[1] << "  "
                                  << state.back()[2] << endl;
        cout << "Final v [m/s] = " << state.back()[3] << "  "
                                  << state.back()[4] << "  "
                                  << state.back()[5] << endl;

        kclear_c();
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        kclear_c();
        return 1;
    }

    return 0;
}
