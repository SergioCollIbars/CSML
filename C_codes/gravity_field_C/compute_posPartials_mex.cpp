/*
 * compute_posPartials_analytic_mex.cpp
 *
 * Analytical MEX implementation of compute_posPartials.m.
 * No finite differences are used. The code evaluates the third-order
 * gravity-potential derivatives from the same closed-form expressions used
 * in the MATLAB function, rotates the third-order tensor to the body frame,
 * and returns the 9-by-3 sensitivity matrix Hk.
 *
 * MATLAB usage:
 *   Hk = compute_posPartials_analytic_mex(n_max, normalized, C_mat, S_mat, ...
 *                                         R, GM, r, ACAF_ACI, ACAF_B)
 *
 * Compile:
 *   mex CXXFLAGS="$CXXFLAGS -std=c++17" compute_posPartials_analytic_mex.cpp
 */

#include "mex.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;
using Matrix = std::vector<std::vector<double>>;

struct Tvals {
    double T1, T2, T3, T4, T5, T6;
};

struct Nvals {
    double N1, N2, N3, N4, N5, N6, N7;
};

static double norm3(const Vec3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static Mat3 zeros3() {
    Mat3 A{};
    for (auto& row : A) row.fill(0.0);
    return A;
}

static Vec3 matVec3(const Mat3& A, const Vec3& x) {
    Vec3 y{0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            y[i] += A[i][j] * x[j];
    return y;
}

static Mat3 transpose3(const Mat3& A) {
    Mat3 AT = zeros3();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            AT[i][j] = A[j][i];
    return AT;
}

static Vec3 matlabVectorToVec3(const mxArray* A, const char* name) {
    if (!mxIsDouble(A) || mxIsComplex(A) || mxGetNumberOfElements(A) != 3) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:invalidInput",
                          "%s must be a real double vector with 3 elements.", name);
    }
    const double* p = mxGetPr(A);
    return {p[0], p[1], p[2]};
}

static Mat3 matlabMatrixToMat3(const mxArray* A, const char* name) {
    if (!mxIsDouble(A) || mxIsComplex(A) || mxGetM(A) != 3 || mxGetN(A) != 3) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:invalidInput",
                          "%s must be a real double 3-by-3 matrix.", name);
    }
    const double* p = mxGetPr(A);
    Mat3 M = zeros3();
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
            M[i][j] = p[i + j*3];
    return M;
}

static Matrix matlabMatrixToCoeff(const mxArray* A, int n_max, const char* name) {
    if (!mxIsDouble(A) || mxIsComplex(A)) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:invalidInput",
                          "%s must be a real double matrix.", name);
    }
    const mwSize nr = mxGetM(A);
    const mwSize nc = mxGetN(A);
    if (nr < static_cast<mwSize>(n_max + 1) || nc < static_cast<mwSize>(n_max + 1)) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:invalidSize",
                          "%s must be at least (n_max+1)-by-(n_max+1).", name);
    }

    const double* p = mxGetPr(A);
    Matrix M(n_max + 1, std::vector<double>(n_max + 1, 0.0));
    for (int n = 0; n <= n_max; ++n)
        for (int m = 0; m <= n_max; ++m)
            M[n][m] = p[n + m*nr];
    return M;
}

static void getB_unnormalized(int n_max, double x, double y, double z,
                              double r_n, double Re,
                              Matrix& b_real, Matrix& b_imag) {
    const int N = n_max + 4;
    b_real.assign(N, std::vector<double>(N, 0.0));
    b_imag.assign(N, std::vector<double>(N, 0.0));

    for (int m = 0; m <= n_max + 3; ++m) {
        for (int n = m; n <= n_max + 3; ++n) {
            if (m == n) {
                if (m == 0) {
                    b_real[n][m] = Re / r_n;
                    b_imag[n][m] = 0.0;
                } else {
                    b_real[n][n] = (2.0*n - 1.0) * Re / r_n *
                                   ((x/r_n)*b_real[n-1][n-1] - (y/r_n)*b_imag[n-1][n-1]);
                    b_imag[n][n] = (2.0*n - 1.0) * Re / r_n *
                                   ((y/r_n)*b_real[n-1][n-1] + (x/r_n)*b_imag[n-1][n-1]);
                }
            } else {
                if (n >= 2) {
                    b_real[n][m] = (2.0*n - 1.0)/(n - m) * (Re*z)/(r_n*r_n) * b_real[n-1][m]
                                 - (n + m - 1.0)/(n - m) * std::pow(Re/r_n, 2) * b_real[n-2][m];
                    b_imag[n][m] = (2.0*n - 1.0)/(n - m) * (Re*z)/(r_n*r_n) * b_imag[n-1][m]
                                 - (n + m - 1.0)/(n - m) * std::pow(Re/r_n, 2) * b_imag[n-2][m];
                } else {
                    b_real[n][m] = (2.0*n - 1.0)/(n - m) * (Re*z)/(r_n*r_n) * b_real[n-1][m];
                    b_imag[n][m] = (2.0*n - 1.0)/(n - m) * (Re*z)/(r_n*r_n) * b_imag[n-1][m];
                }
            }
        }
    }
}

static void getB_normalized(int n_max, double x, double y, double z,
                            double r_n, double Re,
                            Matrix& b_real, Matrix& b_imag) {
    const int N = n_max + 4;
    b_real.assign(N, std::vector<double>(N, 0.0));
    b_imag.assign(N, std::vector<double>(N, 0.0));

    for (int m = 0; m <= n_max + 3; ++m) {
        for (int n = m; n <= n_max + 3; ++n) {
            if (m == n) {
                if (m == 0) {
                    b_real[n][m] = Re / r_n;
                    b_imag[n][m] = 0.0;
                } else {
                    const double delta_1_n = (n == 1) ? 1.0 : 0.0;
                    const double fac = std::sqrt((1.0 + delta_1_n) * (2.0*n + 1.0) / (2.0*n)) * Re / r_n;
                    b_real[n][n] = fac * ((x/r_n)*b_real[n-1][n-1] - (y/r_n)*b_imag[n-1][n-1]);
                    b_imag[n][n] = fac * ((y/r_n)*b_real[n-1][n-1] + (x/r_n)*b_imag[n-1][n-1]);
                }
            } else {
                const double fac1 = std::sqrt((4.0*n*n - 1.0)/(n*n - m*m)) * (Re*z)/(r_n*r_n);
                if (n >= 2) {
                    const double fac2 = std::sqrt((2.0*n + 1.0) * ((n - 1.0)*(n - 1.0) - m*m) /
                                                  ((2.0*n - 3.0) * (n*n - m*m))) * std::pow(Re/r_n, 2);
                    b_real[n][m] = fac1*b_real[n-1][m] - fac2*b_real[n-2][m];
                    b_imag[n][m] = fac1*b_imag[n-1][m] - fac2*b_imag[n-2][m];
                } else {
                    b_real[n][m] = fac1*b_real[n-1][m];
                    b_imag[n][m] = fac1*b_imag[n-1][m];
                }
            }
        }
    }
}

static Tvals getT(int n, int m) {
    Tvals T;
    T.T1 = (n - m + 1.0);
    T.T2 = (n - m + 2.0)*(n - m + 1.0);
    T.T3 = (n - m + 3.0)*(n - m + 2.0)*(n - m + 1.0);
    T.T4 = (n - m + 4.0)*(n - m + 3.0)*(n - m + 2.0)*(n - m + 1.0);
    T.T5 = (n - m + 5.0)*(n - m + 4.0)*(n - m + 3.0)*(n - m + 2.0)*(n - m + 1.0);
    T.T6 = (n - m + 6.0)*(n - m + 5.0)*(n - m + 4.0)*(n - m + 3.0)*(n - m + 2.0)*(n - m + 1.0);
    return T;
}

static Nvals getN(int n, int m) {
    const double delta0 = (m == 0) ? 1.0 : 0.0;
    const double delta1 = (m == 1) ? 1.0 : 0.0;
    const double delta2 = (m == 2) ? 1.0 : 0.0;
    const double delta3 = (m == 3) ? 1.0 : 0.0;
    Nvals N;
    N.N1 = std::sqrt(((2.0-delta0)*(2.0*n + 1.0)*(n+m+6.0)*(n+m+5.0)*(n+m+4.0)*(n+m+3.0)*(n+m+2.0)*(n+m+1.0)) /
                     (2.0*(2.0*n + 7.0)));
    N.N2 = std::sqrt(((2.0-delta0)*(2.0*n + 1.0)*(n+m+5.0)*(n+m+4.0)*(n+m+3.0)*(n+m+2.0)*(n+m+1.0)) /
                     (2.0*(2.0*n + 7.0)*(n-m+1.0)));
    N.N3 = std::sqrt(((2.0-delta0)*(2.0*n + 1.0)*(n+m+4.0)*(n+m+3.0)*(n+m+2.0)*(n+m+1.0)) /
                     (2.0*(2.0*n + 7.0)*(n-m+2.0)*(n-m+1.0)));
    N.N4 = std::sqrt(((2.0*n + 1.0)*(n+m+3.0)*(n+m+2.0)*(n+m+1.0)) /
                     ((2.0*n + 7.0)*(n-m+3.0)*(n-m+2.0)*(n-m+1.0)));
    N.N5 = std::sqrt(((2.0-delta0)*(2.0*n + 1.0)*(n+m+2.0)*(n+m+1.0)) /
                     ((2.0-delta1)*(2.0*n + 7.0)*(n-m+4.0)*(n-m+3.0)*(n-m+2.0)*(n-m+1.0)));
    N.N6 = std::sqrt(((2.0-delta0)*(2.0*n + 1.0)*(n+m+1.0)) /
                     ((2.0-delta2)*(2.0*n + 7.0)*(n-m+5.0)*(n-m+4.0)*(n-m+3.0)*(n-m+2.0)*(n-m+1.0)));
    N.N7 = std::sqrt(((2.0-delta0)*(2.0*n + 1.0)) /
                     ((2.0-delta3)*(2.0*n + 7.0)*(n-m+6.0)*(n-m+5.0)*(n-m+4.0)*(n-m+3.0)*(n-m+2.0)*(n-m+1.0)));
    return N;
}

static void fillOutput(double* Hk, const double T[3][3][3]) {
    // MATLAB output is 9-by-3, column major.
    int row = 0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            for (int c = 0; c < 3; ++c) {
                Hk[row + c*9] = T[a][b][c];
            }
            ++row;
        }
    }
}

static void computeAnalytic(const Matrix& C_mat, const Matrix& S_mat,
                            int n_max, int normalized,
                            double R, double GM, const Vec3& r,
                            const Mat3& ACAF_ACI, const Mat3& ACAF_B,
                            double* Hk) {
    // MATLAB: ACI_ACAF = ACAF_ACI'; B_ACAF = ACAF_B'; r_pos = ACI_ACAF' * r = ACAF_ACI * r.
    const Mat3 B_ACAF = transpose3(ACAF_B);
    const Vec3 r_pos = matVec3(ACAF_ACI, r);

    const double r_n = norm3(r_pos);
    if (r_n <= 0.0) throw std::runtime_error("Position vector norm must be positive.");

    const double x = r_pos[0];
    const double y = r_pos[1];
    const double z = r_pos[2];

    Matrix b_real, b_imag;
    if (normalized == 1) getB_normalized(n_max, x, y, z, r_n, R, b_real, b_imag);
    else                 getB_unnormalized(n_max, x, y, z, r_n, R, b_real, b_imag);

    double dddU_dddx   = 0.0;
    double dddU_ddxdy  = 0.0;
    double dddU_ddxdz  = 0.0;
    double dddU_dxddy  = 0.0;
    double dddU_dxdydz = 0.0;
    double dddU_dddy   = 0.0;
    double dddU_ddydz  = 0.0;
    double dddU_dxddz  = 0.0;
    double dddU_dyddz  = 0.0;
    double dddU_dddz   = 0.0;

    const double R4 = R*R*R*R;
    const double GM_R4_8 = GM/(8.0*R4);
    const double GM_R4_4 = GM/(4.0*R4);
    const double GM_R4_2 = GM/(2.0*R4);
    const double GM_R4   = GM/R4;

    for (int n = 0; n <= n_max; ++n) {
        for (int m = 0; m <= n; ++m) {
            const double C = C_mat[n][m];
            const double S = S_mat[n][m];
            const Tvals T = getT(n, m);
            Nvals Nv;
            if (normalized == 1) Nv = getN(n, m);
            else Nv = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};

            // MATLAB b_real(N+3,M+3) maps to C++ b_real[n+3][m+3].
            auto br = [&](int dn, int dm) -> double { return b_real[n + dn][m + dm]; };
            auto bi = [&](int dn, int dm) -> double { return b_imag[n + dn][m + dm]; };

            if (m > 2) {
                dddU_dddx += GM_R4_8 * (
                    C*(-Nv.N1*br(3,3) + 3*T.T2*Nv.N3*br(3,1) - 3*T.T4*Nv.N5*br(3,-1) + T.T6*Nv.N7*br(3,-3))
                  + S*(-Nv.N1*bi(3,3) + 3*T.T2*Nv.N3*bi(3,1) - 3*T.T4*Nv.N5*bi(3,-1) + T.T6*Nv.N7*bi(3,-3)) );
                dddU_ddxdy += GM_R4_8 * (
                    S*(T.T6*Nv.N7*br(3,-3) - T.T4*Nv.N5*br(3,-1) - T.T2*Nv.N3*br(3,1) + Nv.N1*br(3,3))
                  - C*(T.T6*Nv.N7*bi(3,-3) - T.T4*Nv.N5*bi(3,-1) - T.T2*Nv.N3*bi(3,1) + Nv.N1*bi(3,3)) );
                dddU_ddxdz += GM_R4_4 * (
                    C*(-T.T1*Nv.N2*br(3,2) + 2*T.T3*Nv.N4*br(3,0) - T.T5*Nv.N6*br(3,-2))
                  + S*(-T.T1*Nv.N2*bi(3,2) + 2*T.T3*Nv.N4*bi(3,0) - T.T5*Nv.N6*bi(3,-2)) );
                dddU_dxddy += GM_R4_8 * (
                    C*(-T.T6*Nv.N7*br(3,-3) - T.T4*Nv.N5*br(3,-1) + T.T2*Nv.N3*br(3,1) + Nv.N1*br(3,3))
                  + S*(-T.T6*Nv.N7*bi(3,-3) - T.T4*Nv.N5*bi(3,-1) + T.T2*Nv.N3*bi(3,1) + Nv.N1*bi(3,3)) );
                dddU_dxdydz += GM_R4_4 * (
                    S*(T.T1*Nv.N2*br(3,2) - T.T5*Nv.N6*br(3,-2))
                  - C*(T.T1*Nv.N2*bi(3,2) - T.T5*Nv.N6*bi(3,-2)) );
                dddU_dddy += GM_R4_8 * (
                    S*(-Nv.N1*br(3,3) - 3*T.T2*Nv.N3*br(3,1) - 3*T.T4*Nv.N5*br(3,-1) - T.T6*Nv.N7*br(3,-3))
                  - C*(-Nv.N1*bi(3,3) - 3*T.T2*Nv.N3*bi(3,1) - 3*T.T4*Nv.N5*bi(3,-1) - T.T6*Nv.N7*bi(3,-3)) );
                dddU_ddydz += GM_R4_4 * (
                    C*(T.T1*Nv.N2*br(3,2) + 2*T.T3*Nv.N4*br(3,0) + T.T5*Nv.N6*br(3,-2))
                  + S*(T.T1*Nv.N2*bi(3,2) + 2*T.T3*Nv.N4*bi(3,0) + T.T5*Nv.N6*bi(3,-2)) );
                dddU_dxddz += GM_R4_2 * (
                    C*(-T.T2*Nv.N3*br(3,1) + T.T4*Nv.N5*br(3,-1))
                  + S*(-T.T2*Nv.N3*bi(3,1) + T.T4*Nv.N5*bi(3,-1)) );
                dddU_dyddz += GM_R4_2 * (
                    S*(T.T2*Nv.N3*br(3,1) + T.T4*Nv.N5*br(3,-1))
                  - C*(T.T2*Nv.N3*bi(3,1) + T.T4*Nv.N5*bi(3,-1)) );
            } else if (m == 2) {
                const double A = (n+2.0)*(n+1.0)*n*(n-1.0);
                dddU_dddx += GM_R4_8 * (C*(-Nv.N1*br(3,3)+3*T.T2*Nv.N3*br(3,1)-3*T.T4*Nv.N5*br(3,-1)-A*Nv.N5*br(3,-1)) + S*(-Nv.N1*bi(3,3)+3*T.T2*Nv.N3*bi(3,1)-3*T.T4*Nv.N5*bi(3,-1)+A*Nv.N5*bi(3,-1)));
                dddU_ddxdy += GM_R4_8 * (S*(-A*Nv.N5*br(3,-1)-T.T4*Nv.N5*br(3,-1)-T.T2*Nv.N3*br(3,1)+Nv.N1*br(3,3)) - C*(A*Nv.N5*bi(3,-1)-T.T4*Nv.N5*bi(3,-1)-T.T2*Nv.N3*bi(3,1)+Nv.N1*bi(3,3)));
                dddU_ddxdz += GM_R4_4 * (C*(-T.T1*Nv.N2*br(3,2)+2*T.T3*Nv.N4*br(3,0)-T.T5*Nv.N6*br(3,-2)) + S*(-T.T1*Nv.N2*bi(3,2)+2*T.T3*Nv.N4*bi(3,0)-T.T5*Nv.N6*bi(3,-2)));
                dddU_dxddy += GM_R4_8 * (C*(A*Nv.N5*br(3,-1)-T.T4*Nv.N5*br(3,-1)+T.T2*Nv.N3*br(3,1)+Nv.N1*br(3,3)) + S*(-A*Nv.N5*bi(3,-1)-T.T4*Nv.N5*bi(3,-1)+T.T2*Nv.N3*bi(3,1)+Nv.N1*bi(3,3)));
                dddU_dxdydz += GM_R4_4 * (S*(T.T1*Nv.N2*br(3,2)-T.T5*Nv.N6*br(3,-2)) - C*(T.T1*Nv.N2*bi(3,2)-T.T5*Nv.N6*bi(3,-2)));
                dddU_dddy += GM_R4_8 * (S*(-Nv.N1*br(3,3)-3*T.T2*Nv.N3*br(3,1)-3*T.T4*Nv.N5*br(3,-1)+A*Nv.N5*br(3,-1)) - C*(-Nv.N1*bi(3,3)-3*T.T2*Nv.N3*bi(3,1)-3*T.T4*Nv.N5*bi(3,-1)-A*Nv.N5*bi(3,-1)));
                dddU_ddydz += GM_R4_4 * (C*(T.T1*Nv.N2*br(3,2)+2*T.T3*Nv.N4*br(3,0)+T.T5*Nv.N6*br(3,-2)) + S*(T.T1*Nv.N2*bi(3,2)+2*T.T3*Nv.N4*bi(3,0)+T.T5*Nv.N6*bi(3,-2)));
                dddU_dxddz += GM_R4_2 * (C*(-T.T2*Nv.N3*br(3,1)+T.T4*Nv.N5*br(3,-1)) + S*(-T.T2*Nv.N3*bi(3,1)+T.T4*Nv.N5*bi(3,-1)));
                dddU_dyddz += GM_R4_2 * (S*(T.T2*Nv.N3*br(3,1)+T.T4*Nv.N5*br(3,-1)) - C*(T.T2*Nv.N3*bi(3,1)+T.T4*Nv.N5*bi(3,-1)));
            } else if (m == 1) {
                const double A1 = (n+1.0)*n;
                const double A2 = (n+2.0)*(n+1.0)*n;
                dddU_dddx += GM_R4_8 * (C*(-Nv.N1*br(3,3)+3*T.T2*Nv.N3*br(3,1)-3*T.T4*Nv.N5*br(3,-1)+A1*Nv.N3*br(3,1)) + S*(-Nv.N1*bi(3,3)+3*T.T2*Nv.N3*bi(3,1)-3*T.T4*Nv.N5*bi(3,-1)-A1*Nv.N3*bi(3,1)));
                dddU_ddxdy += GM_R4_8 * (S*(A1*Nv.N3*br(3,1)-T.T4*Nv.N5*br(3,-1)-T.T2*Nv.N3*br(3,1)+Nv.N1*br(3,3)) - C*(-A1*Nv.N3*bi(3,1)-T.T4*Nv.N5*bi(3,-1)-T.T2*Nv.N3*bi(3,1)+Nv.N1*bi(3,3)));
                dddU_ddxdz += GM_R4_4 * (C*(-T.T1*Nv.N2*br(3,2)+2*T.T3*Nv.N4*br(3,0)+A2*Nv.N4*br(3,0)) + S*(-T.T1*Nv.N2*bi(3,2)+2*T.T3*Nv.N4*bi(3,0)-A2*Nv.N4*bi(3,0)));
                dddU_dxddy += GM_R4_8 * (C*(-A1*Nv.N3*br(3,1)-T.T4*Nv.N5*br(3,-1)+T.T2*Nv.N3*br(3,1)+Nv.N1*br(3,3)) + S*(A1*Nv.N3*bi(3,1)-T.T4*Nv.N5*bi(3,-1)+T.T2*Nv.N3*bi(3,1)+Nv.N1*bi(3,3)));
                dddU_dxdydz += GM_R4_4 * (S*(T.T1*Nv.N2*br(3,2)+A2*Nv.N4*br(3,0)) - C*(T.T1*Nv.N2*bi(3,2)-A2*Nv.N4*bi(3,0)));
                dddU_dddy += GM_R4_8 * (S*(-Nv.N1*br(3,3)-3*T.T2*Nv.N3*br(3,1)-3*T.T4*Nv.N5*br(3,-1)-A1*Nv.N3*br(3,1)) - C*(-Nv.N1*bi(3,3)-3*T.T2*Nv.N3*bi(3,1)-3*T.T4*Nv.N5*bi(3,-1)+A1*Nv.N3*bi(3,1)));
                dddU_ddydz += GM_R4_4 * (C*(T.T1*Nv.N2*br(3,2)+2*T.T3*Nv.N4*br(3,0)-A2*Nv.N4*br(3,0)) + S*(T.T1*Nv.N2*bi(3,2)+2*T.T3*Nv.N4*bi(3,0)+A2*Nv.N4*bi(3,0)));
                dddU_dxddz += GM_R4_2 * (C*(-T.T2*Nv.N3*br(3,1)+T.T4*Nv.N5*br(3,-1)) + S*(-T.T2*Nv.N3*bi(3,1)+T.T4*Nv.N5*bi(3,-1)));
                dddU_dyddz += GM_R4_2 * (S*(T.T2*Nv.N3*br(3,1)+T.T4*Nv.N5*br(3,-1)) - C*(T.T2*Nv.N3*bi(3,1)+T.T4*Nv.N5*bi(3,-1)));
            } else { // m == 0
                const double A1 = (n+2.0)*(n+1.0);
                const double A2 = (n+1.0);
                dddU_dddx += GM_R4_8 * (C*(-Nv.N1*br(3,3)+3*T.T2*Nv.N3*br(3,1)+3*A1*Nv.N3*br(3,1)-Nv.N1*br(3,3)) + S*(-Nv.N1*bi(3,3)+3*T.T2*Nv.N3*bi(3,1)-3*A1*Nv.N3*bi(3,1)+Nv.N1*bi(3,3)));
                dddU_ddxdy += GM_R4_8 * (S*(-Nv.N1*br(3,3)+A1*Nv.N3*br(3,1)-T.T2*Nv.N3*br(3,1)+Nv.N1*br(3,3)) - C*(Nv.N1*bi(3,3)-A1*Nv.N3*bi(3,1)-T.T2*Nv.N3*bi(3,1)+Nv.N1*bi(3,3)));
                dddU_ddxdz += GM_R4_4 * (C*(-T.T1*Nv.N2*br(3,2)+2*T.T3*Nv.N4*br(3,0)-A2*Nv.N2*br(3,2)) + S*(-T.T1*Nv.N2*bi(3,2)+2*T.T3*Nv.N4*bi(3,0)+A2*Nv.N2*bi(3,2)));
                dddU_dxddy += GM_R4_8 * (C*(Nv.N1*br(3,3)+A1*Nv.N3*br(3,1)+T.T2*Nv.N3*br(3,1)+Nv.N1*br(3,3)) + S*(-Nv.N1*bi(3,3)-A1*Nv.N3*bi(3,1)+T.T2*Nv.N3*bi(3,1)+Nv.N1*bi(3,3)));
                dddU_dxdydz += GM_R4_4 * (S*(T.T1*Nv.N2*br(3,2)-A2*Nv.N2*br(3,2)) - C*(T.T1*Nv.N2*bi(3,2)+A2*Nv.N2*bi(3,2)));
                dddU_dddy += GM_R4_8 * (S*(-Nv.N1*br(3,3)-3*T.T2*Nv.N3*br(3,1)+3*A1*Nv.N3*br(3,1)+Nv.N1*br(3,3)) - C*(-Nv.N1*bi(3,3)-3*T.T2*Nv.N3*bi(3,1)-3*A1*Nv.N3*bi(3,1)-Nv.N1*bi(3,3)));
                dddU_ddydz += GM_R4_4 * (C*(T.T1*Nv.N2*br(3,2)+2*T.T3*Nv.N4*br(3,0)+A2*Nv.N2*br(3,2)) + S*(T.T1*Nv.N2*bi(3,2)+2*T.T3*Nv.N4*bi(3,0)-A2*Nv.N2*bi(3,2)));
                dddU_dxddz += GM_R4_2 * (C*(-T.T2*Nv.N3*br(3,1)-A1*Nv.N3*br(3,1)) + S*(-T.T2*Nv.N3*bi(3,1)+A1*Nv.N3*bi(3,1)));
                dddU_dyddz += GM_R4_2 * (S*(T.T2*Nv.N3*br(3,1)-A1*Nv.N3*br(3,1)) - C*(T.T2*Nv.N3*bi(3,1)+A1*Nv.N3*bi(3,1)));
            }

            dddU_dddz -= GM_R4 * T.T3 * (C*Nv.N4*br(3,0) + S*Nv.N4*bi(3,0));
        }
    }

    double T_acaf[3][3][3]{};
    T_acaf[0][0][0] = dddU_dddx;   T_acaf[0][0][1] = dddU_ddxdy;  T_acaf[0][0][2] = dddU_ddxdz;
    T_acaf[0][1][0] = dddU_ddxdy;  T_acaf[0][1][1] = dddU_dxddy;  T_acaf[0][1][2] = dddU_dxdydz;
    T_acaf[0][2][0] = dddU_ddxdz;  T_acaf[0][2][1] = dddU_dxdydz; T_acaf[0][2][2] = dddU_dxddz;

    T_acaf[1][0][0] = dddU_ddxdy;  T_acaf[1][0][1] = dddU_dxddy;  T_acaf[1][0][2] = dddU_dxdydz;
    T_acaf[1][1][0] = dddU_dxddy;  T_acaf[1][1][1] = dddU_dddy;   T_acaf[1][1][2] = dddU_ddydz;
    T_acaf[1][2][0] = dddU_dxdydz; T_acaf[1][2][1] = dddU_ddydz;  T_acaf[1][2][2] = dddU_dyddz;

    T_acaf[2][0][0] = dddU_ddxdz;  T_acaf[2][0][1] = dddU_dxdydz; T_acaf[2][0][2] = dddU_dxddz;
    T_acaf[2][1][0] = dddU_dxdydz; T_acaf[2][1][1] = dddU_ddydz;  T_acaf[2][1][2] = dddU_dyddz;
    T_acaf[2][2][0] = dddU_dxddz;  T_acaf[2][2][1] = dddU_dyddz;  T_acaf[2][2][2] = dddU_dddz;

    double T_body[3][3][3]{};
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            for (int c = 0; c < 3; ++c)
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        for (int k = 0; k < 3; ++k)
                            T_body[a][b][c] += B_ACAF[a][i]*B_ACAF[b][j]*B_ACAF[c][k]*T_acaf[i][j][k];

    fillOutput(Hk, T_body);
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 9) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:nrhs",
            "Usage: Hk = compute_posPartials_analytic_mex(n_max, normalized, C_mat, S_mat, R, GM, r, ACAF_ACI, ACAF_B)");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:nlhs", "Too many output arguments.");
    }

    const int n_max = static_cast<int>(mxGetScalar(prhs[0]));
    const int normalized = static_cast<int>(mxGetScalar(prhs[1]));
    if (n_max < 0) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:invalidInput", "n_max must be nonnegative.");
    }
    if (!(normalized == 0 || normalized == 1)) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:invalidInput", "normalized must be 0 or 1.");
    }

    Matrix C_mat = matlabMatrixToCoeff(prhs[2], n_max, "C_mat");
    Matrix S_mat = matlabMatrixToCoeff(prhs[3], n_max, "S_mat");
    const double R  = mxGetScalar(prhs[4]);
    const double GM = mxGetScalar(prhs[5]);
    const Vec3 r = matlabVectorToVec3(prhs[6], "r");
    const Mat3 ACAF_ACI = matlabMatrixToMat3(prhs[7], "ACAF_ACI");
    const Mat3 ACAF_B   = matlabMatrixToMat3(prhs[8], "ACAF_B");

    if (R <= 0.0 || GM == 0.0) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:invalidInput", "R must be positive and GM must be nonzero.");
    }

    plhs[0] = mxCreateDoubleMatrix(9, 3, mxREAL);
    double* Hk = mxGetPr(plhs[0]);

    try {
        computeAnalytic(C_mat, S_mat, n_max, normalized, R, GM, r, ACAF_ACI, ACAF_B, Hk);
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("compute_posPartials_analytic_mex:runtime", "%s", e.what());
    }
}
