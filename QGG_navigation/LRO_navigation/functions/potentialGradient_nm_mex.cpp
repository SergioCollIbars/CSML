/*
 * potentialGradient_nm_mex.cpp
 *
 * MATLAB MEX wrapper for the potentialGradient_nm function
 * extracted/adapted from potentialGradinet_nm.m
 *
 * Usage:
 *   [U, dU, ddU] = potentialGradient_nm_mex(C, S, n_max, r, Re, GM, normalized)
 *
 * Inputs:
 *   C, S       : coefficient matrices with C(n+1,m+1), S(n+1,m+1)
 *                for degree n and order m. Must contain at least rows/cols
 *                up to n_max+1 in MATLAB indexing.
 *   n_max      : maximum degree/order, scalar integer
 *   r          : 3x1 or 1x3 position vector [m], in body-fixed frame
 *   Re         : reference radius [m]
 *   GM         : gravitational parameter [m^3/s^2]
 *   normalized : 1 for normalized coefficients, 0 for unnormalized
 *
 * Outputs:
 *   U    : potential [m^2/s^2]
 *   dU   : 3x1 gradient/acceleration-like vector [m/s^2]
 *   ddU  : 3x3 Hessian / gravity-gradient tensor [1/s^2]
 *
 * Compile in MATLAB:
 *   mex CXXFLAGS="$CXXFLAGS -std=c++17" potentialGradient_nm_mex.cpp
 */

#include "mex.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;
using Matrix = std::vector<std::vector<double>>;

struct GravityOutput {
    double U{};
    Vec3 dU{};
    Mat3 ddU{};
};

struct GSOutput {
    double g1{}, g2{}, g3{};
    double s1{}, s2{}, s3{}, s4{}, s5{}, s6{};
};

static double norm3(const Vec3& a) {
    return std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

static Mat3 zeros3() {
    Mat3 A{};
    for (auto& row : A) row.fill(0.0);
    return A;
}

static void getB_unnormalized(int n_max, double x, double y, double z,
                              double r_n, double Re,
                              Matrix& b_real, Matrix& b_imag) {
    int N = n_max + 3;
    b_real.assign(N, std::vector<double>(N, 0.0));
    b_imag.assign(N, std::vector<double>(N, 0.0));

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
                                 - (n + m - 1.0) / (n - m) * std::pow(Re/r_n, 2) * b_real[n-2][m];
                    b_imag[n][m] = (2.0*n - 1.0) / (n - m) * (Re*z) / (r_n*r_n) * b_imag[n-1][m]
                                 - (n + m - 1.0) / (n - m) * std::pow(Re/r_n, 2) * b_imag[n-2][m];
                } else {
                    b_real[n][m] = (2.0*n - 1.0) / (n - m) * (Re*z) / (r_n*r_n) * b_real[n-1][m];
                    b_imag[n][m] = (2.0*n - 1.0) / (n - m) * (Re*z) / (r_n*r_n) * b_imag[n-1][m];
                }
            }
        }
    }
}

static void getB_normalized(int n_max, double x, double y, double z,
                            double r_n, double Re,
                            Matrix& b_real, Matrix& b_imag) {
    int N = n_max + 3;
    b_real.assign(N, std::vector<double>(N, 0.0));
    b_imag.assign(N, std::vector<double>(N, 0.0));

    for (int m = 0; m <= n_max + 2; ++m) {
        for (int n = m; n <= n_max + 2; ++n) {
            if (m == n) {
                if (m == 0) {
                    b_real[n][m] = Re / r_n;
                    b_imag[n][m] = 0.0;
                } else {
                    double delta_1_n = (n == 1) ? 1.0 : 0.0;
                    double fac = std::sqrt((1.0 + delta_1_n) * (2.0*n + 1.0) / (2.0*n)) * Re / r_n;
                    b_real[n][n] = fac * ((x/r_n) * b_real[n-1][n-1] - (y/r_n) * b_imag[n-1][n-1]);
                    b_imag[n][n] = fac * ((y/r_n) * b_real[n-1][n-1] + (x/r_n) * b_imag[n-1][n-1]);
                }
            } else {
                double fac1 = std::sqrt((4.0*n*n - 1.0) / (n*n - m*m)) * (Re*z) / (r_n*r_n);
                if (n >= 2) {
                    double fac2 = std::sqrt((2.0*n + 1.0) * ((n - 1.0)*(n - 1.0) - m*m) /
                                      ((2.0*n - 3.0) * (n*n - m*m))) * std::pow(Re/r_n, 2);
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

static GSOutput getGS_unnormalized(int n, int m) {
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

static GSOutput getGS_normalized(int n, int m) {
    GSOutput out;
    double d0 = (m == 0) ? 1.0 : 0.0;
    double d1 = (m == 1) ? 1.0 : 0.0;
    double d2 = (m == 2) ? 1.0 : 0.0;

    out.g1 = 0.5 * (1.0 + d0) * std::sqrt((2.0 - d0) * (2.0*n + 1.0) *
                                     (n + m + 2.0) * (n + m + 1.0) /
                                     (2.0 * (2.0*n + 3.0)));
    out.g2 = 0.5 * std::sqrt((2.0 - d0) * (2.0*n + 1.0) *
                        (n - m + 2.0) * (n - m + 1.0) /
                        ((2.0 - d1) * (2.0*n + 3.0)));
    out.g3 = -std::sqrt((2.0*n + 1.0) * (n + m + 1.0) *
                   (n - m + 1.0) / (2.0*n + 3.0));

    out.s1 = std::sqrt(0.5 * (2.0 - d0) * (2.0*n + 1.0) *
                  (n + m + 4.0) * (n + m + 3.0) *
                  (n + m + 2.0) * (n + m + 1.0) / (2.0*n + 5.0));
    out.s2 = std::sqrt((2.0*n + 1.0) * (n + m + 2.0) * (n + m + 1.0) *
                  (n - m + 2.0) * (n - m + 1.0) / (2.0*n + 5.0));
    out.s3 = std::sqrt((2.0 - d0) * (2.0*n + 1.0) *
                  (n - m + 4.0) * (n - m + 3.0) *
                  (n - m + 2.0) * (n - m + 1.0) /
                  ((2.0 - d2) * (2.0*n + 5.0)));
    out.s4 = std::sqrt(0.5 * (2.0 - d0) * (2.0*n + 1.0) *
                  (n + m + 3.0) * (n + m + 2.0) *
                  (n + m + 1.0) * (n - m + 1.0) / (2.0*n + 5.0));
    out.s5 = std::sqrt((2.0 - d0) * (2.0*n + 1.0) *
                  (n + m + 1.0) * (n - m + 3.0) *
                  (n - m + 2.0) * (n - m + 1.0) /
                  ((2.0 - d1) * (2.0*n + 5.0)));
    out.s6 = std::sqrt((2.0*n + 1.0) * (n + m + 2.0) * (n + m + 1.0) *
                  (n - m + 2.0) * (n - m + 1.0) / (2.0*n + 5.0));

    return out;
}

static GravityOutput potentialGradient_nm(const Matrix& C_mat,
                                          const Matrix& S_mat,
                                          int n_max,
                                          const Vec3& r,
                                          double Re,
                                          double GM,
                                          int normalized) {
    GravityOutput out;
    out.ddU = zeros3();

    double r_n = norm3(r);
    if (r_n <= 0.0) {
        throw std::runtime_error("Position vector norm must be positive.");
    }

    double x = r[0], y = r[1], z = r[2];

    Matrix b_real, b_imag;
    if (normalized == 1) getB_normalized(n_max, x, y, z, r_n, Re, b_real, b_imag);
    else getB_unnormalized(n_max, x, y, z, r_n, Re, b_real, b_imag);

    double ddU_ddx = 0.0, ddU_ddz = 0.0;
    double ddU_dxdy = 0.0, ddU_dxdz = 0.0, ddU_dydz = 0.0;
    double Re2 = Re * Re;
    double Re3 = Re * Re * Re;

    for (int n = 0; n <= n_max; ++n) {
        for (int m = 0; m <= n; ++m) {
            double C = C_mat[n][m];
            double S = S_mat[n][m];

            out.U += GM / Re * (b_real[n][m] * C + b_imag[n][m] * S);

            GSOutput gs = (normalized == 1) ? getGS_normalized(n, m) : getGS_unnormalized(n, m);

            if (m == 0) {
                out.dU[0] += GM / Re2 * gs.g1 * (-C * b_real[n+1][m+1] - S * b_imag[n+1][m+1]);
                out.dU[1] += GM / Re2 * gs.g1 * ( S * b_real[n+1][m+1] - C * b_imag[n+1][m+1]);
            } else {
                out.dU[0] += GM / Re2 *
                    (gs.g1 * (-C*b_real[n+1][m+1] - S*b_imag[n+1][m+1]) +
                     gs.g2 * ( C*b_real[n+1][m-1] + S*b_imag[n+1][m-1]));
                out.dU[1] += GM / Re2 *
                    (gs.g1 * ( S*b_real[n+1][m+1] - C*b_imag[n+1][m+1]) +
                     gs.g2 * ( S*b_real[n+1][m-1] - C*b_imag[n+1][m-1]));
            }

            out.dU[2] += GM / Re2 * gs.g3 * (C*b_real[n+1][m] + S*b_imag[n+1][m]);

            if (m == 0) {
                ddU_ddx += GM / Re3 * 0.5 * C_mat[n][0] *
                           (gs.s1*b_real[n+2][2] - gs.s2*b_real[n+2][0]);
                ddU_dxdy += GM / Re3 * 0.5 * C_mat[n][0] * gs.s1 * b_imag[n+2][2];
            } else if (m == 1) {
                ddU_ddx += GM / Re3 * 0.25 *
                           (C_mat[n][1] * (gs.s1*b_real[n+2][3] - 3.0*gs.s2*b_real[n+2][1]) +
                            S_mat[n][1] * (gs.s1*b_imag[n+2][3] - gs.s2*b_imag[n+2][1]));
                ddU_dxdy += GM / Re3 * 0.25 *
                            (-S_mat[n][1] * (gs.s1*b_real[n+2][3] + gs.s2*b_real[n+2][1]) +
                              C_mat[n][1] * (gs.s1*b_imag[n+2][3] - gs.s2*b_imag[n+2][1]));
            } else {
                ddU_ddx += GM / Re3 * 0.25 *
                           (C * (gs.s1*b_real[n+2][m+2] - 2.0*gs.s2*b_real[n+2][m] + gs.s3*b_real[n+2][m-2]) +
                            S * (gs.s1*b_imag[n+2][m+2] - 2.0*gs.s2*b_imag[n+2][m] + gs.s3*b_imag[n+2][m-2]));
                ddU_dxdy += GM / Re3 * 0.25 *
                            (-S * (gs.s1*b_real[n+2][m+2] - gs.s3*b_real[n+2][m-2]) +
                              C * (gs.s1*b_imag[n+2][m+2] - gs.s3*b_imag[n+2][m-2]));
            }

            if (m == 0) {
                ddU_dxdz += GM / Re3 * gs.s4 * C_mat[n][0] * b_real[n+2][1];
                ddU_dydz += GM / Re3 * gs.s4 * C_mat[n][0] * b_imag[n+2][1];
            } else {
                ddU_dxdz += GM / Re3 * 0.5 *
                            (C * (gs.s4*b_real[n+2][m+1] - gs.s5*b_real[n+2][m-1]) +
                             S * (gs.s4*b_imag[n+2][m+1] - gs.s5*b_imag[n+2][m-1]));
                ddU_dydz += GM / Re3 * 0.5 *
                            (-S * (gs.s4*b_real[n+2][m+1] + gs.s5*b_real[n+2][m-1]) +
                              C * (gs.s4*b_imag[n+2][m+1] + gs.s5*b_imag[n+2][m-1]));
            }

            ddU_ddz += GM / Re3 * gs.s6 * (C*b_real[n+2][m] + S*b_imag[n+2][m]);
        }
    }

    double ddU_ddy = -ddU_ddx - ddU_ddz;
    out.ddU = {{{ddU_ddx,  ddU_dxdy, ddU_dxdz},
                {ddU_dxdy, ddU_ddy,  ddU_dydz},
                {ddU_dxdz, ddU_dydz, ddU_ddz}}};

    return out;
}

static Matrix matlabMatrixToCoefficientMatrix(const mxArray* A, int n_max, const char* name) {
    if (!mxIsDouble(A) || mxIsComplex(A)) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:invalidInput", "%s must be a real double matrix.", name);
    }

    mwSize nr = mxGetM(A);
    mwSize nc = mxGetN(A);
    if (nr < static_cast<mwSize>(n_max + 1) || nc < static_cast<mwSize>(n_max + 1)) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:invalidSize",
                          "%s must be at least (n_max+1)-by-(n_max+1).", name);
    }

    const double* p = mxGetPr(A);
    Matrix M(n_max + 1, std::vector<double>(n_max + 1, 0.0));

    for (int n = 0; n <= n_max; ++n) {
        for (int m = 0; m <= n_max; ++m) {
            M[n][m] = p[n + m * nr];
        }
    }
    return M;
}

static Vec3 matlabVectorToVec3(const mxArray* A) {
    if (!mxIsDouble(A) || mxIsComplex(A) || mxGetNumberOfElements(A) != 3) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:invalidInput", "r must be a real double vector with 3 elements.");
    }
    const double* p = mxGetPr(A);
    return {p[0], p[1], p[2]};
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 7) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:nrhs",
            "Usage: [U,dU,ddU] = potentialGradient_nm_mex(C,S,n_max,r,Re,GM,normalized)");
    }
    if (nlhs > 3) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:nlhs", "Too many output arguments.");
    }

    int n_max = static_cast<int>(mxGetScalar(prhs[2]));
    if (n_max < 0) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:invalidInput", "n_max must be nonnegative.");
    }

    Matrix C = matlabMatrixToCoefficientMatrix(prhs[0], n_max, "C");
    Matrix S = matlabMatrixToCoefficientMatrix(prhs[1], n_max, "S");
    Vec3 r = matlabVectorToVec3(prhs[3]);
    double Re = mxGetScalar(prhs[4]);
    double GM = mxGetScalar(prhs[5]);
    int normalized = static_cast<int>(mxGetScalar(prhs[6]));

    if (!(normalized == 0 || normalized == 1)) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:invalidInput", "normalized must be 0 or 1.");
    }
    if (Re <= 0.0 || GM == 0.0) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:invalidInput", "Re must be positive and GM must be nonzero.");
    }

    try {
        GravityOutput out = potentialGradient_nm(C, S, n_max, r, Re, GM, normalized);

        plhs[0] = mxCreateDoubleScalar(out.U);

        if (nlhs >= 2) {
            plhs[1] = mxCreateDoubleMatrix(3, 1, mxREAL);
            double* dU = mxGetPr(plhs[1]);
            for (int i = 0; i < 3; ++i) dU[i] = out.dU[i];
        }

        if (nlhs >= 3) {
            plhs[2] = mxCreateDoubleMatrix(3, 3, mxREAL);
            double* ddU = mxGetPr(plhs[2]);
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    ddU[i + j*3] = out.ddU[i][j];
                }
            }
        }
    } catch (const std::exception& e) {
        mexErrMsgIdAndTxt("potentialGradient_nm_mex:runtime", "%s", e.what());
    }
}
