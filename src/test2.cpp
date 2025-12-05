// joukowsky_panel_with_constant_gamma_chord.cpp
// Build with: g++ -std=c++17 -O2 joukowsky_panel_with_constant_gamma_chord.cpp -o joukowsky_panel
// Example run: ./joukowsky_panel

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Vec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;

constexpr double PI = 3.14159265358979323846;
constexpr double RadToDegree = 180.0 / PI;
constexpr double DegreeToRad = PI / 180.0;

// ----------------------- Linear algebra (from previous tidy up) -----------------------
Mat inverseLU_v2(const Mat& A)
{
    int n = (int)A.size();
    if (n == 0) return {};
    for (auto const& r : A) if ((int)r.size() != n) throw std::invalid_argument("inverseLU_v2: A must be square");

    Mat L(n, Vec(n, 0.0)), U(n, Vec(n, 0.0));
    for (int i = 0; i < n; ++i) L[i][i] = 1.0;

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) sum += L[i][k] * U[k][j];
            U[i][j] = A[i][j] - sum;
        }
        if (std::abs(U[i][i]) < 1e-14) {
            std::ostringstream ss; ss << "zero pivot at U[" << i << "][" << i << "]";
            throw std::runtime_error(ss.str());
        }
        for (int j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < i; ++k) sum += L[j][k] * U[k][i];
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    }

    // invL (forward)
    Mat invL(n, Vec(n, 0.0));
    for (int i = 0; i < n; ++i) {
        invL[i][i] = 1.0;
        for (int j = 0; j < i; ++j) {
            double s = 0.0;
            for (int k = j; k < i; ++k) s += L[i][k] * invL[k][j];
            invL[i][j] = -s;
        }
    }

    // invU (backwards)
    Mat invU(n, Vec(n, 0.0));
    for (int i = n - 1; i >= 0; --i) {
        invU[i][i] = 1.0 / U[i][i];
        for (int j = i + 1; j < n; ++j) {
            double s = 0.0;
            for (int k = i + 1; k <= j; ++k) s += U[i][k] * invU[k][j];
            invU[i][j] = -s / U[i][i];
        }
    }

    Mat invA(n, Vec(n, 0.0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < n; ++k) s += invU[i][k] * invL[k][j];
            invA[i][j] = s;
        }

    return invA;
}

Vec multiplyMatVec(const Mat& A, const Vec& x)
{
    int n = (int)A.size();
    if (n == 0) return {};
    if ((int)x.size() != n) throw std::invalid_argument("multiplyMatVec: size mismatch");
    Vec y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        if ((int)A[i].size() != n) throw std::invalid_argument("multiplyMatVec: A must be square");
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += A[i][j] * x[j];
        y[i] = s;
    }
    return y;
}

// ----------------------- Joukowsky transform geometry -----------------------
std::complex<double> joukowskyTransform(std::complex<double> z, double shape)
{
    // shape acts as 'a' parameter (radius)
    return z + std::pow(shape,2)/z;
}

// Generate Joukowsky aerofoil (close loop) from a circle displaced by (disX, disY)
// backFat controls scaling of the sampling radius; resolution = number of points on the contour
// Returns vectors px, py of size = resolution (closed loop)
void genJoukowsky(double disX, double disY, double backFat, int resolution, std::vector<double>& px, std::vector<double>& py)
{
    px.clear(); py.clear();
    px.reserve(resolution);
    py.reserve(resolution);

    // circle radius and mapping parameter
    double r = 1.0;          // base circle radius
    double shape = 1.0;      // 'a' used inside transform (keeps units consistent)
    std::complex<double> center(disX, disY);

    // choose a sampling radius around center
    double extraR = std::sqrt(std::pow(r - std::abs(center.real()), 2.0) + std::pow(center.imag(), 2.0));

    // sample angle around circle (we will produce a closed loop of resolution points)
    for (int i = 0; i < resolution; ++i) {
        double theta = -PI + 2.0 * PI * double(i) / double(resolution); // -pi..pi
        std::complex<double> z = std::complex<double>(extraR * std::cos(theta) / backFat + center.real(),
                                                      extraR * std::sin(theta) / backFat + center.imag());
        std::complex<double> air = joukowskyTransform(z, shape);
        // scale & shift to make chord roughly order 1
        double X = 0.25 * (air.real() + 2.0);
        double Y = 0.25 * air.imag();
        px.push_back(X);
        py.push_back(Y);
    }
}

// ----------------------- Panel geometry and influence integrals -----------------------
void buildPanelsAndInfluences(const std::vector<double>& pointX, const std::vector<double>& pointY,
                              Mat& I, Mat& J, Mat& K, Mat& L,
                              std::vector<double>& midX, std::vector<double>& midY, std::vector<double>& phi, std::vector<double>& Sj)
{
    int N = (int)pointX.size();
    if (N < 2) throw std::invalid_argument("buildPanelsAndInfluences: need at least 2 points");

    int numPanels = N;
    midX.assign(numPanels, 0.0);
    midY.assign(numPanels, 0.0);
    phi.assign(numPanels, 0.0);
    Sj.assign(numPanels, 0.0);

    for (int i = 0; i < numPanels; ++i) {
        int j = (i + 1) % N;
        double dx = pointX[j] - pointX[i];
        double dy = pointY[j] - pointY[i];
        phi[i] = std::atan2(dy, dx); // panel angle (tangent)
        midX[i] = pointX[i] + 0.5 * dx;
        midY[i] = pointY[i] + 0.5 * dy;
        Sj[i] = std::sqrt(dx*dx + dy*dy);
    }

    // allocate matrices
    I.assign(numPanels, Vec(numPanels, 0.0));
    J.assign(numPanels, Vec(numPanels, 0.0));
    K.assign(numPanels, Vec(numPanels, 0.0));
    L.assign(numPanels, Vec(numPanels, 0.0));

    // compute influence integrals
    for (int i = 0; i < numPanels; ++i) {
        for (int j = 0; j < numPanels; ++j) {
            double A = -(midX[i] - pointX[j]) * std::cos(phi[j]) - (midY[i] - pointY[j]) * std::sin(phi[j]);
            double B = std::pow(midX[i] - pointX[j], 2.0) + std::pow(midY[i] - pointY[j], 2.0);
            double Cn = std::sin(phi[i] - phi[j]);
            double Dn = (midY[i] - pointY[j]) * std::cos(phi[i]) - (midX[i] - pointX[j]) * std::sin(phi[i]);
            double Ct = -std::cos(phi[i] - phi[j]);
            double Dt = (midX[i] - pointX[j]) * std::cos(phi[i]) + (midY[i] - pointY[j]) * std::sin(phi[i]);
            double E2 = B - A*A;
            double E = (E2 <= 0.0) ? 0.0 : std::sqrt(E2);

            double tmp = std::pow(Sj[j], 2.0) + 2.0*A*Sj[j] + B;
            if (B <= 0.0) {
                I[i][j] = 0.0;
                J[i][j] = 0.0;
            } else {
                if (E == 0.0) {
                    I[i][j] = 0.0;
                    J[i][j] = 0.0;
                } else {
                    I[i][j] = 0.5 * Cn * std::log(tmp / B) + (Dn - A*Cn) * ( std::atan2(Sj[j] + A, E) - std::atan2(A, E) ) / E;
                    J[i][j] = 0.5 * Ct * std::log(tmp / B) + (Dt - A*Ct) * ( std::atan2(Sj[j] + A, E) - std::atan2(A, E) ) / E;
                }
            }

            double Cn2 = -std::cos(phi[i] - phi[j]);
            double Dn2 = (midX[i] - pointX[j]) * std::cos(phi[i]) + (midY[i] - pointY[j]) * std::sin(phi[i]);
            double Ct2 = std::sin(phi[j] - phi[i]);
            double Dt2 = -(midY[i] - pointY[j]) * std::cos(phi[i]) + (midX[i] - pointX[j]) * std::sin(phi[i]);
            double Ekl = (B - A*A <= 0.0) ? 0.0 : std::sqrt(std::max(B - A*A, 0.0));
            if (Ekl == 0.0) {
                K[i][j] = 0.0;
                L[i][j] = 0.0;
            } else {
                K[i][j] = 0.5 * Cn2 * std::log(tmp / B) + (Dn2 - A*Cn2) * ( std::atan2(Sj[j] + A, Ekl) - std::atan2(A, Ekl) ) / Ekl;
                L[i][j] = 0.5 * Ct2 * std::log(tmp / B) + (Dt2 - A*Ct2) * ( std::atan2(Sj[j] + A, Ekl) - std::atan2(A, Ekl) ) / Ekl;
            }

            if (i == j) {
                I[i][j] = 0.0;
                J[i][j] = 0.0;
                K[i][j] = 0.0;
                L[i][j] = 0.0;
            }

            if (!std::isfinite(I[i][j])) I[i][j] = 0.0;
            if (!std::isfinite(J[i][j])) J[i][j] = 0.0;
            if (!std::isfinite(K[i][j])) K[i][j] = 0.0;
            if (!std::isfinite(L[i][j])) L[i][j] = 0.0;
        }
    }
}

// ----------------------- Panel solver (source + bound vortex) -----------------------
double solvePanelMethod(const std::vector<double>& pointX, const std::vector<double>& pointY,
                        double AoA_deg, double Vinf,
                        std::vector<double>& Cp_out,
                        int debug_write_files = 1)
{
    int Npoints = (int)pointX.size();
    if (Npoints < 3) throw std::invalid_argument("solvePanelMethod: need >=3 points");
    int numPanels = Npoints;

    Mat I, J, K, L;
    std::vector<double> midX, midY, phi, Sj;
    buildPanelsAndInfluences(pointX, pointY, I, J, K, L, midX, midY, phi, Sj);

    int M = numPanels + 1;
    Mat A(M, Vec(M, 0.0));
    Vec B(M, 0.0);

    double AoA = AoA_deg * DegreeToRad;

    for (int i = 0; i < numPanels; ++i) {
        B[i] = -2.0 * Vinf * PI * std::cos(phi[i] + PI/2.0 - AoA);
        for (int j = 0; j < numPanels; ++j) {
            A[i][j] = I[i][j];
            A[i][numPanels] += -K[i][j];
        }
        A[i][i] = PI;
    }

    B[numPanels] = -Vinf * 2.0 * PI * ( std::sin(phi[0] + PI/2.0 - AoA) + std::sin(phi[numPanels - 1] + PI/2.0 - AoA) );
    double sumL = 0.0;
    for (int j = 0; j < numPanels; ++j) {
        sumL += L[0][j] + L[numPanels - 1][j];
        A[numPanels][j] = J[0][j] + J[numPanels - 1][j];
    }
    A[numPanels][numPanels] = 2.0 * PI - sumL;

    Mat invA = inverseLU_v2(A);

    Vec sol(M, 0.0);
    for (int i = 0; i < M; ++i) {
        double s = 0.0;
        for (int j = 0; j < M; ++j) s += invA[i][j] * B[j];
        sol[i] = s;
    }

    Vec lambda(numPanels, 0.0);
    for (int i = 0; i < numPanels; ++i) lambda[i] = sol[i];
    double aeroGamma = sol[numPanels];

    Vec Vt(numPanels, 0.0);
    Cp_out.assign(numPanels, 0.0);
    for (int i = 0; i < numPanels; ++i) {
        double sumV = 0.0;
        for (int j = 0; j < numPanels; ++j) {
            sumV += lambda[j] * J[i][j] / (2.0 * PI) - aeroGamma * L[i][j] / (2.0 * PI);
        }
        double Vtan = Vinf * std::sin(phi[i] + PI / 2.0 - AoA) + 0.5 * aeroGamma + sumV;
        Vt[i] = Vtan;
        Cp_out[i] = 1.0 - (Vtan*Vtan) / (Vinf*Vinf);
    }

    double sumLength = 0.0;
    for (int i = 0; i < numPanels; ++i) sumLength += Sj[i];
    double Cl = sumLength * 2.0 * aeroGamma;

    if (debug_write_files) {
        std::ofstream shape_out("output_shape.dat");
        std::ofstream cp_out("output_cp.dat");
        if (shape_out) {
            shape_out << std::setprecision(12);
            for (int i = 0; i < Npoints; ++i) shape_out << pointX[i] << " " << pointY[i] << "\n";
        }
        if (cp_out) {
            cp_out << std::setprecision(12);
            for (int i = 0; i < numPanels; ++i) cp_out << midX[i] << " " << Cp_out[i] << "\n";
        }
    }

    std::cout << "Computed aeroGamma = " << std::setprecision(8) << aeroGamma << "\n";
    return Cl;
}

// ----------------------- Lifting-line: compute chord for constant Gamma -----------------------
/*
  computeChordForConstantGamma

  r0, r1 : hub and tip radii (m)
  Ncoll  : number of collocation points (and number of odd Fourier terms kept)
  Vinf, rho : flow properties (m/s, kg/m^3)
  Lprime_const_Nperm : desired sectional lift per unit span (N/m) -> Gamma0 = L'/(rho*Vinf)
  alpha_root_deg, alpha_tip_deg : prescribed geometric angle distribution (deg) at root & tip
  c_min, c_max : chord bounds (m)
  outputs: r_out (radial positions), chord_out (m), alpha_out (deg)
*/
void computeChordForConstantGamma(double r0, double r1, int Ncoll,
                                  double Vinf, double rho,
                                  double Lprime_const_Nperm,
                                  double alpha_root_deg, double alpha_tip_deg,
                                  double c_min, double c_max,
                                  std::vector<double>& r_out,
                                  std::vector<double>& chord_out,
                                  std::vector<double>& alpha_out)
{
    if (Ncoll < 3) throw std::invalid_argument("computeChordForConstantGamma: Ncoll must be >=3");
    if (r1 <= r0) throw std::invalid_argument("computeChordForConstantGamma: r1 must be > r0");

    // semi-span used in lifting-line formulas: we treat the blade as a "half-wing"
    double s = r1 - r0;         // treat semi-span = blade length
    double b = 2.0 * s;         // full-span equivalent (for formulas)

    // desired constant circulation
    double Gamma0 = Lprime_const_Nperm / (rho * Vinf);

    // collocation theta points (avoid 0 and pi)
    int N = Ncoll;
    std::vector<double> theta(N);
    for (int m = 0; m < N; ++m) theta[m] = (m + 1) * PI / (N + 1);

    // map theta -> radial coordinate r in [r0, r1]
    // mapping chosen: r(theta) = r0 + s * (1 - cos(theta)) / 2
    r_out.assign(N, 0.0);
    for (int m = 0; m < N; ++m) {
        double th = theta[m];
        r_out[m] = r0 + s * (1.0 - std::cos(th)) * 0.5;
    }

    // prescribed alpha distribution (linear between root & tip by default)
    alpha_out.assign(N, 0.0);
    for (int m = 0; m < N; ++m) {
        double frac = (r_out[m] - r0) / (r1 - r0); // 0..1
        double a_deg = alpha_root_deg + (alpha_tip_deg - alpha_root_deg) * frac;
        alpha_out[m] = a_deg;
    }

    // form odd-n list and analytic A_n coefficients for constant Gamma:
    std::vector<int> nlist(N);
    for (int k = 0; k < N; ++k) nlist[k] = 2 * k + 1; // 1,3,5,...

    std::vector<double> A(N);
    for (int k = 0; k < N; ++k) {
        int n = nlist[k];
        A[k] = 2.0 * Gamma0 / (PI * double(n) * b * Vinf);
    }

    // compute chord at collocation points
    const double sin_eps = 1e-8;
    const double denom_eps = 1e-9;

    chord_out.assign(N, 0.0);

    for (int m = 0; m < N; ++m) {
        double th = theta[m];
        double sinth = std::sin(th);
        double sinth_safe = (std::abs(sinth) < sin_eps) ? ( (sinth < 0) ? -sin_eps : sin_eps ) : sinth;

        // S1 = sum A_n * sin(n theta)
        double S1 = 0.0;
        double S2 = 0.0;
        for (int k = 0; k < N; ++k) {
            double sn = std::sin(nlist[k] * th);
            S1 += A[k] * sn;
            S2 += A[k] * sn * ( double(nlist[k]) / sinth_safe );
        }

        double alpha_rad = alpha_out[m] * DegreeToRad;
        double alpha_L0 = 0.0; // assume zero-lift angle = 0; could be extended to array

        // denominator of mu (mu = c / b)
        double denom = 2.0 * (alpha_rad - alpha_L0) + S2;
        if (std::abs(denom) < denom_eps) {
            // avoid blow-up; set chord to a large value but clamp afterwards
            denom = (denom >= 0.0) ? denom_eps : -denom_eps;
        }

        double mu = (4.0 * S1) / denom; // mu = c/b
        double c_local = b * mu;

        // clamp to manufacturable bounds
        if (std::isnan(c_local) || !std::isfinite(c_local)) c_local = c_min;
        c_local = std::max(c_min, std::min(c_max, c_local));

        chord_out[m] = c_local;
    }

    // OPTIONAL: simple smoothing (3-point moving average) to remove tiny oscillations
    std::vector<double> smooth = chord_out;
    if (N >= 5) {
        for (int iter = 0; iter < 2; ++iter) {
            for (int m = 1; m < N-1; ++m) {
                smooth[m] = 0.25 * chord_out[m-1] + 0.5 * chord_out[m] + 0.25 * chord_out[m+1];
            }
            chord_out = smooth;
        }
    }

    // write to file
    std::ofstream ofs("constant_gamma_chord.dat");
    ofs << std::setprecision(10);
    ofs << "# r(m)   chord(m)   alpha_deg\n";
    for (int m = 0; m < N; ++m) {
        ofs << r_out[m] << "   " << chord_out[m] << "   " << alpha_out[m] << "\n";
    }
    ofs.close();

    std::cout << "Wrote constant_gamma_chord.dat (" << N << " stations) with chord clamped to [" 
              << c_min << ", " << c_max << "] m\n";
}

// ----------------------- Main demonstrating the computation -----------------------
int main()
{
    try {
        // ---------------- user-changeable parameters ----------------
        int resolution = 500;      // number of sample points around circle -> approximate panels
        double disX = 0.04;        // x-offset of circle center (negative/positive shapes camber)
        double disY = 0.10;        // y-offset of circle center
        double backFat = 1.04;     // 'backFat' scaling parameter used in generation
        double AoA_deg = 2.0;      // angle of attack in degrees (used as alpha root default)
        double Vinf = 1.0;         // free-stream speed / axial inflow (m/s)

        // blade radial / lifting-line parameters
        double r0 = 0.08;          // hub radius (m)
        double r1 = 0.14;          // tip radius (m)
        int nRadial = 50;          // number of radial collocation stations for lifting-line
        double rpm = 2400.0;       // rotational speed (kept for potential later use)

        // desired sectional lift (N/m) -> constant across span
        // You can replace this with L_total/(2*s) if you prefer to specify total lift
        double Lprime_const_N_per_m = 500.0;   // N/m (example). Change as needed.

        // prescribed geometric alpha distribution (deg) at root and tip (you can change)
        double alpha_root_deg = AoA_deg;  // use global AoA as root angle
        double alpha_tip_deg = 2.5;       // deg at tip

        // chord bounds (manufacturable)
        double c_min = 0.08;   // minimum chord (m)
        double c_max = 0.16;    // maximum chord (m)
        // ----------------------------------------------------------

        // generate shape (Joukowsky 2D airfoil coords) - remains available for section evaluation
        std::vector<double> px, py;
        genJoukowsky(disX, disY, backFat, resolution, px, py);

        // Solve panel method and compute Cl for the base section
        std::vector<double> Cp;
        double Cl = solvePanelMethod(px, py, AoA_deg, Vinf, Cp, /*write files*/ 1);

        std::cout << "Estimated sectional lift coefficient Cl (2D panel) = " << std::setprecision(8) << Cl << "\n";
        std::cout << "Wrote output_shape.dat (coords) and output_cp.dat (midX, Cp per panel midpoint)\n";

        // ---------------- compute chord & twist from lifting-line (constant Gamma) ----------------
        std::vector<double> r_stations, chord_stations, alpha_stations_deg;
        computeChordForConstantGamma(r0, r1, nRadial, Vinf, /*rho*/1.225,
                                     Lprime_const_N_per_m,
                                     alpha_root_deg, alpha_tip_deg,
                                     c_min, c_max,
                                     r_stations, chord_stations, alpha_stations_deg);

        // Write twist (alpha) vs r (for clarity)
        std::ofstream tfs("constant_gamma_twist.dat");
        tfs << std::setprecision(10);
        tfs << "# r(m)   alpha_deg\n";
        for (size_t i = 0; i < r_stations.size(); ++i) tfs << r_stations[i] << "   " << alpha_stations_deg[i] << "\n";
        tfs.close();
        std::cout << "Wrote constant_gamma_twist.dat\n";

        std::cout << "Done.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}

