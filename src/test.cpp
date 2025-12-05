// joukowsky_panel_with_blade_twist.cpp
// Build with: g++ -std=c++17 -O2 joukowsky_panel_with_blade_twist.cpp -o joukowsky_panel
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

    // Ensure closed loop: last point equals first (for safety; panels use consecutive indices)
    if (px.size() >= 2 && (px.front() != px.back() || py.front() != py.back())) {
        // We'll keep as open closed loop by replicating first as last if needed
        // but panels expect points [0..N-1] with last panel between N-1 and 0.
    }
}

// ----------------------- Panel geometry and influence integrals -----------------------
// Given boundary points pointX, pointY (size N), build panels between i and (i+1)%N (so numPanels = N)
// We create arrays: midX/midY (control points), phi (panel angle), Sj (panel length)
// and compute influence matrices I,J,K,L (numPanels x numPanels) following formulas used earlier.
void buildPanelsAndInfluences(const std::vector<double>& pointX, const std::vector<double>& pointY,
                              Mat& I, Mat& J, Mat& K, Mat& L,
                              std::vector<double>& midX, std::vector<double>& midY, std::vector<double>& phi, std::vector<double>& Sj)
{
    int N = (int)pointX.size();
    if (N < 2) throw std::invalid_argument("buildPanelsAndInfluences: need at least 2 points");

    int numPanels = N; // using N points with panels between i and (i+1)%N; this matches your original approach where panel count = resolution
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
            // point j panel from pointX[j] to pointX[j+1]
            double A = -(midX[i] - pointX[j]) * std::cos(phi[j]) - (midY[i] - pointY[j]) * std::sin(phi[j]);
            double B = std::pow(midX[i] - pointX[j], 2.0) + std::pow(midY[i] - pointY[j], 2.0);
            double Cn = std::sin(phi[i] - phi[j]);
            double Dn = (midY[i] - pointY[j]) * std::cos(phi[i]) - (midX[i] - pointX[j]) * std::sin(phi[i]);
            double Ct = -std::cos(phi[i] - phi[j]);
            double Dt = (midX[i] - pointX[j]) * std::cos(phi[i]) + (midY[i] - pointY[j]) * std::sin(phi[i]);
            double E2 = B - A*A;
            double E = (E2 <= 0.0) ? 0.0 : std::sqrt(E2);

            double tmp = std::pow(Sj[j], 2.0) + 2.0*A*Sj[j] + B;
            // protect logs and divide by zero in atan differences
            if (B <= 0.0) {
                I[i][j] = 0.0;
                J[i][j] = 0.0;
            } else {
                if (E == 0.0) {
                    // limiting form when E -> 0 (panel very close / singular), set to 0 for stability
                    I[i][j] = 0.0;
                    J[i][j] = 0.0;
                } else {
                    I[i][j] = 0.5 * Cn * std::log(tmp / B) + (Dn - A*Cn) * ( std::atan2(Sj[j] + A, E) - std::atan2(A, E) ) / E;
                    J[i][j] = 0.5 * Ct * std::log(tmp / B) + (Dt - A*Ct) * ( std::atan2(Sj[j] + A, E) - std::atan2(A, E) ) / E;
                }
            }

            // K,L (different sign choices for tangential influence)
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
                // self-influence singular terms — set to zero as in original snippet
                I[i][j] = 0.0;
                J[i][j] = 0.0;
                K[i][j] = 0.0;
                L[i][j] = 0.0;
            }

            // sanitize NaN/Inf
            if (!std::isfinite(I[i][j])) I[i][j] = 0.0;
            if (!std::isfinite(J[i][j])) J[i][j] = 0.0;
            if (!std::isfinite(K[i][j])) K[i][j] = 0.0;
            if (!std::isfinite(L[i][j])) L[i][j] = 0.0;
        }
    }
}

// ----------------------- Panel solver (source + bound vortex) -----------------------
// Uses system size (numPanels+1) to solve for lambda[0..numPanels-1] (source strengths) and aeroGamma (= last unknown).
// AoA in degrees, Vinf velocity
// Returns Cp vector per panel and Cl estimate.
double solvePanelMethod(const std::vector<double>& pointX, const std::vector<double>& pointY,
                        double AoA_deg, double Vinf,
                        std::vector<double>& Cp_out,
                        int debug_write_files = 1)
{
    int Npoints = (int)pointX.size();
    if (Npoints < 3) throw std::invalid_argument("solvePanelMethod: need >=3 points");
    int numPanels = Npoints;

    // build geometry and influence matrices
    Mat I, J, K, L;
    std::vector<double> midX, midY, phi, Sj;
    buildPanelsAndInfluences(pointX, pointY, I, J, K, L, midX, midY, phi, Sj);

    int M = numPanels + 1; // system size
    Mat A(M, Vec(M, 0.0));
    Vec B(M, 0.0);

    double AoA = AoA_deg * DegreeToRad;

    // build equations
    for (int i = 0; i < numPanels; ++i) {
        // RHS: -2 Vinf * PI * cos( phi[i] + PI/2 - AoA )
        B[i] = -2.0 * Vinf * PI * std::cos(phi[i] + PI/2.0 - AoA);
        for (int j = 0; j < numPanels; ++j) {
            A[i][j] = I[i][j];
            A[i][numPanels] += -K[i][j];
        }
        A[i][i] = PI; // diagonal
    }

    // Kutta condition (last equation)
    B[numPanels] = -Vinf * 2.0 * PI * ( std::sin(phi[0] + PI/2.0 - AoA) + std::sin(phi[numPanels - 1] + PI/2.0 - AoA) );
    double sumL = 0.0;
    for (int j = 0; j < numPanels; ++j) {
        sumL += L[0][j] + L[numPanels - 1][j];
        A[numPanels][j] = J[0][j] + J[numPanels - 1][j];
    }
    A[numPanels][numPanels] = 2.0 * PI - sumL;

    // invert A
    Mat invA = inverseLU_v2(A);

    // multiply invA * B -> solution vector (size M)
    Vec sol(M, 0.0);
    for (int i = 0; i < M; ++i) {
        double s = 0.0;
        for (int j = 0; j < M; ++j) s += invA[i][j] * B[j];
        sol[i] = s;
    }

    // first numPanels entries are lambda (source strengths), last entry is aeroGamma
    Vec lambda(numPanels, 0.0);
    for (int i = 0; i < numPanels; ++i) lambda[i] = sol[i];
    double aeroGamma = sol[numPanels];

    // compute tangential velocity Vt and Cp on each panel
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

    // compute Cl estimate: follow original formula aeroCl = sumLength * 2 * aeroGamma
    double sumLength = 0.0;
    for (int i = 0; i < numPanels; ++i) sumLength += Sj[i];
    double Cl = sumLength * 2.0 * aeroGamma;

    // write output files if requested
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

// ----------------------- Blade twist calculation (new) -----------------------
/*
  Compute local geometric pitch angle theta(r) for blade elements between r0 (hub) and r1 (tip)
  using simple blade-element inflow geometry:

    phi(r)   = atan( Vaxial / ( omega * r ) )      [inflow angle, radians]
    theta(r) = phi(r) + alpha_opt                   [geometric pitch, radians]

  Where:
    - Vaxial is axial inflow (m/s) (use Vinf)
    - omega = 2*pi*rpm/60 (rad/s)
    - alpha_opt is target local angle of attack (radians)

  Outputs file blade_twist.dat with columns:
      r    phi_deg    alpha_opt_deg    theta_deg
*/
void computeBladeTwist(double r0, double r1, int nRadial,
                       double Vaxial, double rpm, double alpha_opt_deg,
                       const std::string& outfn = "blade_twist.dat")
{
    if (nRadial <= 0) throw std::invalid_argument("computeBladeTwist: nRadial must be positive");
    if (r1 <= r0) throw std::invalid_argument("computeBladeTwist: r1 must be > r0");

    double alpha_opt = alpha_opt_deg * DegreeToRad;
    double omega = 2.0 * PI * rpm / 60.0;

    std::ofstream ofs(outfn);
    ofs << std::setprecision(8);
    ofs << "# r (m)    phi_deg    alpha_opt_deg    theta_deg\n";

    // avoid exact r==0; start slightly above r0 if r0==0
    double eps = 1e-6 * std::max(1.0, r1);

    for (int i = 0; i <= nRadial; ++i) {
        double t = double(i) / double(nRadial);
        double r = r0 + t * (r1 - r0);
        if (r <= 0.0) r = eps;

        // local inflow angle phi
        double denom = omega * r;
        double phi = 0.0;
        if (denom <= 0.0) {
            // if omega==0, flow is purely axial -> phi = pi/2 (90 deg). set phi to pi/2.
            phi = PI / 2.0;
        } else {
            phi = std::atan( Vaxial / denom );
            // ensure phi is in [0, pi/2)
            if (phi < 0.0) phi = 0.0;
        }

        double theta = phi + alpha_opt;

        ofs << r << "    " << (phi * RadToDegree) << "    " << alpha_opt_deg << "    " << (theta * RadToDegree) << "\n";
    }

    ofs.close();
    std::cout << "Wrote blade twist distribution to " << outfn << "\n";
}

// ----------------------- Main demonstrating the computation -----------------------
int main()
{
    try {
        // Parameters (user can change)
        int resolution = 200;      // number of sample points around circle -> approximate panels
        double disX = 0.04;        // x-offset of circle center (negative/positive shapes camber)
        double disY = 0.08;        // y-offset of circle center
        double backFat = 1.02;     // 'backFat' scaling parameter used in generation
        double AoA_deg = 4.0;      // angle of attack in degrees (used as alpha_opt in blade twist)
        double Vinf = 1.0;         // free-stream speed / axial inflow (m/s)

        // Blade twist parameters (new)
        double r0 = 0.10;          // hub radius (m) -- change to your hub radius
        double r1 = 1.00;          // tip radius (m) -- blade length
        int nRadial = 50;          // number of radial stations
        double rpm = 2400.0;       // rotational speed in RPM
        // Note: AoA_deg used as alpha_opt (target local angle of attack)

        // generate shape
        std::vector<double> px, py;
        genJoukowsky(disX, disY, backFat, resolution, px, py);

        // Solve panel method and compute Cl
        std::vector<double> Cp;
        double Cl = solvePanelMethod(px, py, AoA_deg, Vinf, Cp, /*write files*/ 1);

        std::cout << "Estimated sectional lift coefficient Cl = " << std::setprecision(8) << Cl << "\n";
        std::cout << "Wrote output_shape.dat (coords) and output_cp.dat (midX, Cp per panel midpoint)\n";

        // Compute blade twist distribution using the same axial inflow Vinf and AoA_deg as alpha_opt
        computeBladeTwist(r0, r1, nRadial, Vinf, rpm, AoA_deg, "blade_twist.dat");

        std::cout << "Done.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}

