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
#include <array>

using Vec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;

constexpr double PI = 3.14159265358979323846;
constexpr double RadToDegree = 180.0 / PI;
constexpr double DegreeToRad = PI / 180.0;

//Geometry

using Vec3 = std::array<double,3>;

std::vector<std::vector<Vec3>> buildBladeGeometry(
        const std::vector<double>& r,
        const std::vector<double>& chord,
        const std::vector<double>& twist,
        const std::vector<std::array<double,2>>& baseAirfoil  // (x,y) on unit chord
        ){
    int NR = r.size();
    int NP = baseAirfoil.size();

    std::vector<std::vector<Vec3>> rings(NR, std::vector<Vec3>(NP));

    for(int i=0; i<NR; ++i){
        double c = chord[i];
        double t = -twist[i];
        double ct = cos(t), st = sin(t);

        for(int j=0; j<NP; ++j){
            double x = baseAirfoil[j][0] * c;
            double y = baseAirfoil[j][1] * c;

            // rotate around Z for twist
            double xr = x * ct - y * st;
            double yr = x * st + y * ct;

            // translate to radius r[i], blade lies in x-y plane (r-axis along x)
            rings[i][j] = { r[i], xr, yr };
        }
    }

    return rings;
}


// Write ASCII STL
void writeSTL(
        const std::string& filename,
        const std::vector<std::vector<Vec3>>& bladeRings   // bladeRings[i][j] = point j of section i
        ){
    std::ofstream out(filename);
    out << "solid blade\n";

    int NR = bladeRings.size();
    int NP = bladeRings[0].size();  // points per section

    for(int i=0; i<NR-1; ++i){
        for(int j=0; j<NP-1; ++j){
            // quad split into two triangles
            Vec3 p1 = bladeRings[i][j];
            Vec3 p2 = bladeRings[i+1][j];
            Vec3 p3 = bladeRings[i+1][j+1];
            Vec3 p4 = bladeRings[i][j+1];

            auto writeTri = [&](Vec3 a, Vec3 b, Vec3 c){
                // Compute normal (not normalized)
                Vec3 u = {b[0]-a[0],b[1]-a[1],b[2]-a[2]};
                Vec3 v = {c[0]-a[0],c[1]-a[1],c[2]-a[2]};
                Vec3 n = { u[1]*v[2] - u[2]*v[1],
                    u[2]*v[0] - u[0]*v[2],
                    u[0]*v[1] - u[1]*v[0] };
                out << "  facet normal " << n[0] << " " << n[1] << " " << n[2] << "\n";
                out << "    outer loop\n";
                out << "      vertex " << a[0] << " " << a[1] << " " << a[2] << "\n";
                out << "      vertex " << b[0] << " " << b[1] << " " << b[2] << "\n";
                out << "      vertex " << c[0] << " " << c[1] << " " << c[2] << "\n";
                out << "    endloop\n";
                out << "  endfacet\n";
            };

            writeTri(p1,p2,p3);
            writeTri(p1,p3,p4);
        }
    }

    out << "endsolid blade\n";
    out.close();
}


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
        //double Vtan = Vinf * std::sin(phi[i] + PI / 2.0 - AoA) + 0.5 * aeroGamma + sumV;
        double Vtan = Vinf * std::sin(phi[i] + PI / 2.0 - AoA) + 0.5 * aeroGamma + sumV;
        Vt[i] = Vtan;
        Cp_out[i] = 1.0 - (Vtan*Vtan) / (Vinf*Vinf);
    }

    double sumLength = 0.0;
    for (int i = 0; i < numPanels; ++i) sumLength += Sj[i];
    //double Cl = sumLength * 2.0 * aeroGamma;
    double Cl = 2.0 * aeroGamma / Vinf;

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
void computeChordForConstantGamma(
        double r0, double r1, int Ncoll,
        double Vinf, double rho,
        double Lprime_root, double Lprime_tip, // linear lift endpoints (used to compute total lift)
        double alpha_root_deg, double alpha_tip_deg,
        double c_min, double c_max,
        double alpha_L0_rad,
        std::vector<double>& r_out,
        std::vector<double>& chord_out,
        std::vector<double>& alpha_out,
        double rpm, std::vector<double> px, std::vector<double> py)
{
    if (Ncoll < 3) throw std::invalid_argument("Ncoll must be >=3");
    if (r1 <= r0) throw std::invalid_argument("r1 must be > r0");

    int N = Ncoll;
    double span = r1 - r0;

    // -------------------------
    // Cosine-spaced radial stations
    // -------------------------
    r_out.resize(N);
    for (int m = 0; m < N; ++m) {
        double theta = (m + 1) * M_PI / (N + 1);
        r_out[m] = r0 + 0.5 * span * (1.0 - std::cos(theta));
    }

    // -------------------------
    // Linear geometric twist
    // -------------------------
    alpha_out.resize(N);
    for (int m = 0; m < N; ++m) {
        double frac = (r_out[m] - r0) / span;
        alpha_out[m] = alpha_root_deg + (alpha_tip_deg - alpha_root_deg) * frac;
    }
    // -------------------------
    // Chord calculation using PANEL-METHOD Cl at each station
    // -------------------------
    chord_out.resize(N);
    std::vector<double> Vlocal(N);
    const double eps = 1e-6;
    double omega = rpm * 2.0 * M_PI / 60.0; // RPM -> rad/s

    // -------------------------
    // Build polynomial L'(r) shape (Option 3)
    // We create a shape function P(s) = s^n with s in [0,1] (s=0 at r=r0).
    // Choose exponent n (2 => quadratic root taper). Then scale so that the
    // total lift per blade equals the total lift that the original linear
    // profile (Lprime_root..Lprime_tip) would have produced.
    // -------------------------
    int poly_n = 1.6; // quadratic root taper (change to 3 or 4 for stronger root decay)

    // compute the target total lift per blade implied by the linear endpoints
    // (integral of a linear function between r0 and r1)
    double Ltotal_from_linear = 0.5 * (Lprime_root + Lprime_tip) * span;

    // normalization: integral_0^1 s^n ds = 1/(n+1)
    double norm = 1.0 / double(poly_n + 1);

    // constant C such that: integral_{r0}^{r1} C * s^n dr = Ltotal_from_linear
    // => span * C * norm = Ltotal_from_linear  => C = Ltotal_from_linear / (span * norm)
    double C = (span * norm > 0.0) ? (Ltotal_from_linear / (span * norm)) : 0.0;

    for (int m = 0; m < N; ++m)
    {
        double r = r_out[m];
        double frac = (r - r0) / span;
        if (frac < 0.0) frac = 0.0;
        if (frac > 1.0) frac = 1.0;

        // ---- POLYNOMIAL SECTIONAL LIFT TARGET ----
        double s = frac;
        double Lprime_r = C * std::pow(s, double(poly_n)); // L'(r) = C * s^n

        // Local relative velocity (hover inflow + rotation)
        double Vloc = std::sqrt(Vinf * Vinf + std::pow(omega * r, 2));
        Vlocal[m] = Vloc;

        // Local geometric AoA
        double alpha_geo_deg = alpha_out[m];
        double alpha_eff_rad = alpha_geo_deg * DegreeToRad - alpha_L0_rad;
        double alpha_eff_deg = alpha_eff_rad * RadToDegree;

        // ---- PANEL METHOD CL ----
        std::vector<double> Cp_dummy;
        double Cl_local = solvePanelMethod(px, py, alpha_eff_deg, Vloc, Cp_dummy, 0);

        // handle singularities and extremely small Cl
        if (!std::isfinite(Cl_local) || std::fabs(Cl_local) < eps)
            Cl_local = eps;

        // ---- BLADE ELEMENT CHORD ----
        double dyn = 0.5 * rho * Vloc * Vloc;
        double c_local = (dyn > 0.0) ? (Lprime_r / (dyn * Cl_local)) : c_max;

        // clamp for manufacturability
        if (c_local < c_min && r < r0 + (r1-r0) * 0.2) c_local = c_min;
        if (c_local > c_max && r < r0 + (r1-r0) * 0.2) c_local = c_max;

        chord_out[m] = c_local;

        std::cout << "r=" << r
            << " s=" << s
            << " Vloc=" << Vloc
            << " alpha_eff_deg=" << alpha_eff_deg
            << " Cl=" << Cl_local
            << " L'=" << Lprime_r
            << " chord=" << c_local
            << std::endl;
    }

    // -------------------------
    // Write output
    // -------------------------
    std::ofstream ofs("linear_lift_chord.dat");
    ofs << std::setprecision(10);
    ofs << "# r(m)   chord(m)   alpha_deg\n";
    for (int m = 0; m < N; ++m)
        ofs << r_out[m] << "   " << chord_out[m] << "   " << alpha_out[m] << "\n";
    ofs.close();

    std::cout << "Wrote linear_lift_chord.dat (" << N
        << " stations) using local inflow V(r)\n";

}

double computeZeroLiftAoAPanel(
        const std::vector<double>& px,
        const std::vector<double>& py,
        double Vinf)
{
    // pick 3 sample AoA values near expected linear range
    double a1 = -15.0;
    double a2 = 0.0;
    double a3 = +15.0;

    std::vector<double> Cp_dummy;

    double Cl1 = solvePanelMethod(px, py, a1, Vinf, Cp_dummy, 0);
    double Cl2 = solvePanelMethod(px, py, a2, Vinf, Cp_dummy, 0);
    double Cl3 = solvePanelMethod(px, py, a3, Vinf, Cp_dummy, 0);

    // linear fit Cl = A * alpha + B  (alpha in radians)
    double x1 = a1 * DegreeToRad;
    double x2 = a2 * DegreeToRad;
    double x3 = a3 * DegreeToRad;

    // least-squares fit
    double Sx  = x1 + x2 + x3;
    double Sy  = Cl1 + Cl2 + Cl3;
    double Sxx = x1*x1 + x2*x2 + x3*x3;
    double Sxy = x1*Cl1 + x2*Cl2 + x3*Cl3;

    double denom = 3.0*Sxx - Sx*Sx;
    double A = (3.0*Sxy - Sx*Sy) / denom;
    double B = (Sy - A*Sx) / 3.0;

    // zero-lift angle: Cl = 0 → A*alpha_0 + B = 0
    double alpha0_rad = -B / A;

    std::cout << "Panel-derived zero-lift angle α_L0 = "
        << alpha0_rad * RadToDegree << " deg\n";

    return alpha0_rad;
}

// ----------------------- Main demonstrating the computation -----------------------
int main()
{
    try {
        // ---------------- user-changeable parameters ----------------
        int resolution = 500;      // number of sample points around circle -> approximate panels
        double disX = 0.02;        // x-offset of circle center (negative/positive shapes camber)
        double disY = 0.14;        // y-offset of circle center
        double backFat = 1.06;     // 'backFat' scaling parameter used in generation
        double AoA_deg = 0.0;    // angle of attack in degrees (used as alpha root default)

        // blade radial / lifting-line parameters
        double r0 = 0.03;          // hub radius (m)
        double r1 = 0.16;           // tip radius (m)
        int nRadial = 40;           // number of radial collocation stations for lifting-line
        double rpm = 4500.0;        // rotational speed (kept for potential later use)

        double rho = 1.225;         // air density
        double T = 10.0;            // rotor thrust in N (example)
        double R = r1;              // rotor tip radius from your code
        int B = 2;                   // number of blades (adjust as needed)

        double rotorArea = 3.14159265358979323846 * R * R;
        double Vinf = std::sqrt(T / (2.0 * rho * rotorArea)); // hover inflow

        // desired sectional lift (N/m) -> constant across span
        double Lprime_const_N_per_m = T / (B * (r1 - r0)); // distribute total lift across blades
        double k = 0.9; // tip/root ratio
        double Lprime_avg = Lprime_const_N_per_m; // from previous constant lift
        double Lprime_root = 2.0 * Lprime_avg / (1.0 + k);
        double Lprime_tip = k * Lprime_root;


        // prescribed geometric alpha distribution (deg) at root and tip
        double alpha_root_deg = 16.0;
        double alpha_tip_deg = 10.0; 

        // chord bounds (manufacturable)
        double c_min = 0.02;
        double c_max = 0.126;
        // ----------------------------------------------------------

        // generate shape (Joukowsky 2D airfoil coords)
        std::vector<double> px, py;
        genJoukowsky(disX, disY, backFat, resolution, px, py);

        // Solve panel method and compute Cl for the base section
        std::vector<double> Cp;
        double Cl = solvePanelMethod(px, py, AoA_deg, Vinf, Cp, 1);

        std::cout << "Estimated sectional lift coefficient Cl (2D panel) = " << std::setprecision(8) << Cl << "\n";
        std::cout << "Wrote output_shape.dat and output_cp.dat\n";

        // compute zero-lift AoA from panel method
        double alpha_L0_rad = computeZeroLiftAoAPanel(px, py, Vinf);

        // lifting-line chord solver with REAL airfoil α_L0
        std::vector<double> r_stations, chord_stations, alpha_stations_deg;
        computeChordForConstantGamma(
                r0, r1, nRadial,
                Vinf, rho,
                Lprime_root, Lprime_tip,
                alpha_root_deg, alpha_tip_deg,
                c_min, c_max,
                alpha_L0_rad,                     // <<< add this
                r_stations, chord_stations, alpha_stations_deg,
                rpm, px, py
                );

        // Write twist (alpha) vs r
        std::ofstream tfs("constant_gamma_twist.dat");
        tfs << std::setprecision(10);
        tfs << "# r(m)   alpha_deg\n";
        for (size_t i = 0; i < r_stations.size(); ++i) tfs << r_stations[i] << "   " << alpha_stations_deg[i] << "\n";
        tfs.close();
        std::cout << "Wrote constant_gamma_twist.dat\n";

        // ---------------- compute total lift ----------------
        double L_total_per_blade = Lprime_const_N_per_m * (r1 - r0);
        double L_total_rotor = L_total_per_blade * B;

        std::cout << "Estimated total lift per blade = " << L_total_per_blade << " N\n";
        std::cout << "Estimated total rotor lift = " << L_total_rotor << " N\n";

        // ----------------------------------------------------------
        // Build 3D blade geometry → STL
        std::vector<std::array<double,2>> airfoilPts;
        airfoilPts.reserve(px.size());
        for (size_t i = 0; i < px.size(); ++i)
            airfoilPts.push_back({px[i], py[i]});

        std::vector<double> twist_rad(alpha_stations_deg.size());
        for (size_t i = 0; i < alpha_stations_deg.size(); ++i)
            twist_rad[i] = alpha_stations_deg[i] * DegreeToRad;

        auto rings = buildBladeGeometry(
                r_stations,
                chord_stations,
                twist_rad,
                airfoilPts
                );

        writeSTL("blade.stl", rings);
        std::cout << "Wrote blade.stl (3D blade geometry)\n";

        std::cout << "Done.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

}


