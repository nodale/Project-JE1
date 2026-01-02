// ============================================================
//  PHYSICALLY CONSISTENT ROTOR / AIRFOIL SOLVER
//  - Joukowski airfoil
//  - Vortex panel method
//  - Prandtl lifting-line theory
//  - Induced velocity
//  - STL output
//
//  Author: (rewritten for correctness)
// ============================================================

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

constexpr double PI = 3.141592653589793;
constexpr double DEG2RAD = PI / 180.0;

// ============================================================
// JOUKOWSKI AIRFOIL
// ============================================================

void generateJoukowski(
        int N,
        double x0,
        double y0,
        std::vector<double>& x,
        std::vector<double>& y)
{
    const double a = 1.0;
    x.resize(N);
    y.resize(N);

    for (int i = 0; i < N; ++i)
    {
        double th = 2.0 * PI * i / (N - 1);
        double zr = a * cos(th) + x0;
        double zi = a * sin(th) + y0;

        double r2 = zr*zr + zi*zi;
        double xr = zr + a*a * zr / r2;
        double yi = zi - a*a * zi / r2;

        x[i] = xr;
        y[i] = yi;
    }
}
double solvePanelCl(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double alpha)
{
    const int N = x.size() - 1;   // closed contour
    const double Vinf = 1.0;

    std::vector<double> xc(N), yc(N), phi(N), S(N);

    for (int i = 0; i < N; ++i)
    {
        double dx = x[i+1] - x[i];
        double dy = y[i+1] - y[i];
        S[i] = std::sqrt(dx*dx + dy*dy);
        phi[i] = std::atan2(dy, dx);
        xc[i] = 0.5 * (x[i+1] + x[i]);
        yc[i] = 0.5 * (y[i+1] + y[i]);
    }

    // Influence matrix
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
    std::vector<double> b(N, 0.0);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j)
            {
                A[i][j] = 0.5;   // SELF INFLUENCE (CRITICAL)
                continue;
            }

            double dx = xc[i] - x[j];
            double dy = yc[i] - y[j];
            double r2 = dx*dx + dy*dy;

            A[i][j] = -(1.0 / (2 * M_PI)) *
                      (dx * sin(phi[j]) - dy * cos(phi[j])) / r2;
        }

        b[i] = -Vinf * sin(phi[i] - alpha);
    }

    // Enforce Kutta condition:
    // Replace last equation
    for (int j = 0; j < N; ++j)
        A[N-1][j] = 0.0;

    A[N-1][0] = 1.0;
    A[N-1][N-1] = 1.0;
    b[N-1] = 0.0;

    // Solve by Gaussian elimination
    for (int i = 0; i < N; ++i)
    {
        double pivot = A[i][i];
        if (std::abs(pivot) < 1e-10)
            throw std::runtime_error("Singular matrix");

        for (int j = i; j < N; ++j)
            A[i][j] /= pivot;
        b[i] /= pivot;

        for (int k = 0; k < N; ++k)
        {
            if (k == i) continue;
            double f = A[k][i];
            for (int j = i; j < N; ++j)
                A[k][j] -= f * A[i][j];
            b[k] -= f * b[i];
        }
    }

    // Circulation
    double Gamma = 0.0;
    for (int i = 0; i < N; ++i)
        Gamma += b[i] * S[i];

    double chord =
        *std::max_element(x.begin(), x.end()) -
        *std::min_element(x.begin(), x.end());

    return 2.0 * Gamma / (Vinf * chord);
}


void liftingLine(
        int N,
        double R,
        double V,
        double alpha0,
        std::vector<double>& r,
        std::vector<double>& chord)
{
    r.resize(N);
    chord.resize(N);

    std::vector<double> A(N,0.0);

    for (int i = 0; i < N; ++i)
    {
        double theta = PI * (i + 1) / (N + 1);
        r[i] = 0.5 * R * (1 - cos(theta));

        A[i] = (alpha0) / (PI * (i+1));
    }

    for (int i = 0; i < N; ++i)
    {
        double Gamma = 2 * R * V * A[i];
        chord[i] = Gamma / (0.5 * V * V);
    }
}

struct BladeSection {
    double r;
    double chord;
    double twist_deg;
};

std::vector<BladeSection> generatePropellerBlade(
    int N,
    double R,
    double rpm,
    double rho,
    double thrust,
    double Cl_design,
    double alpha0_rad)
{
    std::vector<BladeSection> blade(N);

    double omega = rpm * 2.0 * M_PI / 60.0;
    double vi = std::sqrt(thrust / (2.0 * rho * M_PI * R * R));

    // Compute circulation scale
    double Gamma0 = thrust / (rho * 2.0 * M_PI * R);

    for (int i = 0; i < N; ++i)
    {
        double mu = (i + 1.0) / (N + 1.0);
        double r = mu * R;

        // Prandtl optimal circulation
        double Gamma = Gamma0 * std::sqrt(1.0 - mu * mu);

        double V = std::sqrt((omega * r)*(omega * r) + vi*vi);

        // Chord from lift
        double chord = 2.0 * Gamma / (V * Cl_design);

        // Inflow angle
        double phi = std::atan2(vi, omega * r);

        // Twist = inflow + alpha
        double twist = phi + alpha0_rad;

        blade[i] = {
            r,
            chord,
            twist * 180.0 / M_PI
        };
    }

    return blade;
}
struct Vec3 {
    double x,y,z;
};

std::vector<std::vector<Vec3>>
buildBlade(
    const std::vector<double>& r,
    const std::vector<double>& chord,
    const std::vector<double>& twist_deg,
    const std::vector<double>& airfoil_x,
    const std::vector<double>& airfoil_y)
{
    int Nr = r.size();
    int Np = airfoil_x.size();

    std::vector<std::vector<Vec3>> blade(Nr,
        std::vector<Vec3>(Np));

    for (int i = 0; i < Nr; ++i)
    {
        double c = chord[i];
        double t = twist_deg[i] * M_PI / 180.0;

        double ct = cos(t);
        double st = sin(t);

        for (int j = 0; j < Np; ++j)
        {
            double x = airfoil_x[j] * c;
            double y = airfoil_y[j] * c;

            // rotate by twist
            double xr =  x * ct - y * st;
            double zr =  x * st + y * ct;

            blade[i][j] = {
                r[i],
                xr,
                zr
            };
        }
    }
    return blade;
}
void writeSTL(
    const std::string& filename,
    const std::vector<std::vector<Vec3>>& blade)
{
    std::ofstream out(filename);
    out << "solid blade\n";

    int N = blade.size();
    int M = blade[0].size();

    auto writeTri = [&](Vec3 a, Vec3 b, Vec3 c)
    {
        Vec3 u{b.x-a.x, b.y-a.y, b.z-a.z};
        Vec3 v{c.x-a.x, c.y-a.y, c.z-a.z};
        Vec3 n{
            u.y*v.z - u.z*v.y,
            u.z*v.x - u.x*v.z,
            u.x*v.y - u.y*v.x
        };

        out << " facet normal " << n.x << " " << n.y << " " << n.z << "\n";
        out << "  outer loop\n";
        out << "   vertex " << a.x << " " << a.y << " " << a.z << "\n";
        out << "   vertex " << b.x << " " << b.y << " " << b.z << "\n";
        out << "   vertex " << c.x << " " << c.y << " " << c.z << "\n";
        out << "  endloop\n endfacet\n";
    };

    for (int i = 0; i < N-1; ++i)
        for (int j = 0; j < M-1; ++j)
        {
            writeTri(blade[i][j], blade[i+1][j], blade[i+1][j+1]);
            writeTri(blade[i][j], blade[i+1][j+1], blade[i][j+1]);
        }

    out << "endsolid blade\n";
}


int main()
{
    // --- Airfoil
    std::vector<double> ax, ay;
    generateJoukowski(120, 0.08, 0.1, ax, ay);

    // --- Blade physics
    auto blade = generatePropellerBlade(
        40,        // radial stations
        0.14,      // radius
        6500,      // RPM
        1.225,     // rho
        36.0,      // thrust
        0.8,       // Cl
        -4.0 * M_PI/180.0
    );

    std::vector<double> r, c, t;
    for (auto& b : blade) {
        r.push_back(b.r);
        c.push_back(b.chord);
        t.push_back(b.twist_deg);
    }

    auto geom = buildBlade(r, c, t, ax, ay);

    writeSTL("propeller.stl", geom);

    std::cout << "STL written: propeller.stl\n";
}

