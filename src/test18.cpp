// propeller_bemt_full.cpp
// Build: g++ -std=c++17 -O2 propeller_bemt_full.cpp -o propeller_bemt_full
#include <cmath>
#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <complex>

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI/180.0;
constexpr double RAD2DEG = 180.0/PI;

using Vec = std::vector<double>;
using Vec3 = std::array<double,3>;

// ---------------------- Airfoil Model ----------------------
struct Airfoil {
    static double Cl(double alpha) {
        // thin airfoil approx, alpha in rad
        return 2.0 * PI * alpha;
    }
    static double Cd(double alpha) {
        return 0.008 + 0.02*alpha*alpha;
    }
};

// ---------------------- Joukowsky Airfoil ----------------------
void genJoukowsky(double disX, double disY, double backFat, int resolution, Vec& px, Vec& py) {
    px.clear(); py.clear(); px.reserve(resolution); py.reserve(resolution);
    double r = 1.0;
    double shape = 1.0;
    std::complex<double> center(disX, disY);

    double extraR = std::sqrt(std::pow(r - std::abs(center.real()),2.0) + std::pow(center.imag(),2.0));

    for(int i=0;i<resolution;++i){
        double theta = -PI + 2.0*PI*i/(resolution);
        std::complex<double> z(extraR*cos(theta)/backFat + center.real(),
                               extraR*sin(theta)/backFat + center.imag());
        std::complex<double> air = z + std::pow(shape,2)/z;
        double X = 0.25*(air.real() + 2.0);
        double Y = 0.25*air.imag();
        px.push_back(X); py.push_back(Y);
    }
}

// ---------------------- Blade Geometry ----------------------
std::vector<std::vector<Vec3>> buildBladeGeometry(
        const Vec& r, const Vec& chord, const Vec& twist,
        const std::vector<std::array<double,2>>& baseAirfoil,
        double offset_root_frac, double offset_tip_frac)
{
    int NR = r.size();
    int NP = baseAirfoil.size();
    std::vector<std::vector<Vec3>> rings(NR,std::vector<Vec3>(NP));

    for(int i=0;i<NR;++i){
        double c = chord[i];
        double t = -twist[i];
        double ct = std::cos(t), st = std::sin(t);

        double frac = (r[i]-r.front())/(r.back()-r.front());
        double cof = offset_root_frac + (offset_tip_frac - offset_root_frac)*frac;
        double xshift = cof*c;

        for(int j=0;j<NP;++j){
            double x = baseAirfoil[j][0]*c - xshift;
            double y = baseAirfoil[j][1]*c;

            double xr = x*ct - y*st;
            double yr = x*st + y*ct;

            double theta = 0.0;
            double cr = std::cos(theta), sr = std::sin(theta);

            Vec3 e_r{cr,sr,0.0}, e_t{-sr,cr,0.0}, e_z{0.0,0.0,1.0};
            Vec3 P{ r[i]*e_r[0] + xr*e_t[0] + yr*e_z[0],
                    r[i]*e_r[1] + xr*e_t[1] + yr*e_z[1],
                    r[i]*e_r[2] + xr*e_t[2] + yr*e_z[2] };
            rings[i][j] = P;
        }
    }
    return rings;
}

// ---------------------- STL Writing ----------------------
void writeSTL(const std::string& filename,const std::vector<std::vector<Vec3>>& bladeRings){
    std::ofstream out(filename);
    out << "solid blade\n";
    int NR = bladeRings.size();
    int NP = bladeRings[0].size();

    auto writeTri=[&](Vec3 a,Vec3 b,Vec3 c){
        Vec3 u{b[0]-a[0],b[1]-a[1],b[2]-a[2]};
        Vec3 v{c[0]-a[0],c[1]-a[1],c[2]-a[2]};
        Vec3 n{u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]};
        out<<"  facet normal "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"\n";
        out<<"    outer loop\n";
        out<<"      vertex "<<a[0]<<" "<<a[1]<<" "<<a[2]<<"\n";
        out<<"      vertex "<<b[0]<<" "<<b[1]<<" "<<b[2]<<"\n";
        out<<"      vertex "<<c[0]<<" "<<c[1]<<" "<<c[2]<<"\n";
        out<<"    endloop\n  endfacet\n";
    };

    for(int i=0;i<NR-1;++i){
        for(int j=0;j<NP;++j){
            int jn = (j+1)%NP;
            writeTri(bladeRings[i][j],bladeRings[i+1][j],bladeRings[i+1][jn]);
            writeTri(bladeRings[i][j],bladeRings[i+1][jn],bladeRings[i][jn]);
        }
    }
    out << "endsolid blade\n";
    out.close();
}

// ---------------------- BEMT Chord & Twist Solver ----------------------
struct BladeSection {
    double r;
    double chord;
    double twist;
};

// ---------------------- BEMT Chord & Twist Solver (Revised) ----------------------
void computeChordTwistBEMT(
        double R, double Rhub, int B, double T, double rho, double rpm, int N,
        Vec& r_out, Vec& chord_out, Vec& alpha_out)
{
    double omega = rpm * 2.0 * PI / 60.0;  // rad/s
    r_out.resize(N); chord_out.resize(N); alpha_out.resize(N);

    // ---------------- Radial stations ----------------
    for(int i=0; i<N; ++i){
        double mu = (i + 0.5)/N;
        r_out[i] = Rhub + mu * (R - Rhub);
    }

    // ---------------- Compute total lift per blade ----------------
    // Linear lift distribution from hub (0) to tip (max)
    double L_total_per_blade = T / B;

    // maximum sectional lift at tip
    double Lprime_tip = 2.0 * L_total_per_blade / (R - Rhub); // rough scaling

    // ---------------- Compute chord and twist per station ----------------
    for(int i=0; i<N; ++i){
        double r = r_out[i];
        double frac = (r - Rhub) / (R - Rhub); // 0..1 along blade

        // ---------------- Lift per unit span ----------------
        // Linear from near zero at hub to max at tip
        double Lprime_r = Lprime_tip * frac;

        // Small relative velocity due to rotation
        double Vrel = omega * r;       // assume hover, Vax ~ 0
        double phi = 0.0;              // initial guess
        double theta = 20.0*DEG2RAD;   // geometric pitch initial guess

        // ---------------- Compute local chord ----------------
        double alpha = theta - phi;    // rad
        double Cl = Airfoil::Cl(alpha);
        double c = (2.0*Lprime_r/(rho*Vrel*Vrel*Cl));

        // ---------------- Clamp chord ----------------

        // ---------------- Compute geometric pitch (twist) ----------------
        // phi = inflow angle, for hover a ~ 0 → phi ~ 0, so geometric pitch ~ angle of attack
        double alpha_desired = 4.0*DEG2RAD; // small positive angle of attack
        double theta_geom = phi + alpha_desired;

        r_out[i] = r;
        chord_out[i] = c;
        alpha_out[i] = theta_geom*RAD2DEG;
    }
}

// ---------------------- MAIN ----------------------
int main(){
    // ---------------------- User Inputs ----------------------
    double R = 0.12;
    double Rhub = 0.04;
    int B = 3;
    double T = 36.0;
    double rho = 1.225;
    double rpm = 6500;
    int N = 40;

    double disX = 0.04, disY = 0.20, backFat=1.12;
    double offset_root_frac = 0.0, offset_tip_frac = -0.15;

    // ---------------------- Airfoil ----------------------
    std::vector<double> px, py;
    genJoukowsky(disX, disY, backFat, 80, px, py);
    std::vector<std::array<double,2>> airfoilPts;
    for(size_t i=0;i<px.size();++i) airfoilPts.push_back({px[i],py[i]});

    // ---------------------- BEMT ----------------------
    std::vector<double> r_stations, chord_stations, alpha_stations;
    computeChordTwistBEMT(R,Rhub,B,T,rho,rpm,N,r_stations,chord_stations,alpha_stations);

    // ---------------------- Build STL ----------------------
    Vec twist_rad(alpha_stations.size());
    for(size_t i=0;i<alpha_stations.size();++i) twist_rad[i] = alpha_stations[i]*DEG2RAD;

    auto rings = buildBladeGeometry(r_stations,chord_stations,twist_rad,airfoilPts,
                                    offset_root_frac,offset_tip_frac);
    writeSTL("blade.stl",rings);

    // ---------------------- Output Data ----------------------
    std::ofstream ofs("constant_gamma_twist.dat");
    ofs<<std::setprecision(10)<<"# r(m) alpha_deg chord(m)\n";
    for(int i=0;i<N;++i)
        ofs<<r_stations[i]<<" "<<alpha_stations[i]<<" "<<chord_stations[i]<<"\n";
    ofs.close();

    std::cout<<"✔ Blade STL generated: blade.stl\n";
    std::cout<<"✔ Chord & twist data: constant_gamma_twist.dat\n";
    return 0;
}

