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

    // ---------------- Bell curve parameters ----------------
    double r_peak = Rhub + 0.7*(R - Rhub);
    double sigma = 0.4*(R - Rhub);

    // Approximate lift distribution scaling
    double area = 0.0;
    for(int i=0;i<N;++i){
        double r = r_out[i];
        area += std::exp(-0.5*std::pow((r - r_peak)/sigma,2));
    }
    double L_scale = (T / B) / (area * (R - Rhub)/N);  // N/m per blade

    // ---------------- Compute chord and twist per station ----------------
    for(int i=0; i<N; ++i){
        double r = r_out[i];
        double Lprime_r = L_scale * std::exp(-0.5*std::pow((r - r_peak)/sigma,2));

        // Initial guesses
        double a = 0.04, ap = 0.01;
        double theta = 20.0 * DEG2RAD;

        for(int iter=0; iter<250; ++iter){
            double Vax = omega * r * a + 1e-6;
            double Vtan = omega * r * (1.0 - ap);
            double Vrel = std::sqrt(Vax*Vax + Vtan*Vtan);
            double phi = std::atan2(Vax, Vtan);
            double alpha = theta - phi;

            double Cl = Airfoil::Cl(alpha);
            double Cd = Airfoil::Cd(alpha);

            double L = 0.5*rho*Vrel*Vrel*Cl;  // per unit chord
            double D = 0.5*rho*Vrel*Vrel*Cd;

            double Fn = L*std::cos(phi) - D*std::sin(phi);
            double Ft = L*std::sin(phi) + D*std::cos(phi);

            double sigma = B / (2*PI*r) * chord_out[i];
            double a_new = Fn / (4*PI*r*rho*Vrel*Vrel + 1e-6);
            double ap_new = Ft / (4*PI*r*rho*Vrel*Vtan + 1e-6);

            a = 0.7*a + 0.3*a_new;
            ap = 0.7*ap + 0.3*ap_new;
        }

        // Compute chord to match local lift
        double Vax_final = omega * r * a + 1e-6;
        double Vtan_final = omega * r * (1.0 - ap);
        double Vrel_final = std::sqrt(Vax_final*Vax_final + Vtan_final*Vtan_final);
        double phi_final = std::atan2(Vax_final, Vtan_final);

        double alpha_desired = 16.0 * DEG2RAD; // target AoA
        double theta_geom = phi_final + alpha_desired;

        double Cl_final = Airfoil::Cl(alpha_desired);
        double c = 2.0 * Lprime_r / (rho * Vrel_final * Vrel_final * Cl_final + 1e-6);

        chord_out[i] = c;  // m
        alpha_out[i] = theta_geom * RAD2DEG;
    }
}


// ---------------------- MAIN ----------------------
int main(){
    // ---------------------- User Inputs ----------------------
    double R = 0.14;
    double Rhub = 0.04;
    int B = 4;
    double T = 18.0;
    double rho = 1.225;
    double rpm = 6500;
    int N = 80;

    double disX = 0.04, disY = 0.18, backFat=1.08;
    double offset_root_frac = 0.0, offset_tip_frac = -0.1;

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

