#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
using namespace orbit;
using f = Interactions<NewtonianGrav>;
using Solver = methods::DefaultMethod<f>;
using Particle = Solver::Particle;
using Vector = Particle::Vector;

double V_INF = 50_kms;
// Q_max: max closest approach
double calc_max_impact_parameter(double Q_max, double v_inf, double M_tot) {
    return Q_max * sqrt(1 + 2 * M_tot / v_inf / v_inf / Q_max);
}

void job(size_t scattering_num, double a_bh, bool retro, double a_smbh_r, double M_prime, double q) {
    double v_inf = V_INF;
    double r_start = 10 * a_bh;

    std::fstream post_flyby_file("TDE-" + std::to_string(M_prime) + "-" + std::to_string(q) + "-" +
                                     std::to_string(a_smbh_r) + "-" + std::to_string(a_bh) + "-" +
                                     std::to_string(v_inf / 1_kms) + ".txt",
                                 std::ios::out);

    print(post_flyby_file, std::setprecision(16));

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle SMBH{1e8_Ms};
        Particle BH1{M_prime};
        Particle BH2{M_prime * q};
        Particle star{1_Ms};

        double R_TDE1 = pow(BH1.mass / star.mass, 1.0 / 3) * pow(star.mass, 0.75) * 1_Rs;
        double R_TDE2 = pow(BH2.mass / star.mass, 1.0 / 3) * pow(star.mass, 0.75) * 1_Rs;

        double a_smbh = a_smbh_r * pow(SMBH.mass / (BH1.mass + BH2.mass), 1.0 / 3) * a_bh;

        auto smbh_orb = Elliptic(SMBH.mass, M_tot(BH1, BH2, star), a_smbh, 0.0, 0.0, 0.0, 0.0, -consts::pi * 0.5);

        auto phi = random::Uniform(0, 2 * consts::pi);

        auto bh_orb = Elliptic(BH1.mass, BH2.mass, a_bh, 0.0, consts::pi * double(retro), 0.0, 0.0, phi);

        double b_max = calc_max_impact_parameter(a_bh * 4, v_inf, M_tot(BH1, BH2, star));

        double eps = 1e-6;

        auto b = random::Uniform(eps, b_max);

        auto incident_orb = Hyperbolic(BH1.mass + BH2.mass, star.mass, v_inf, b, consts::pi * double(b < 0), 0.0, 0.0,
                                       r_start, orbit::Hyper::in);

        move_particles(bh_orb, BH2);

        move_to_COM_frame(BH1, BH2);

        move_particles(incident_orb, star);

        move_particles(smbh_orb, BH1, BH2, star);

        move_to_COM_frame(SMBH, BH1, BH2, star);

        double scattering_t_end = 2 * time_to_periapsis(incident_orb);

        Solver sim{0, SMBH, BH1, BH2, star};

        Vector v0 = star.vel - SMBH.vel;

        Solver::RunArgs args;

        args.atol = 1e-13;

        int stat_flag = 0;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_condition([&](auto& ptc, auto h) {
            double d1 = norm(ptc.pos(1) - ptc.pos(3));
            double d2 = norm(ptc.pos(2) - ptc.pos(3));
            if (d1 < R_TDE1) {
                stat_flag = 1;
                return true;
            } else if (d2 < R_TDE2) {
                stat_flag = 2;
                return true;
            } else {
                return false;
            }
        });

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            auto [a, e] = calc_a_e(ptc.mass(1) + ptc.mass(2), ptc.pos(1) - ptc.pos(2), ptc.vel(1) - ptc.vel(2));

            auto dr = ptc.pos(3) - ptc.pos(0);

            auto dv = ptc.vel(3) - ptc.vel(0);

            double cos1 = dot(dr, dv) / (norm(dr) * norm(dv));

            double cos2 = dot(v0, dv) / (norm(v0) * norm(dv));

            print(post_flyby_file, b, ',', phi, ',', a / a_bh, ',', e, ',', cos1, ',', cos2, ',', dr, ',', dv, ',',
                  M_prime, ',', q, ',', stat_flag, '\n');
        });

        sim.run(args);
    }
}

int main(int argc, char** argv) {
    size_t n = 100000;

    tf::Executor executor;

    V_INF = 3000_kms;

    for (double M = 5; M < 60; M += 5) {
        for (double q = 1; q * M >= 5; q -= 0.1) {
            executor.silent_async(job, n, 1_AU, true, 100, M, q);
            executor.silent_async(job, n, 1_AU, false, 100, M, q);
        }
    }

    executor.wait_for_all();

    return 0;
}
