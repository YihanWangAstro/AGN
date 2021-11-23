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

double MBH3 = 30;
double V_INF = 50_kms;
// Q_max: max closest approach
double calc_max_impact_parameter(double Q_max, double v_inf, double M_tot) {
    return Q_max * sqrt(1 + 2 * M_tot / v_inf / v_inf / Q_max);
}

void job(size_t scattering_num, double a_bh, std::string fname, bool retro, double a_smbh_r) {
    double v_inf = V_INF;
    double r_start = 10 * a_bh;

    std::fstream post_flyby_file("ae-" + std::to_string(int(MBH3)) + "-" + fname + std::to_string(a_smbh_r) + "-" +
                                     std::to_string(a_bh) + ".txt",
                                 std::ios::out);

    double a_smbh = a_smbh_r * pow(1e8 / 60.0, 1.0 / 3) * a_bh;

    print(post_flyby_file, std::setprecision(16));

    print(std::cout, "starting job on thread ", thread_id, " \n");

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle SMBH{1e8_Ms};
        Particle BH1{30_Ms};
        Particle BH2{30_Ms};
        Particle BH3{MBH3};

        auto smbh_orb = Elliptic(SMBH.mass, M_tot(BH1, BH2, BH3), a_smbh, 0.0, 0.0, 0.0, 0.0, -consts::pi * 0.5);

        auto phi = random::Uniform(0, 2 * consts::pi);

        auto bh_orb = Elliptic(BH1.mass, BH2.mass, a_bh, 0.0, consts::pi * double(retro), 0.0, 0.0, phi);

        double b_max = calc_max_impact_parameter(a_bh * 4, v_inf, M_tot(BH1, BH2, BH3));

        double eps = 1e-6;

        auto b = random::Uniform(eps, b_max);

        auto incident_orb = Hyperbolic(BH1.mass + BH2.mass, BH3.mass, v_inf, b, consts::pi * double(b < 0), 0.0, 0.0,
                                       r_start, orbit::Hyper::in);

        move_particles(bh_orb, BH2);

        move_to_COM_frame(BH1, BH2);

        move_particles(incident_orb, BH3);

        move_particles(smbh_orb, BH1, BH2, BH3);

        move_to_COM_frame(SMBH, BH1, BH2, BH3);

        double scattering_t_end = 2 * time_to_periapsis(incident_orb);

        Solver sim{0, SMBH, BH1, BH2, BH3};

        Vector v0 = BH3.vel - SMBH.vel;

        Solver::RunArgs args;

        args.atol = 1e-13;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            auto [a, e] = calc_a_e(ptc.mass(1) + ptc.mass(2), ptc.pos(1) - ptc.pos(2), ptc.vel(1) - ptc.vel(2));

            auto dr = ptc.pos(3) - ptc.pos(0);

            auto dv = ptc.vel(3) - ptc.vel(0);

            double cos1 = dot(dr, dv) / (norm(dr) * norm(dv));

            double cos2 = dot(v0, dv) / (norm(v0) * norm(dv));

            print(post_flyby_file, b, ',', phi, ',', a / a_bh, ',', e, ',', cos1, ',', cos2, ',', dr, ',', dv, '\n');
        });

        // args.add_start_point_operation([&](auto& ptc, auto h) { print(std::cout, b, ',', retro, ',', ptc, '\n'); });

        sim.run(args);
    }
}

int main(int argc, char** argv) {
    size_t n = 100000;

    tf::Executor executor;

    MBH3 = std::stod(argv[1]);
    V_INF = std::stod(argv[2]) * 1_kms;

    executor.silent_async(job, n, 10_AU, "smbh_pro", false, 100);

    executor.silent_async(job, n, 10_AU, "smbh_retro", true, 100);

    executor.silent_async(job, n, 10_AU, "smbh_pro", false, 50);

    executor.silent_async(job, n, 10_AU, "smbh_retro", true, 50);

    executor.silent_async(job, n, 10_AU, "smbh_pro", false, 30);

    executor.silent_async(job, n, 10_AU, "smbh_retro", true, 30);

    executor.silent_async(job, n, 10_AU, "smbh_pro", false, 20);

    executor.silent_async(job, n, 10_AU, "smbh_retro", true, 20);

    executor.silent_async(job, n, 10_AU, "smbh_pro", false, 10);

    executor.silent_async(job, n, 10_AU, "smbh_retro", true, 10);

    executor.silent_async(job, n, 1_AU, "smbh_pro", false, 100);

    executor.silent_async(job, n, 1_AU, "smbh_retro", true, 100);

    executor.silent_async(job, n, 1_AU, "smbh_pro", false, 50);

    executor.silent_async(job, n, 1_AU, "smbh_retro", true, 50);

    executor.silent_async(job, n, 1_AU, "smbh_pro", false, 30);

    executor.silent_async(job, n, 1_AU, "smbh_retro", true, 30);

    executor.silent_async(job, n, 1_AU, "smbh_pro", false, 20);

    executor.silent_async(job, n, 1_AU, "smbh_retro", true, 20);

    executor.silent_async(job, n, 1_AU, "smbh_pro", false, 10);

    executor.silent_async(job, n, 1_AU, "smbh_retro", true, 10);

    executor.silent_async(job, n, 0.1_AU, "smbh_pro", false, 100);

    executor.silent_async(job, n, 0.1_AU, "smbh_retro", true, 100);

    executor.silent_async(job, n, 0.1_AU, "smbh_pro", false, 50);

    executor.silent_async(job, n, 0.1_AU, "smbh_retro", true, 50);

    executor.silent_async(job, n, 0.1_AU, "smbh_pro", false, 30);

    executor.silent_async(job, n, 0.1_AU, "smbh_retro", true, 30);

    executor.silent_async(job, n, 0.1_AU, "smbh_pro", false, 20);

    executor.silent_async(job, n, 0.1_AU, "smbh_retro", true, 20);

    executor.silent_async(job, n, 0.1_AU, "smbh_pro", false, 10);

    executor.silent_async(job, n, 0.1_AU, "smbh_retro", true, 10);

    executor.wait_for_all();

    return 0;
}
