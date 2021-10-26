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

double MBH3 = 30;
// Q_max: max closest approach
double calc_max_impact_parameter(double Q_max, double v_inf, double M_tot) {
    return Q_max * sqrt(1 + 2 * M_tot / v_inf / v_inf / Q_max);
}

void job(size_t thread_id, size_t scattering_num, double a_bh, std::string fname, bool retro, double a_smbh_r) {
    double v_inf = 50_kms;
    double r_start = 20 * a_bh;

    std::fstream post_flyby_file("ae-" + std::to_string(int(MBH3)) + "-" + fname + std::to_string(a_smbh_r) + "-" +
                                     std::to_string(a_bh) + "-" + std::to_string(thread_id) + ".txt",
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

        auto b = random::Uniform(-b_max, b_max);

        double eps = 1e-6;

        auto incident_orb = Hyperbolic(BH1.mass + BH2.mass, BH3.mass, v_inf, fabs(b) + eps, consts::pi * double(b < 0),
                                       0.0, 0.0, r_start, orbit::Hyper::in);

        move_particles(bh_orb, BH2);

        move_to_COM_frame(BH1, BH2);

        move_particles(incident_orb, BH3);

        move_particles(smbh_orb, BH1, BH2, BH3);

        move_to_COM_frame(SMBH, BH1, BH2, BH3);

        double scattering_t_end = 4 * time_to_periapsis(incident_orb);

        Solver sim{0, SMBH, BH1, BH2, BH3};

        Solver::RunArgs args;

        args.atol = 1e-13;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            auto [a, e] = calc_a_e(ptc.mass(1) + ptc.mass(2), ptc.pos(1) - ptc.pos(2), ptc.vel(1) - ptc.vel(2));

            print(post_flyby_file, b, ',', phi, ',', a / a_bh, ',', e, '\n');
        });

        // args.add_start_point_operation([&](auto& ptc, auto h) { print(std::cout, b, ',', retro, ',', ptc, '\n'); });

        sim.run(args);
    }
}

int main(int argc, char** argv) {
    size_t n = 1000000;
    size_t job_num = 40;

    tf::Executor executor;

    MBH3 = std::stod(argv[1]);

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "smbh_pro", false, 16);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "smbh_retro", true, 16);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "smbh_pro", false, 8);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "smbh_retro", true, 8);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "smbh_pro", false, 4);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "smbh_retro", true, 4);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "smbh_pro", false, 2);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "smbh_retro", true, 2);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "smbh_pro", false, 16);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "smbh_retro", true, 16);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "smbh_pro", false, 8);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "smbh_retro", true, 8);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "smbh_pro", false, 4);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "smbh_retro", true, 4);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "smbh_pro", false, 2);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "smbh_retro", true, 2);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "smbh_pro", false, 16);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "smbh_retro", true, 16);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "smbh_pro", false, 8);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "smbh_retro", true, 8);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "smbh_pro", false, 4);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "smbh_retro", true, 4);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "smbh_pro", false, 2);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "smbh_retro", true, 2);
    }
    executor.wait_for_all();

    return 0;
}
