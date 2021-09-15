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

// Q_max: max closest approach
double calc_max_impact_parameter(double Q_max, double v_inf, double M_tot) {
    return Q_max * sqrt(1 + 2 * M_tot / v_inf / v_inf / Q_max);
}

void job(size_t thread_id, size_t scattering_num, double a_bh, std::string fname, bool retro) {
    double v_inf = 50_kms;
    double r_start = 10 * a_bh;
    // double a_bh = 10_AU;

    std::fstream post_flyby_file(fname + std::to_string(a_bh) + "-" + std::to_string(thread_id) + ".txt",
                                 std::ios::out);

    print(post_flyby_file, std::setprecision(16));

    print(std::cout, "starting job on thread ", thread_id, " \n");

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle BH1{30_Ms};
        Particle BH2{30_Ms};
        Particle BH3{30_Ms};

        // create planetary system with two giant planets
        double inc = 0;

        auto phi = random::Uniform(0, 2 * consts::pi);

        auto bh_orb = Elliptic(BH1.mass, BH2.mass, a_bh, 0.0, inc, 0.0, 0.0, phi);

        double b_max = calc_max_impact_parameter(a_bh * 4, v_inf, M_tot(BH1, BH2, BH3));

        auto b = random::Uniform(0, b_max);

        auto incident_orb = Hyperbolic(BH1.mass + BH2.mass, BH3.mass, v_inf, b, consts::pi * double(retro), 0.0, 0.0,
                                       r_start, orbit::Hyper::in);

        move_particles(bh_orb, BH2);

        move_to_COM_frame(BH1, BH2);

        move_particles(incident_orb, BH3);

        double scattering_t_end = 4 * time_to_periapsis(group(BH1, BH2), BH3);

        Solver sim{0, BH1, BH2, BH3};

        Solver::RunArgs args;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            auto [a, e] = calc_a_e(ptc.mass(0) + ptc.mass(1), ptc.pos(0) - ptc.pos(1), ptc.vel(0) - ptc.vel(1));

            int fate = 0;

            if (a > a_bh) {
                fate = 0;
            } else if (a < 0) {
                fate = -1;
            } else {
                fate = 2;
            }

            print(post_flyby_file, b, ',', phi, ',', fate, '\n');
        });

        sim.run(args);
    }
}

int main() {
    size_t n = 1000000;  // total scattering number
    size_t job_num = 40;

    tf::Executor executor;

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "pro", false);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 10_AU, "retro", true);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "pro", false);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 1_AU, "retro", true);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "pro", false);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, 0.1_AU, "retro", true);
    }
    executor.wait_for_all();

    return 0;
}
