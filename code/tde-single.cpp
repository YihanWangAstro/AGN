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
// Q_max: max closest approach
double calc_max_impact_parameter(double Q_max, double v_inf, double M_tot) {
    return Q_max * sqrt(1 + 2 * M_tot / v_inf / v_inf / Q_max);
}

void job(size_t thread_id, size_t scattering_num, std::string fname, bool retro, double a_smbh_r) {
    double v_inf = 50_kms;
    double r_start = 10_AU;

    std::fstream post_flyby_file("single-" + std::to_string(int(MBH3)) + "-" + fname + std::to_string(a_smbh_r) + "-" +
                                     std::to_string(thread_id) + ".txt",
                                 std::ios::out);

    double a_smbh = a_smbh_r * pow(1e8 / 60.0, 1.0 / 3) * a_bh;

    print(post_flyby_file, std::setprecision(16));

    print(std::cout, "starting job on thread ", thread_id, " \n");

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle SMBH{1e8_Ms};
        Particle BH1{30_Ms};

        Particle BH3{MBH3};

        auto smbh_orb = Elliptic(SMBH.mass, M_tot(BH1, BH3), a_smbh, 0.0, 0.0, 0.0, 0.0, -consts::pi * 0.5);

        double b_max = calc_max_impact_parameter(a_bh * 4, v_inf, M_tot(BH1, BH2, BH3));

        double eps = 1e-6;

        auto b = random::Uniform(eps, b_max);

        auto incident_orb =
            Hyperbolic(BH1.mass, BH3.mass, v_inf, b, consts::pi * double(b < 0), 0.0, 0.0, r_start, orbit::Hyper::in);

        move_particles(incident_orb, BH3);

        move_particles(smbh_orb, BH1, BH3);

        move_to_COM_frame(SMBH, BH1, BH3);

        double scattering_t_end = 2 * time_to_periapsis(incident_orb);

        Solver sim{0, SMBH, BH1, BH3};

        Vector v0 = BH3.vel - SMBH.vel;

        Solver::RunArgs args;

        args.atol = 1e-13;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            auto dr = ptc.pos(2) - ptc.pos(0);

            auto dv = ptc.vel(2) - ptc.vel(0);

            double cos1 = dot(dr, dv) / (norm(dr) * norm(dv));

            double cos2 = dot(v0, dv) / (norm(v0) * norm(dv));

            print(post_flyby_file, b, ',', cos1, ',', cos2, ',', dr, ',', dv, '\n');
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
        executor.silent_async(job, i, n / job_num, "smbh_pro", false, 10);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_retro", true, 10);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_pro", false, 5);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_retro", true, 5);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_pro", false, 3);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_retro", true, 3);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_pro", false, 2);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_retro", true, 2);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_pro", false, 1);
    }
    executor.wait_for_all();

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, "smbh_retro", true, 1);
    }
    executor.wait_for_all();
    return 0;
}
