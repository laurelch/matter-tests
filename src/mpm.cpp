// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#include "tools.hpp"
#include "simulation/simulation.hpp"
#include "sampling/sampling_particles.hpp"

#include "objects/object_bump.hpp"
#include "objects/object_gate.hpp"
#include "objects/object_ramp.hpp"
#include "objects/object_box.hpp"
#include "objects/object_plate.hpp"
#include <random>
#include <cmath>
#include <ctime>
#include <chrono>
#include <iomanip>
#include "./data_structures.hpp"




////////////////////////////////////////////////////////
// Setup environment for simulation
////////////////////////////////////////////////////////

void setup_environment_with_bounding_box(Simulation& sim) {
    sim.plates.push_back(std::make_unique<ObjectPlate>(-1.0, PlateType::left,   BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(1.0, PlateType::right,  BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(-1.0, PlateType::back,   BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(1.0, PlateType::front,  BC::NoSlip));

    TV center(0, 1, 0);
    TV half_extents(0.4, 0.4, 0.4);
    TM R = TM::Identity();
    T mytheta = M_PI / 4;
    R(0, 0) = std::cos(mytheta); R(0, 1) = -std::sin(mytheta);
    R(1, 0) = std::sin(mytheta); R(1, 1) =  std::cos(mytheta);

    auto obox = std::make_unique<ObjectBoxRotated>(BC::NoSlip, 0.0, center, half_extents, R);
    sim.objects.push_back(std::move(obox));
}

void setup_environment_without_bounding_box(Simulation& sim) {
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));

    TV center(0, 1, 0);
    TV half_extents(0.4, 0.4, 0.4);
    TM R = TM::Identity();
    T mytheta = M_PI / 4;
    R(0, 0) = std::cos(mytheta); R(0, 1) = -std::sin(mytheta);
    R(1, 0) = std::sin(mytheta); R(1, 1) =  std::cos(mytheta);

    auto obox = std::make_unique<ObjectBoxRotated>(BC::NoSlip, 0.0, center, half_extents, R);
    sim.objects.push_back(std::move(obox));
}

void setup_environment_for_disney(Simulation& sim) {
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));

    TV center(0.5, 0.6828-(2.0/64.0), 0.5);
    TV half_extents(0.2, 0.2, 1);
    TM R = TM::Identity();
    T mytheta = M_PI / 4;
    R(0, 0) = std::cos(mytheta); R(0, 1) = -std::sin(mytheta);
    R(1, 0) = std::sin(mytheta); R(1, 1) =  std::cos(mytheta);

    auto obox = std::make_unique<ObjectBoxRotated>(BC::NoSlip, 0.0, center, half_extents, R);
    sim.objects.push_back(std::move(obox));
}

void setup_environment_only_bounding_box(Simulation& sim) {
    sim.plates.push_back(std::make_unique<ObjectPlate>(-1.0, PlateType::left,   BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(1.0, PlateType::right,  BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(0.0, PlateType::bottom, BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(-1.0, PlateType::back,   BC::NoSlip));
    sim.plates.push_back(std::make_unique<ObjectPlate>(1.0, PlateType::front,  BC::NoSlip));
}

////////////////////////////////////////////////////////
// Setup shared simulation factors
////////////////////////////////////////////////////////

void setup_simulation_with_particles(int particle_count, Simulation& sim, int grid_resolution = 64) {
    T theta_deg = 0;
    T theta = theta_deg * M_PI / 180;
    sim.gravity = TV::Zero();
    sim.gravity[0] = +10.0 * std::sin(theta);
    sim.gravity[1] = -10.0 * std::cos(theta);

    sim.Lx = 1;
    sim.Ly = 1;
    T k_rad;
    switch (particle_count) {
        case 20000: k_rad = 0.031; break;
        case 10000: k_rad = 0.039; break;
        case 5000:  k_rad = 0.050; break;
        case 2000:  k_rad = 0.068; break;
        case 3000:  k_rad = 0.059; break;
        default:    k_rad = 0.05;  break;
    }

    #ifdef THREEDIM
        sim.Lz = 1;
    #endif
    T ppc = 8;
    sim.dx = 1.0 / static_cast<T>(grid_resolution);  // e.g., 1/32 or 1/64
    sim.particle_volume = sim.dx * sim.dx * sim.dx / ppc; // = Lx*Ly*Lz / T(square_samples.size())
    sim.particle_mass = sim.rho * sim.particle_volume;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dist(0.0, 1.0);
    
    sim.Np = particle_count;
    sim.particles = Particles(sim.Np);
    for(int p = 0; p < sim.Np; p++){
        double x = dist(gen)*0.8-0.4;
        double y = dist(gen)*0.4-0.2;
        double z = dist(gen)*0.4-0.2;
        sim.particles.x[p](0)=x;
        sim.particles.x[p](1)=y;
        sim.particles.x[p](2)=z;
    }
    sim.grid_reference_point = TV::Zero();
}

////////////////////////////////////////////////////////
// Benchmark data structure
////////////////////////////////////////////////////////

struct BenchmarkData {
    int num_particles;
    int grid_resolution;  // e.g., 32 for 32^3 grid, 64 for 64^3 grid
    int threads;
    std::vector<int> steps_per_frame;  // Steps at each frame completion
    int total_steps;
    int total_time_ms;
};

////////////////////////////////////////////////////////
// Simulations
////////////////////////////////////////////////////////

int sand_collision(int particle_count) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Simulation sim;
    std::string name = "sand_collision_" + std::to_string(particle_count);
    sim.initialize(true, "output/", name);

    sim.save_grid = true;
    sim.end_frame = 60;
    sim.fps = 30;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    setup_simulation_with_particles(particle_count, sim);

    for (int p = 0; p < sim.Np; p++) {
        sim.particles.x[p][0]  -= 0.5;
        sim.particles.x[p][1]  += 3.0;
        sim.particles.x[p][2]  -= 0.5;
    }

    setup_environment_with_bounding_box(sim);

    sim.plastic_model = PlasticModel::DPVisc;
    sim.use_pradhana = true;
    sim.use_mises_q = false;
    sim.M = std::tan(30*M_PI/180.0);
    sim.q_cohesion = 0;
    sim.perzyna_exp = 1;
    sim.perzyna_visc = 0;

    sim.simulate();

    auto end = high_resolution_clock::now();
    return static_cast<int>(duration_cast<milliseconds>(end - start).count());
}

int bouncy_cube(int particle_count) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Simulation sim;
    std::string name = "bouncy_cube_" + std::to_string(particle_count);
    sim.initialize(true, "output/", name);

    sim.save_grid = true;
    sim.end_frame = 60;
    sim.fps = 30;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = 0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e5;
    sim.nu = 0.45;
    sim.rho = 1000;

    setup_simulation_with_particles(particle_count, sim);

    for (int p = 0; p < sim.Np; p++) {
        sim.particles.x[p][0]  -= 0.5;
        sim.particles.x[p][1]  += 3.0;
        sim.particles.x[p][2]  -= 0.5;
    }

    setup_environment_with_bounding_box(sim);

    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.hardening_law = HardeningLaw::NoHard;
    sim.use_pradhana = false;
    sim.use_mises_q = false;
    sim.M = 0;
    sim.q_cohesion = 0;
    sim.perzyna_exp = 1;
    sim.perzyna_visc = 3000;

    sim.simulate();

    auto end = high_resolution_clock::now();
    return static_cast<int>(duration_cast<milliseconds>(end - start).count());
}

int no_plasticity(int particle_count) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Simulation sim;
    std::string name = "no_plasticity_" + std::to_string(particle_count);
    sim.initialize(true, "output/", name);

    sim.save_grid = true;
    sim.end_frame = 60;
    sim.fps = 30;
    sim.n_threads = 8;
    sim.cfl = 0.5;
    sim.flip_ratio = -0.95;

    sim.elastic_model = ElasticModel::Hencky;
    sim.E = 1e6;
    sim.nu = 0.3;
    sim.rho = 1000;

    setup_simulation_with_particles(particle_count, sim);

    for (int p = 0; p < sim.Np; p++) {
        sim.particles.x[p][0]  -= 0.5;
        sim.particles.x[p][1]  += 3.0;
        sim.particles.x[p][2]  -= 0.5;
    }

    setup_environment_without_bounding_box(sim);

    sim.plastic_model = PlasticModel::NoPlasticity;
    sim.use_pradhana = true;
    sim.use_mises_q = false;
    sim.M = std::tan(30*M_PI/180.0);
    sim.q_cohesion = 0;
    sim.perzyna_exp = 1;
    sim.perzyna_visc = 0;

    sim.simulate();

    auto end = high_resolution_clock::now();
    return static_cast<int>(duration_cast<milliseconds>(end - start).count());
}

// New benchmarking version of snow that returns detailed data
BenchmarkData snow(int particle_count, int end_frame, int threads, int grid_resolution = 64) {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Simulation sim;
    std::string name = "snow_" + std::to_string(particle_count);
    sim.initialize(false, "output/", name);  // Set to false to disable file output

    sim.save_grid = false;  // Disable grid saving for benchmarking
    sim.end_frame = end_frame;
    sim.fps = 24;
    sim.n_threads = threads;
    sim.cfl = 0.1;
    sim.flip_ratio = 0.95;
    sim.reduce_verbose = true;  // Reduce console output

    sim.elastic_model = ElasticModel::Hencky;
    sim.plastic_model = PlasticModel::Snow;

    sim.E = 360;
    sim.nu = 0.3;
    sim.rho = 3;

    setup_simulation_with_particles(particle_count, sim, grid_resolution);

    for (int p = 0; p < sim.Np; p++) {
        sim.particles.x[p][0]  += 0.5;
        sim.particles.x[p][1]  += 1.5-(2.0/static_cast<T>(grid_resolution));
        sim.particles.x[p][2]  += 0.5;
    }

    setup_environment_for_disney(sim);

    sim.simulate();

    auto end = high_resolution_clock::now();
    int duration_ms = static_cast<int>(duration_cast<milliseconds>(end - start).count());
    
    // Collect benchmark data
    BenchmarkData data;
    data.num_particles = sim.Np;
    data.grid_resolution = grid_resolution;
    data.threads = threads;
    data.total_steps = sim.getCurrentTimeStep();
    
    // Get exact steps per frame from simulation
    const auto& frame_steps = sim.getStepsPerFrame();
    data.steps_per_frame.resize(frame_steps.size());
    for (size_t i = 0; i < frame_steps.size(); i++) {
        data.steps_per_frame[i] = static_cast<int>(frame_steps[i]);
    }
    
    data.total_time_ms = duration_ms;
    
    return data;
}

////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////

std::string get_timestamp_filename() {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_now = *std::localtime(&time_t_now);
    
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%y%m%d_%H%M", &tm_now);
    return std::string(buffer) + ".csv";
}

void record_to_csv(const BenchmarkData& data) {
    static std::string csv_filename = get_timestamp_filename();
    
    std::ofstream out(csv_filename, std::ios::app);
    if (!out) {
        std::cerr << "Failed to open " << csv_filename << "\n";
        return;
    }
    
    // Write data: Particles, Grid Size, Threads, Frame 01-20 steps, Total Steps, Time
    out << data.num_particles << "," 
        << data.grid_resolution << ","
        << data.threads;
    
    // Write individual frame steps (up to frame 20)
    int frames_to_write = std::min(20, static_cast<int>(data.steps_per_frame.size()));
    for (int i = 0; i < frames_to_write; i++) {
        out << "," << data.steps_per_frame[i];
    }
    
    // If we have fewer than 20 frames, pad with zeros
    for (int i = frames_to_write; i < 20; i++) {
        out << ",0";
    }
    
    out << "," << data.total_steps << "," 
        << data.total_time_ms << "\n" << std::flush;
}

int main() {
    std::cout << "Starting snow benchmark suite..." << std::endl;
    std::cout << "Results will be saved to: " << get_timestamp_filename() << std::endl;
    std::cout << "Format: Particles, Grid Size, Threads, Steps in Frame 01-20, Total Steps, Explicit C++ T(ms)" << std::endl;
    std::cout << std::endl;

    // 32 Grid Size (1 thread) sweep
    std::cout << "Running 32 grid size (1 thread) sweep..." << std::endl;
    record_to_csv(snow(10000, 20, 1, 32));
    record_to_csv(snow(20000, 20, 1, 32));
    record_to_csv(snow(30000, 20, 1, 32));
    record_to_csv(snow(40000, 20, 1, 32));
    record_to_csv(snow(50000, 20, 1, 32));
    record_to_csv(snow(100000, 20, 1, 32));
    record_to_csv(snow(200000, 20, 1, 32));

    // 64 Grid Size (1 thread) sweep
    std::cout << "Running 64 grid size (1 thread) sweep..." << std::endl;
    record_to_csv(snow(10000, 20, 1, 64));
    record_to_csv(snow(20000, 20, 1, 64));
    record_to_csv(snow(30000, 20, 1, 64));
    record_to_csv(snow(40000, 20, 1, 64));
    record_to_csv(snow(50000, 20, 1, 64));
    record_to_csv(snow(100000, 20, 1, 64));
    record_to_csv(snow(200000, 20, 1, 64));

    // 32 Grid multi-thread sweep
    std::cout << "Running 32 grid size (10k, 50k, 100k, 200k) multi-thread (1, 2, 4, 8, 16, 32, 64) sweep..." << std::endl;
    std::cout << "    Running 32 grid size 10k multi-thread sweep..." << std::endl;
    
    record_to_csv(snow(10000, 20, 2, 32));
    record_to_csv(snow(10000, 20, 4, 32));
    record_to_csv(snow(10000, 20, 8, 32));
    record_to_csv(snow(10000, 20, 16, 32));
    record_to_csv(snow(10000, 20, 32, 32));
    record_to_csv(snow(10000, 20, 64, 32));

    std::cout << "    Running 32 grid size 50k multi-thread sweep..." << std::endl;
    record_to_csv(snow(50000, 20, 1, 32));
    record_to_csv(snow(50000, 20, 2, 32));
    record_to_csv(snow(50000, 20, 4, 32));
    record_to_csv(snow(50000, 20, 8, 32));
    record_to_csv(snow(50000, 20, 16, 32));
    record_to_csv(snow(50000, 20, 32, 32));
    record_to_csv(snow(50000, 20, 64, 32));

    std::cout << "    Running 32 grid size 100k multi-thread sweep..." << std::endl;
    record_to_csv(snow(100000, 20, 1, 32));
    record_to_csv(snow(100000, 20, 2, 32));
    record_to_csv(snow(100000, 20, 4, 32));
    record_to_csv(snow(100000, 20, 8, 32));
    record_to_csv(snow(100000, 20, 16, 32));
    record_to_csv(snow(100000, 20, 32, 32));
    record_to_csv(snow(100000, 20, 64, 32));

    std::cout << "    Running 32 grid size 200k multi-thread sweep..." << std::endl;
    record_to_csv(snow(200000, 20, 1, 32));
    record_to_csv(snow(200000, 20, 2, 32));
    record_to_csv(snow(200000, 20, 4, 32));
    record_to_csv(snow(200000, 20, 8, 32));
    record_to_csv(snow(200000, 20, 16, 32));
    record_to_csv(snow(200000, 20, 32, 32));
    record_to_csv(snow(200000, 20, 64, 32));

    // 64 Grid multi-thread sweep
    std::cout << "Running 64 grid size (10k, 50k, 100k, 200k) multi-thread (1, 2, 4, 8, 16, 32, 64) sweep..." << std::endl;
    std::cout << "    Running 64 grid size 10k multi-thread sweep..." << std::endl;
    record_to_csv(snow(10000, 20, 1, 64));
    record_to_csv(snow(10000, 20, 2, 64));
    record_to_csv(snow(10000, 20, 4, 64));
    record_to_csv(snow(10000, 20, 8, 64));
    record_to_csv(snow(10000, 20, 16, 64));
    record_to_csv(snow(10000, 20, 32, 64));
    record_to_csv(snow(10000, 20, 64, 64));

    std::cout << "    Running 64 grid size 50k multi-thread sweep..." << std::endl;
    record_to_csv(snow(50000, 20, 1, 64));
    record_to_csv(snow(50000, 20, 2, 64));
    record_to_csv(snow(50000, 20, 4, 64));
    record_to_csv(snow(50000, 20, 8, 64));
    record_to_csv(snow(50000, 20, 16, 64));
    record_to_csv(snow(50000, 20, 32, 64));
    record_to_csv(snow(50000, 20, 64, 64));

    std::cout << "    Running 64 grid size 100k multi-thread sweep..." << std::endl;
    record_to_csv(snow(100000, 20, 1, 64));
    record_to_csv(snow(100000, 20, 2, 64));
    record_to_csv(snow(100000, 20, 4, 64));
    record_to_csv(snow(100000, 20, 8, 64));
    record_to_csv(snow(100000, 20, 16, 64));
    record_to_csv(snow(100000, 20, 32, 64));
    record_to_csv(snow(100000, 20, 64, 64));

    std::cout << "    Running 64 grid size 200k multi-thread sweep..." << std::endl;
    record_to_csv(snow(200000, 20, 1, 64));
    record_to_csv(snow(200000, 20, 2, 64));
    record_to_csv(snow(200000, 20, 4, 64));
    record_to_csv(snow(200000, 20, 8, 64));
    record_to_csv(snow(200000, 20, 16, 64));
    record_to_csv(snow(200000, 20, 32, 64));
    record_to_csv(snow(200000, 20, 64, 64));

    std::cout << "\nBenchmark complete! Results saved to: " << get_timestamp_filename() << std::endl;
    return 0;
}
