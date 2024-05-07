#ifndef N_BODY_SIMULATION_H
#define N_BODY_SIMULATION_H

#include <vector>
#include <random>
#include <fstream>

#include "particle.h"

enum class integration_method {
    explicit_euler,
    leapfrog
};

enum class boundary_condition {
    none,
    periodic,
    reflecting
};

template<int dim, typename T,
        boundary_condition Boundary = boundary_condition::none,
        integration_method Integration = integration_method::leapfrog>
class simulation {

    // Compile-time checks on the dimension
    static_assert(dim > 0, "The dimension of the simulation must be greater than 0");
    static_assert(dim < 4, "The dimension of the simulation must be less than 4");

public:

    // Constructor
    // Initialize the simulation with N_particles particles, a time step delta_t, a maximum simulation time time_max,
    // a force softening epsilon, and a list of snapshots times
    simulation(size_t _N_particles, T _delta_t, T _time_max, T _epsilon = 0.01, std::vector<T> snapshots = {}) :
            N_particles{_N_particles},
            particles(_N_particles),
            delta_t{_delta_t},
            time_max{_time_max},
            epsilon{_epsilon},
            snapshots{snapshots} {
        // Allocate random positions and velocities
        std::random_device rd{};
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dis(-1.0, 1.0);

        std::array<T, dim> position{};
        std::array<T, dim> velocity{};

        for (size_t i = 0; i < N_particles; ++i) {
            particles[i].set_mass(1.);

            std::generate(position.begin(), position.end(), [&]() {
                return 0.5 * (dis(gen) + 1.0);
            });
            particles[i].set_position(position);

            std::generate(velocity.begin(), velocity.end(), [&]() {
                return 0.01 * dis(gen);
            });
            particles[i].set_velocity(velocity);
        }
    };

    // Deleting default constructor to avoid uninitialized simulations
    simulation() = delete;

    // Perform the simulation
    // Perform a step using the explicit Euler method
    void run() {
        // Loop over time
        while (time <= time_max) {

            perform_step_on_particles();

            // Add snapshot if the time is bigger than the next snapshot
            if (!snapshots.empty() && time >= snapshots[0]) {
                export_snapshot_to_csv();
                snapshots.erase(snapshots.begin());
            }

            time += delta_t;
        }
    }

private:

    // Particles
    size_t N_particles{};
    std::vector<particle<dim, T>> particles{};

    // Time step
    double delta_t{};
    // Maximum simulation time
    double time_max{};
    // Current simulation time
    double time{};

    // Force softening
    T epsilon{};

    // Snapshots
    std::vector<T> snapshots{};

    // Compute the force between the particles i and j
    void compute_force(size_t i, size_t j) {
        // Compute the force
        auto r_i = particles[i].get_position();
        auto r_j = particles[j].get_position();
        std::array<T, dim> r_ij;
        for (size_t d = 0; d < dim; ++d) {
            r_ij[d] = r_j[d] - r_i[d];
        }
        T r2 = epsilon;
        for (size_t d = 0; d < dim; ++d) {
            r2 += r_ij[d] * r_ij[d];
        }
        T r = std::sqrt(r2);
        T f = particles[i].get_mass() * particles[j].get_mass() / r2;
        for (size_t d = 0; d < dim; ++d) {
            particles[i].set_force(d, particles[i].get_force()[d] + f * r_ij[d] / r);
        }
    }

    // Compute forces between the particle i and all other particles
    void compute_all_forces(size_t i) {
        // Loop over other particles
        for (size_t j = 0; j < N_particles; ++j) {
            if (i == j) {
                continue;
            }
            compute_force(i, j);
        }
    }

    void apply_boundary_conditions(size_t i) {
        // Periodic boundary conditions
        if constexpr (Boundary == boundary_condition::periodic) {
            for (size_t d = 0; d < dim; ++d) {
                if (particles[i].get_position()[d] < 0.0) {
                    particles[i].set_position(d, 1.0 + particles[i].get_position()[d]);
                } else if (particles[i].get_position()[d] > 1.0) {
                    particles[i].set_position(d, particles[i].get_position()[d] - 1.0);
                }
            }
        }
            // Reflecting boundary conditions
        else if constexpr (Boundary == boundary_condition::reflecting) {
            for (size_t d = 0; d < dim; ++d) {
                if (particles[i].get_position()[d] < 0.0) {
                    particles[i].set_position(d, -particles[i].get_position()[d]);
                    particles[i].set_velocity(d, -particles[i].get_velocity()[d]);
                } else if (particles[i].get_position()[d] > 1.0) {
                    particles[i].set_position(d, 2.0 - particles[i].get_position()[d]);
                    particles[i].set_velocity(d, -particles[i].get_velocity()[d]);
                }
            }
        }
    }

    void perform_step_on_particles() {
        // If constexpr is an if at compile-time
        // That is, if Integration == integration_method::explicit_euler then the code inside the if block is compiled
        // Otherwise, the code inside the other block is compiled
        if constexpr (Integration == integration_method::explicit_euler) {
            // Loop over particles
            for (size_t i = 0; i < N_particles; ++i) {
                // Reset the force
                particles[i].reset_force();
                // Compute forces
                compute_all_forces(i);
            }
            // Update the position and velocity
            for (size_t i = 0; i < N_particles; ++i) {
                update_particle(i);
                // Apply boundary conditions
                apply_boundary_conditions(i);
            }
        } else if constexpr (Integration == integration_method::leapfrog) {
            // Update the force at current position
            for (size_t i = 0; i < N_particles; ++i) {
                // Reset the force
                particles[i].reset_force();
                // Compute forces
                compute_all_forces(i);
            }
            // First leapfrog step
            for (size_t i = 0; i < N_particles; ++i) {
                update_particle_leapfrog_first_step(i);
            }

            // Update the force on the half step
            for (size_t i = 0; i < N_particles; ++i) {
                // Reset the force
                particles[i].reset_force();
                // Compute forces
                compute_all_forces(i);
            }

            // Second leapfrog step
            for (size_t i = 0; i < N_particles; ++i) {
                update_particle_leapfrog_second_step(i);
                // Apply boundary conditions
                apply_boundary_conditions(i);
            }
        }
    }

    // Update the position and velocity
    void update_particle(size_t i) {
        auto r = particles[i].get_position();
        auto v = particles[i].get_velocity();
        auto f = particles[i].get_force();

        // Explicit Euler method
        for (size_t d = 0; d < dim; ++d) {
            particles[i].set_position(d, r[d] + delta_t * v[d]);
            particles[i].set_velocity(d, v[d] + delta_t * f[d] / particles[i].get_mass());
        }
    }

    void update_particle_leapfrog_first_step(size_t i) {
        auto r = particles[i].get_position();
        auto v = particles[i].get_velocity();
        auto f = particles[i].get_force();

        std::array<T, dim> v_half;
        for (size_t d = 0; d < dim; ++d) {
            v_half[d] = v[d] + 0.5 * delta_t * f[d] / particles[i].get_mass();
        }
        for (size_t d = 0; d < dim; ++d) {
            particles[i].set_position(d, r[d] + delta_t * v_half[d]);
        }
    }

    void update_particle_leapfrog_second_step(size_t i) {
        auto v = particles[i].get_velocity();
        auto f = particles[i].get_force();

        for (size_t d = 0; d < dim; ++d) {
            particles[i].set_velocity(d, v[d] + 0.5 * delta_t * f[d] / particles[i].get_mass());
        }
    }

    // Print a snapshot on a csv file
    void export_snapshot_to_csv() {
        std::string filename = "snapshot_" + std::to_string(time) + ".csv";
        std::ofstream file(filename);

        // Add information about the simulation in the header
        file << "# N_particles: " << N_particles << "\n";
        file << "# delta_t: " << delta_t << "\n";
        file << "# snapshot time: " << time << "\n";

        file << "# boundary condition: ";
        if constexpr (Boundary == boundary_condition::none) {
            file << "none\n";
        } else if constexpr (Boundary == boundary_condition::periodic) {
            file << "periodic\n";
        } else if constexpr (Boundary == boundary_condition::reflecting) {
            file << "reflecting\n";
        }

        file << "# integration method: ";
        if constexpr (Integration == integration_method::explicit_euler) {
            file << "explicit Euler\n";
        } else if constexpr (Integration == integration_method::leapfrog) {
            file << "leapfrog\n";
        }

        // Adjust the header depending on the dimension
        if constexpr (dim == 2) {
            file << "x,y\n";
        } else if constexpr (dim == 3) {
            file << "x,y,z\n";
        }

        // Write the position of the particles
        for (size_t i = 0; i < N_particles; ++i) {
            auto position = particles[i].get_position();
            for (size_t d = 0; d < dim; ++d) {
                file << position[d];
                if (d != dim - 1) {
                    file << ",";
                }
            }
            file << "\n";
        }
    }

};


#endif //N_BODY_SIMULATION_H
