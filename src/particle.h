#ifndef N_BODY_PARTICLE_H
#define N_BODY_PARTICLE_H

#include <span>
#include <algorithm>

template<int dim, typename T>
class particle {

    // Compile-time checks on the dimension
    static_assert(dim > 0, "The dimension of the particle must be greater than 0");
    static_assert(dim < 4, "The dimension of the particle must be less than 4");

private:
    T mass{};
    T position[dim]{};
    T velocity[dim]{};
    T force[dim]{};

public:

    particle() = default;

    T get_mass() const {
        return mass;
    }

    void set_mass(T _mass) {
        particle::mass = _mass;
    }

    std::span<const T, dim> get_position() const {
        return std::span<const T, dim>(position);
    }

    void set_position(std::span<T, dim> _position) {
        std::copy(_position.begin(), _position.end(), position);
    }

    void set_position(size_t i, T value) {
        position[i] = value;
    }

    std::span<const T, dim> get_velocity() const {
        return std::span<const T, dim>(velocity);
    }

    void set_velocity(std::span<T, dim> _velocity) {
        std::copy(_velocity.begin(), _velocity.end(), velocity);
    }

    void set_velocity(size_t i, T value) {
        velocity[i] = value;
    }

    std::span<const T, dim> get_force() const {
        return std::span<const T, dim>(force);
    }

    void set_force(std::span<T, dim> _force) {
        std::copy(_force.begin(), _force.end(), force);
    }

    void set_force(size_t i, T value) {
        force[i] = value;
    }

    void reset_force() {
        std::span<T, dim> force_span(force);
        std::fill(force_span.begin(), force_span.end(), static_cast<T>(0.));
    }


};


#endif //N_BODY_PARTICLE_H
