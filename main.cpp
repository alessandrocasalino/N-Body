#include <iostream>

#include "src/simulation.h"

int main() {

    simulation<2, double,
            boundary_condition::periodic,
            integration_method::leapfrog>
            sim(1000, 0.01, 10, 0.01, {0.0, 1.0, 2.0, 3.0, 4.0});

    sim.run();

    return 0;
}
