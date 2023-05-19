#include "sa.hh"


int main(int argc, char** argv) {
    std::string gpath = argv[1];
    double t0 = 5;
    double t1 = 0.3;
    size_t active = std::atol(argv[2]);
    size_t seed = std::atol(argv[3]);

    MixedSA annealer = MixedSA(gpath, t0, t1, 1000000, 1000, active, seed);
    printf("Initial Energy: %f %f\n", annealer.energy(), annealer.cut());
    annealer.anneal();
    printf("Final Energy: %f %f\n", annealer.energy(), annealer.cut());
}