#include "sa.hh"


int main(int argc, char** argv) {
    std::string gpath = argv[1];
    double b0 = 0.0001;
    double b1 = 1;
    size_t active = std::atol(argv[2]);
    size_t seed = std::atol(argv[3]);

    MixedSA annealer = MixedSA(gpath, b0, b1, 10000, 1000, active, seed);
    printf("Initial Energy: %f %f\n", annealer.energy(), annealer.cut());
    annealer.anneal();
    printf("Final Energy: %f %f\n", annealer.energy(), annealer.cut());
    annealer.dumplog("lattice_2d.csv");
}