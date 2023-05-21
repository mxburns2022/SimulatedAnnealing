#include "sa.hh"

std::string gpath;
double b0 = 0.0001;
double b1 = 1;
size_t active = 0;
size_t seed = 0;
size_t epochs = 1000;
size_t active_epochs = 1;
bool fixed_beta = false;
bool output = false;
std::string outpath = "output_log.csv";

bool parse_args(int argc, char** argv) {
    for (int i = 0; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-g" || arg == "--graph") {
            assert(i != argc-1);
            gpath = argv[++i];
        }else if (arg == "-b0" || arg == "--beta0") {
            assert(i != argc-1);
            b0 = std::atof(argv[++i]);
        }else if (arg == "-b1" || arg == "--beta1") {
            assert(i != argc-1);
            b1 = std::atof(argv[++i]);
        }else if (arg == "--fixed") {
            assert(i != argc-1);
            fixed_beta = true;
        }else if (arg == "-e" || arg == "--epochs") {
            assert(i != argc-1);
            epochs = std::atol(argv[++i]);
        }else if (arg == "-a" || arg == "--active") {
            assert(i != argc-1);
            active = std::atol(argv[++i]);
        }else if (arg == "-ae" || arg == "--active_epochs") {
            assert(i != argc-1);
            active_epochs = std::atol(argv[++i]);
        }else if (arg == "-o" || arg == "--output") {
            assert(i != argc-1);
            output = true;
            outpath = argv[++i];
        }else if (arg == "-s" || arg == "--seed") {
            assert(i != argc-1);
            seed = std::atol(argv[++i]);
        }
    }
    if (fixed_beta) {
        b1 = b0;
    }
    return true;
}

int main(int argc, char** argv) {
    parse_args(argc, argv);
    MixedSA annealer = MixedSA(gpath, b0, b1, epochs, active_epochs, active, seed);
    annealer.anneal();
    printf("N,active,active_epochs,epochs,Beta0,Beta1,ene,flips,seed\n");
    printf("%ld,%ld,%ld,%ld,%e,%e,%e,%ld,%ld\n", annealer.vcount(), active, active_epochs,epochs, b0, b1,annealer.energy(), annealer.get_flips(),seed);
    // printf("Initial Energy: %f %f\n", annealer.energy(), annealer.cut());
    // 
    // printf("Final Energy: %f %f\n", annealer.energy(), annealer.cut());
    if (output) {
        annealer.dumplog(outpath);
    }
}