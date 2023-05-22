#include "sa.hh"

const double def_b0 = 0.0001;
const double def_b1 = 1;
const size_t def_active = 0;
const size_t def_seed = 0;
const size_t def_epochs = 1000;
const size_t def_active_epochs = 1;
const bool def_fixed_beta = false;

std::string gpath;
double b0 = def_b0;
double b1 = def_b1;
size_t active = def_active;
size_t seed = def_seed;
size_t epochs = def_epochs;
size_t active_epochs = def_active_epochs;
bool fixed_beta = def_fixed_beta;
bool output = false;
std::string outpath = "output_log.csv";


void show_help() {
    std::cout << "Partially Frozen Ising Sampler Usage:" << std::endl <<
                 "\t-g|--graph          <string>:    Path to GSET-formatted edgelist file" << std::endl <<
                 "\t-b0|--beta0         <float>:     Starting value for Beta (Default " << def_b0 << ")" << std::endl <<
                 "\t-b1|--beta1         <float>:     Ending value for Beta (Default: 1, will equal Beta0\n\t\t\t\t\t\tif --fixed) (Default " << def_b1 << ")" << std::endl <<
                 "\t--fixed             <flag>:      If enabled, will only run simulation at fixed Beta0\n\t\t\t\t\t\t(overrides any -b1 with Beta1=Beta0) (Default " << def_fixed_beta << ")" << std::endl <<
                 "\t-e|--epochs         <uint64_t>:  Number of epochs to run (AKA ``sweeps'', consider\n\t\t\t\t\t\tflipping each variable once) (Default " << def_epochs << ")" << std::endl <<
                 "\t-a|--active         <size_t>:    Size of the active node set. If 0, then active is\n\t\t\t\t\t\tset to the problem size. (Default " << def_active << ")" << std::endl <<
                 "\t-ae|--active_epochs <size_t>:    Number of epochs to run before shifting the active\n\t\t\t\t\t\tlist. (Default " << def_active_epochs << ")" << std::endl <<
                 "\t-o|--output         <string>:    Output file (CSV formatted) for log info. If not\n\t\t\t\t\t\tprovided, the data will not be logged" << std::endl <<
                 "\t-s|--seed           <size_t>:    Seed for RNG (Default " << def_seed << ")" << std::endl;
}

bool parse_args(int argc, char** argv) {
    for (int i = 0; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            show_help();
            exit(0);
        }else if (arg == "-g" || arg == "--graph") {
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