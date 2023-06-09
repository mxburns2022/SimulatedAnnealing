#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <string>
#include <list>
#include <fstream>
#include <cassert>
#include <sstream>
#include <unordered_set>
#include <bitset>
/**
 * Experimental implementation of a partially frozen Ising model sampler
 * The energy of a given partition is:
 *  \sum_{i<j, Ci = Cj}J_{ij}*s_i*s_j + \sum_{i<j, Ci \not{=} Cj}J_{ij}*s_i*s_j = 1/2\sum_{i}s_i * (\sum_{j, Ci = Cj}J_{ij}*s_j + J_i)
 *      where J_i = \sum_{j, Ci \not{=} Cj}J_{ij}*s_j is the sum over external couplings with frozen spins.
 * 
 * Only the spins within an active partition can be flipped, and hence will have nonzero covariance
 * After some defined number of epochs, the partitions will be flipped
 * 
 * For now, assume all couplings are visible, but code to support experiments using only **partial information**
 * 
*/
template <typename T>
std::ostream& operator<<(std::ostream& o, std::vector<T> vec) {
    o << "[";
    if (vec.size() > 0) {
        for (size_t i = 0; i < vec.size(); i++) {
            o << vec[i] << ", ";
        }
        o << vec.back();
    }
    o << "]";
    return o;
}
enum class Traversal {
    Topological, // Perform a topological sort on the graph using random seed vertices
    Random, // scramble each time
    Sequential // traverse in variable order (naive)
};
struct SALog {
    double Beta; // Temperature
    double m; // average magnetization per spin
    double ene;
    size_t epoch; // epoch in current sync phase
};
struct SparseEntry {
    size_t u; // index
    double w; // Coupling strength
};
struct Edge {
    size_t u; // index
    size_t v; // index
    double w; // Coupling strength
};
#define ACC2D(vec, i, j, size) (vec).at((i)*(size) + (j))
class MixedSA {
    private:
    /**
     * Private parameters needed:
     * 1. The number of annealing epochs per partition
     * 2. The size of the problem
     * 3. The number of times to run each partition
     * 4. Starting temp
     * 5. Ending temp
     * Private Data Structures Needed
     * 1. Variable state vector
     * 2. Full problem matrix
     * 3. Variable deltaE vector
     * 4. RNG device
     * 5. Partition sizes
     * 6. Partition memberships
     * 7. Visible lists for each partition (naive slice BRIM == all nodes)
     * 8. Partition lists
     * 9. Net magnetization storage
     * 10. Temperature regulation 
    */
        Traversal traversal_type = Traversal::Sequential;
        size_t problem_size;
        size_t active_size;
        size_t epochs;
        size_t active_epochs;
        size_t beta_epochs = 1; // number of epochs to run at fixed beta
        double Beta0;
        double Beta1;
        double Beta;
        double BetaStep;
        size_t flips = 0;
        double M;
        double ene;
        size_t next_active;
        size_t next_fixed;
        size_t active_index;
        bool all_visible; // store whether all spins are visible
        bool all_active; // store whether all spins are active.
        bool sparse;
        std::vector<int8_t> state;
        std::vector<double> J; // if J is dense
        std::vector<std::list<SparseEntry>> Jsparse; // if J is sparse
        std::vector<double> dE;
        std::mt19937 random_gen;
        std::uniform_real_distribution<double> rng;
        std::unordered_set<size_t> activeset;
        std::vector<size_t> activelist;
        std::unordered_set<size_t> fixedset;
        std::vector<size_t> traversal_order;
        std::vector<SALog> logdata;
        std::vector<double> biases;
        std::unordered_map<std::bitset<100>, size_t> sampled_solutions;

        void update_T();
        bool accept(size_t proposed_flip_index);
        void read_graph(std::string gpath);
        void get_next_active();
        void set_next_fixed();
        void flip(size_t index);
        void synchonize();
        void set_bias();
        double _M();

        


    public:
        MixedSA(): problem_size(0),
                   epochs(0),
                   Beta0(0),
                   Beta1(0),
                   Beta(0),
                   next_active(0),
                   next_fixed(0),
                   all_visible(false),
                   all_active(false) {};
        MixedSA(std::string gpath, double _Beta0, double _Beta1, size_t _epochs, size_t _active_epochs, size_t _active = 0, size_t seed = 0):  \
                   problem_size(0),
                   epochs(_epochs),
                   active_epochs(_active_epochs),
                   Beta0(_Beta0),
                   Beta1(_Beta1),
                   Beta(_Beta0),
                   random_gen(seed) {
            // read input data
            read_graph(gpath);
            active_size = (_active == 0) ? problem_size : _active;
            BetaStep = (Beta1-Beta0) / (epochs - 1);
            logdata.reserve(epochs);
            activeset.reserve(active_size);
            fixedset = std::unordered_set<size_t>(problem_size-active_size);
            rng = std::uniform_real_distribution<double>(0, 1);
            traversal_order = std::vector<size_t>(problem_size);
            biases = std::vector<double>(problem_size);

 
           //set_order();
        }
        double energy_active(); // calculate energy (-1 for whole problem)
        double energy(); // calculate energy (-1 for whole problem)
        double anneal(); // anneal, return the energy found at the end
        double cut();
        void dumplog(std::string outfile);
        size_t get_flips() const {
            return flips;
        }
        size_t vcount() const {
            return problem_size;
        }

};