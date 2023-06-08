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
#include <unordered_map>

#ifndef SA
#define SA
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
        for (size_t i = 0; i < vec.size() - 1; i++) {
            o << vec[i] << ", ";
        }
        o << vec.back();
    }
    o << "]";
    return o;
}
template <typename T>
std::ostream& operator<<(std::ostream& o, std::unordered_set<T> set) {
    o << "{";
    size_t setsize = set.size();
    if (setsize > 0) {
        size_t ind = 0;
        for (T i : set) {
            o << i;
            if (ind != setsize - 1) {
                o << ", ";
            }
            ind++;
        }
    }
    o << "}";
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
struct SpinSample {
    size_t epoch;
    double beta;
    double ene;
    std::vector<int8_t> spins;
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
void topological_sort_sparse(const std::vector<std::list<SparseEntry>>& graph, std::vector<size_t>& out, std::mt19937& rng);
void topological_sort_dense(const std::vector<double>& graph, std::vector<size_t>& out, std::mt19937& rng);
#define ACC2D(vec, i, j, size) (vec).at((i)*(size) + (j))
class MixedSA {
    private:
    /**
     * Private parameters needed:
     * 1. The number of annealing sweeps per partition
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
        size_t sweeps;
        size_t active_epochs;
        size_t beta_epochs = 1; // number of sweeps to run at fixed beta
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
        std::vector<SpinSample> samples;
        std::vector<std::vector<size_t>> block_indices;

        bool collect_samples;
        size_t sweep_index;
        bool reshuffle = false;
        bool block = false; // whether to use sliding window or blocking
        size_t block_index;
        void update_T();
        bool accept(size_t proposed_flip_index);
        void read_graph(std::string gpath);
        void get_next_active();
        void get_next_block();
        void set_next_fixed();
        void flip(size_t index);
        void synchonize();
        void set_bias();
        double _M();

        


    public:
        MixedSA(): problem_size(0),
                   sweeps(0),
                   Beta0(0),
                   Beta1(0),
                   Beta(0),
                   next_active(0),
                   next_fixed(0),
                   all_visible(false),
                   all_active(false) {};
        MixedSA(std::string gpath, 
                double _Beta0, 
                double _Beta1, 
                size_t _epochs, 
                size_t _active_epochs,
                size_t _active = 0, 
                size_t seed = 0, 
                Traversal _type = Traversal::Sequential,
                bool blocking = false):  \
                   traversal_type(_type),
                   problem_size(0),
                   sweeps(_epochs),
                   active_epochs(_active_epochs),
                   Beta0(_Beta0),
                   Beta1(_Beta1),
                   Beta(_Beta0),
                   random_gen(seed),
                   sweep_index(0),
                   block(blocking) {
            // read input data
            read_graph(gpath);
            active_size = (_active == 0) ? problem_size : _active;
            sweep_index = active_size;
            BetaStep = (Beta1-Beta0) / (std::max(std::ceil(static_cast<double>(sweeps) / active_epochs) - 1, 1.0));
            logdata.reserve(sweeps);
            activeset.reserve(active_size);
            fixedset = std::unordered_set<size_t>(problem_size-active_size);
            rng = std::uniform_real_distribution<double>(0, 1);
            biases = std::vector<double>(problem_size);
            if (traversal_type == Traversal::Topological) {
                if (sparse) {
                    topological_sort_sparse(Jsparse, traversal_order, random_gen);
                } else {
                    topological_sort_dense(J, traversal_order, random_gen);
                }
            }else {
                traversal_order = std::vector<size_t>(problem_size);
                for (size_t i = 0; i < problem_size; i++) {
                    traversal_order[i] = i;
                }
                if (traversal_type == Traversal::Random) {
                    std::shuffle(traversal_order.begin(), traversal_order.end(), random_gen);
                }
            }
            if (blocking) {
                size_t numblocks = static_cast<size_t>(
                    std::ceil( static_cast<double>(problem_size) / active_size ) );
                block_indices = std::vector<std::vector<size_t>>(numblocks);
                size_t lastblock = problem_size % (numblocks*active_size) ?  problem_size-((numblocks-1)*active_size) : active_size;
                for (size_t i = 0; i < numblocks-1;i++) {
                    block_indices[i] = std::vector<size_t>(
                        traversal_order.begin() + i*active_size, traversal_order.begin() + (i+1)*active_size);
                }   
                block_indices[numblocks-1] = std::vector<size_t>(
                        traversal_order.end() - lastblock, traversal_order.end());
                block_index = 0;
                activelist = block_indices[0];
            } else {
                activelist = std::vector<size_t>(traversal_order.begin(), traversal_order.begin() + active_size);
                next_active = traversal_order[sweep_index % problem_size];
                active_index = 0;
                activeset = std::unordered_set<size_t>(activelist.begin(), activelist.end());
                fixedset = std::unordered_set<size_t>(traversal_order.begin() + active_size, traversal_order.end());

            }
            
 
           //set_order();
        }
        double energy_active(); // calculate energy (-1 for whole problem)
        double energy(); // calculate energy (-1 for whole problem)
        double anneal(size_t sample_count = 0); // anneal, return the energy found at the end
        void dump_samples(std::string outpath);
        double cut();
        void dumplog(std::string outfile);
        size_t get_flips() const {
            return flips;
        }
        size_t vcount() const {
            return problem_size;
        }

};
#endif