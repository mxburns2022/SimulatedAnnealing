#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <list>
#include <fstream>
#include <cassert>
#include <sstream>
#include <unordered_set>
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

struct SALog {
    double T; // Temperature
    double m; // average magnetization per spin
    size_t partition; // index of current partition
    size_t sync_count; // number of previous synchronizations
    size_t epoch; // epoch in current sync phase
};
struct SparseEntry {
    size_t u; // index
    double w; // Coupling strength
};
#define ACC2D(vec, i, j, size) (vec)[(i)*(size) + (j)]
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
        // 
        size_t problem_size;
        size_t partition_epochs;
        size_t synch_count;
        size_t partition_count;
        double T0;
        double T1;
        double T;
        double TSync;
        double TStep;
        bool all_visible; // store whether all spins are visible
        bool all_active; // store whether all spins are active.
        bool sparse;
        std::vector<int8_t> state;
        std::vector<double> J; // if J is dense
        std::vector<std::list<SparseEntry>> Jsparse; // if J is sparse
        std::vector<double> dE;
        std::mt19937 random_gen;
        std::uniform_real_distribution<double> rng;
        std::vector<size_t> partition_sizes;
        std::vector<size_t> membership;
        std::vector<std::vector<size_t>> partitions;
        std::vector<std::unordered_set<size_t>> visible_sets; // only used if !all_visible and !all_active
        std::vector<std::vector<int8_t>> local_states; // store variable states for each partition
        std::vector<SALog> logdata;
        

        void update_T();
        bool accept(size_t proposed_flip_index);
        void read_graph(std::string gpath);
        void partition();
        void flip(size_t index);
        void synchonize();

        


    public:
        MixedSA(): problem_size(0),
                   partition_epochs(0),
                   synch_count(0),
                   partition_count(0),
                   T0(0),
                   T1(0),
                   TSync(0),
                   T(0),
                   all_visible(false),
                   all_active(false) {};
        MixedSA(std::string gpath, double _T0, double _T1, size_t _part_count, size_t _epochs, int64_t max_vis = -1, size_t _synch_count = 1, size_t seed = 0):  \
                   problem_size(0),
                   partition_epochs(_epochs),
                   synch_count(_synch_count),
                   partition_count(_part_count),
                   T0(_T0),
                   T1(_T1),
                   TSync(_T0),
                   T(_T0),
                   all_visible(max_vis == -1),
                   all_active(_part_count == 1),
                   random_gen(seed) {
            // read input data
            read_graph(gpath);
            TStep = (T1-T0) / (partition_epochs * synch_count - 1);
            membership = std::vector<size_t>(problem_size);
            dE = std::vector<double>(problem_size, 0.0);
            logdata = std::vector<SALog>(synch_count * partition_epochs * partition_count);
            partitions = std::vector<std::vector<size_t>>(partition_count);
            visible_sets = std::vector<std::unordered_set<size_t>>(partition_count);
            partition_sizes = std::vector<size_t>(partition_count);
            rng = std::uniform_real_distribution<double>(0, 1);
            local_states = std::vector<std::vector<int8_t>>(partition_count);
            state.reserve(problem_size);
            std::bernoulli_distribution bin;
            // randomly initialize state
            for (size_t i = 0; i < problem_size; i++) {
                int8_t qubov = static_cast<int8_t>(bin(random_gen));
                int8_t spin = 2*qubov - 1;
                state.push_back(spin);
            }
            // copy initial state to each local statevector
            for (size_t i = 0; i < partition_count; i++) {
                local_states[i] = std::vector<int8_t>(state.begin(), state.end());
            }
            partition();
        }
        double energy(int64_t partition); // calculate energy (-1 for whole problem)
        double energy(); // calculate energy (-1 for whole problem)
        double anneal(); // anneal, return the energy found at the end

};