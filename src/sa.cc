#include "sa.hh"


void MixedSA::update_T() {
    // linearly update
    T += TStep;
};
bool MixedSA::accept(size_t proposed_flip_index){
    return (exp(dE[proposed_flip_index] / T) <= rng(random_gen));
};
void MixedSA::read_graph(std::string gpath){
    std::fstream infile;
    infile.open(gpath, std::fstream::in);
    assert(infile.is_open());
    size_t nodes, edges;
    int offset = 0;
    if (gpath.find(".gset") == std::string::npos) {
        offset = -1;
    }
    std::string line;
    std::istringstream linestream;
    std::getline(infile, line);
    linestream = std::istringstream(line);
    linestream >> nodes >> edges;
    problem_size = nodes;
    if ((10 * edges) < (nodes * nodes / 2)) {
        sparse = true;
        Jsparse = std::vector<std::list<SparseEntry>>(nodes);
    } else {
        sparse = false;
        J = std::vector<double> (nodes*nodes);
    }
    size_t u,v;
    double w;
    while(std::getline(infile, line)) {
        linestream = std::istringstream(line);
        linestream >> u >> v;
        if (linestream.eof()) {
            w = 1.0;
        } else {
            linestream >> w;
        }
        u += offset;
        v += offset;
        if (sparse) {
            Jsparse[u].push_back({v, w});
            Jsparse[v].push_back({u, w});
        } else {
            ACC2D(J, u, v, nodes) = w;
            ACC2D(J, v, u, nodes) = w;
        }
        // initialize dE
        dE[u] += 2*state[u]*state[v]*w;
        dE[v] += 2*state[u]*state[v]*w;
    }
};
void MixedSA::partition(){
    // Partition the instance into disjoint subsets of vertices
    // number of nodes for each regular partition
    size_t partition_size_regular = static_cast<size_t>(
        std::floor(static_cast<double>(problem_size) / partition_count));
    // number of nodes for remaining partition
    size_t partition_size_last = problem_size - partition_size_regular;
    // set equal to regular if the partition count was a factor of the problem size
    partition_size_last = (!partition_size_last) ? partition_size_regular : partition_size_last;
    size_t index = 0;
    std::vector<size_t> allnodes (problem_size);
    for (size_t i = 0; i < problem_size; i++) {
        allnodes[i] = i;
    }
    std::shuffle(allnodes.begin(), allnodes.end(), random_gen);
    // for each partition that isn't the last
    for (size_t i = 0; i < partition_count - 1; i++) {
        partition_sizes[i] = partition_size_regular;
        partitions[i].reserve(partition_size_regular);
        for (size_t s = 0; s < partition_size_regular; s++, index++) {
            membership[allnodes[index]] = i;
            partitions[i].push_back(allnodes[index]);
        } 
    }
    // fill the last partition
    size_t last_index = partition_count - 1;
    partition_sizes[last_index] = partition_size_last;
    partitions[last_index].reserve(partition_size_last);
    for (size_t s = 0; s < partition_size_last; s++, index++) {
        membership[allnodes[index]] = last_index;
        partitions[last_index].push_back(allnodes[index]);
    } 
    assert(index == problem_size);
    
};
void MixedSA::flip(size_t index) {
    int8_t si = -state[index];
    if (sparse) {
        for (SparseEntry e: Jsparse[index]) {
            dE[e.u] += 4 * si * state[e.u] * e.w;
        }
    } else {
        int8_t si = -state[index];
        for (int j = 0; j < problem_size; j++) {
            dE[j] += 4 * si * state[j] * ACC2D(J, index, j, problem_size);
        }
    }
    dE[index] = -dE[index];
    state[index] = si;
};

void MixedSA::synchonize() {
    for (size_t i = 0; i < partition_count; i++) {
        for (size_t index: partitions[i]) {
            state[index] = local_states[i][index];
        }
    }
    if (!all_visible) {
        for (size_t i = 0; i < partition_count; i++) {
            local_states[i] = state;
        }
    }
};
 // calculate energy for whole problem
double MixedSA::energy(){
    double ene = 0.0;
    if (sparse) {
        for (size_t i = 0; i < problem_size; i++) {
            int8_t si = state[i];
            // TODO implement partial information
                for (auto j : Jsparse[i]) {
                    ene += j.w * state[j.u] * si;
                }
        }
        ene /= 2;
    }
    else { 
        for (size_t i = 0; i < problem_size; i++) {
            int8_t si = state[i];
            for (size_t j = i+1; j < problem_size; j++) {
                ene += ACC2D(J, i, j, problem_size) * state[j] * si;
            }
        }
    }
    return ene;
};
 // calculate energy (-1 for whole problem)
double MixedSA::energy(int64_t partition){

    std::vector<int8_t>& statevector = local_states[partition];
    std::unordered_set<size_t>& visible = visible_sets[partition];
    std::vector<size_t>& nodelist = partitions[partition];
    double ene = 0.0;
    if (sparse) {
        for (size_t i : nodelist) {
            int8_t si = statevector[i];
            // TODO implement partial information
                for (auto j : Jsparse[i]) {

                    if ((i != j.u) && (!all_visible || visible.find(j.u) != visible.end())) {
                        ene += j.w * si * statevector[j.u];
                    }
                }
        }
        ene /= 2;
    }else{
        for (size_t i : nodelist) {
            int8_t si = statevector[i];
            for (auto j : visible) {
                if (i != j) {
                    ene += ACC2D(J, i, j, problem_size) * statevector[j] * si;
                }
            }
        }
    }
    return ene;
    
};
// anneal, return the energy found at the end
double MixedSA::anneal(){
    size_t log_epoch = 0;
    for (size_t sync_epoch = 0; sync_epoch < synch_count; sync_epoch++) {
        for (size_t part = 0; part < partition_count; part++) {
            std::uniform_int_distribution<size_t> index_sampler(0, partition_sizes[part]) ;
            T = TSync;
            for (size_t epoch = 0; epoch < partition_epochs; epoch++, log_epoch++) {
                size_t index = index_sampler(random_gen);
                if (accept(index)) {
                    flip(index);
                }
                update_T();
            }
        }
        synchonize();
    }
}; 