#include "sa.hh"


void MixedSA::update_T() {
    // linearly update
    T += TStep;
};
bool MixedSA::accept(size_t proposed_flip_index){
    if (dE[proposed_flip_index] < 0) {
        return true;
    }
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
    state.reserve(problem_size);
    std::bernoulli_distribution bin;
    // randomly initialize state
    for (size_t i = 0; i < problem_size; i++) {
        int8_t qubov = static_cast<int8_t>(bin(random_gen));
        int8_t spin = 2*qubov - 1;
        state.push_back(spin);
    }
    if ((10 * edges) < (nodes * nodes / 2)) {
        printf("Sparse\n");
        sparse = true;
        Jsparse = std::vector<std::list<SparseEntry>>(nodes);
        for (size_t i = 0; i < nodes;i++) {
            Jsparse[i] = std::list<SparseEntry>();
        }
    } else {
        printf("Dense\n");
        sparse = false;
        J = std::vector<double> (nodes*nodes);
    }
    size_t u,v;
    double w;
    problem_size = nodes;
    dE = std::vector<double>(nodes, 0.0);
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
        dE[u] -= 2*state[u]*state[v]*w;
        dE[v] -= 2*state[u]*state[v]*w;
    }
};

void MixedSA::flip(size_t index) {
    int8_t si = -state[index];
    if (sparse) {
        for (SparseEntry e: Jsparse[index]) {
            dE[e.u] -= 4 * si * state[e.u] * e.w;
        }
    } else {
        int8_t si = -state[index];
        for (size_t j = 0; j < problem_size; j++) {
            dE[j] -= 4 * si * state[j] * ACC2D(J, index, j, problem_size);
        }
    }
    dE[index] = -dE[index];
    state[index] = si;
};
void MixedSA::get_next_active() {
    std::vector<SparseEntry> result;

    next_fixed = activelist[active_index];
    assert(fixedset.find(next_active) != fixedset.end());
    fixedset.erase(next_active);
    fixedset.insert(activelist[active_index]);
    assert(activeset.find(next_active) == activeset.end());
    activeset.insert(next_active);
    activeset.erase(activelist[active_index]);
    activelist[active_index] = next_active;
    if (sparse) {
        for (auto e: Jsparse[next_active]) {
            if ((activeset.find(e.u) != activeset.end())) {
                biases[e.u] -= e.w * state[next_active];
            }
        }
        for (auto e: Jsparse[next_fixed]) {
            if ((activeset.find(e.u) != activeset.end())) {
                biases[e.u] += e.w * state[next_fixed];
            }
        }
    } else {
        for (auto j : activeset) {
            biases[j] -= ACC2D(J, j, next_active, problem_size) * state[next_active];
            biases[j] += ACC2D(J, j, next_fixed, problem_size) * state[next_fixed];
        }
    }
    biases[next_active] = 0.0;
    if (sparse) {
        for (auto e: Jsparse[next_active]) {
            if ((activeset.find(e.u) == activeset.end())) {
                biases[next_active] += e.w * state[e.u];
            }
        }
    } else {
        for (auto j: fixedset) {
            biases[next_active] += ACC2D(J, j, next_active, problem_size) * state[j];
        }
    }
    std::cout << biases[next_active] << std::endl;
    active_index = ( active_index + 1 ) % active_size;
     if (traversal_type == Traversal::Sequential) {
        next_active = (next_active + 1) % problem_size;
    } else if (traversal_type == Traversal::Random) {
        std::sample(fixedset.begin(), fixedset.end(), &next_active, 1, random_gen);
    }
}

// void MixedSA::update_order() {
//     if (traversal_type == Traversal::Sequential) {
//         for (size_t i = 0; i < )
//     }
// }
 // calculate energy for whole problem
double MixedSA::energy(){
    double ene = 0.0;
    if (sparse) {
        for (size_t i = 0; i < problem_size; i++) {
            int8_t si = state[i];
            // TODO implement partial information
                for (auto j : Jsparse[i]) {
                    ene += j.w * state.at(j.u) * si;
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

void MixedSA::set_bias() {
    for (size_t i : activelist) {
        biases[i] = 0.0;
        if (sparse) {
            for (auto e : Jsparse[i]) {
                if (activeset.find(e.u) == activeset.end()){

                    biases[i] += e.w * state[e.u];
            }}
        } else {
            for (auto j : fixedset) {
                biases[i] += ACC2D(J, i, j, problem_size) * state[j];
            }
        }
    }
}
double MixedSA::cut(){
    double result = 0.0;
    if (sparse) {
        for (size_t i = 0; i < problem_size; i++) {
            for (auto j: Jsparse[i]) {
                if (i < j.u && state[i] != state[j.u])
                    result += j.w;
            }   
        }
    } else {
        for (size_t i = 0; i < problem_size; i++) {
            for (size_t j = i+1; j < problem_size; j++) {
                if (state[i] != state[j])
                    result += ACC2D(J, i, j, problem_size);
            }   
        }
    }
    return result;
}
// anneal, return the energy found at the end
double MixedSA::anneal(){
    std::uniform_int_distribution<size_t> index_sampler(0, active_size-1) ;
    std::vector<size_t> full_list(problem_size);
    for (size_t i = 0; i < problem_size; i++) {
        full_list[i] = i;
    }
    activelist = std::vector<size_t>(full_list.begin(), full_list.begin() + active_size);
    next_active = active_size;
    active_index = 0;
    activeset = std::unordered_set<size_t>(activelist.begin(), activelist.end());
    fixedset = std::unordered_set<size_t>(full_list.begin() + active_size, full_list.end());
    set_bias();
    double ene = energy();
    double best_ene = ene;
    std::vector<int8_t> best_state = state;
    for (size_t e = 0; e < epochs; e+=active_epochs) {
        for (size_t ae = 0; ae < active_epochs; ae++) {
            size_t indval = index_sampler(random_gen);
            size_t index = activelist.at(indval);
            assert(activeset.find(index) != activeset.end());
            
            assert(activeset.size() == active_size);
            if (accept(index)) {
                ene += dE[index];
#ifndef NDEBUG
                state[index] = -state[index];
                std::cout << index << std::endl;
                std::cout << ene << " " << dE[index] << " " << energy() << std::endl;
                assert(ene == energy());
                state[index] = -state[index];
#endif
                flip(index);
                if (ene < best_ene) {
                    best_ene = ene;
                    best_state = state;
                }
            }
            
            update_T();
            std::cout << T << std::endl;
        }
        if (active_size != problem_size)
            get_next_active();
#ifndef NDEBUG
        std::vector<double> temp = biases;
        set_bias();
        for (auto i : activelist) {
            if (std::abs(temp[i] - biases[i]) > 1e-10) {
                std::cout << i <<
                " " << temp[i] <<
                " " << biases[i] << std::endl;
            }
            assert(std::abs(temp[i] - biases[i]) < 1e-10);
        }
#endif
        //exit(0);
    }
    state = best_state;
    return best_ene;
}; 