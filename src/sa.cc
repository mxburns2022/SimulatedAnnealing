#include "sa.hh"


void MixedSA::update_T() {
    // linearly update
    Beta += BetaStep;
};
bool MixedSA::accept(size_t proposed_flip_index){
    if (dE[proposed_flip_index] < 0) {
        return true;
    }
    return exp(-dE[proposed_flip_index] * Beta) > rng(random_gen);
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
        //printf("Sparse\n");
        sparse = true;
        Jsparse = std::vector<std::list<SparseEntry>>(nodes);
        for (size_t i = 0; i < nodes;i++) {
            Jsparse[i] = std::list<SparseEntry>();
        }
    } else {
        //printf("Dense\n");
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
        dE[u] += -2*state[u]*state[v]*w;
        dE[v] += -2*state[u]*state[v]*w;
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
/**
 * Get the next active spin in a sweeping Window scheme
*/
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
    active_index = ( active_index + 1 ) % active_size;
    if (reshuffle) {
        if (sweep_index == problem_size - 1 && traversal_type == Traversal::Random) {
            std::sample(fixedset.begin(), fixedset.end(), traversal_order.begin(), fixedset.size(), random_gen);        
        } 
        if (sweep_index == 0  && traversal_type == Traversal::Random) {
            std::shuffle(traversal_order.begin() + active_size, traversal_order.end(), random_gen);
        }
    }
    sweep_index = (sweep_index + 1) % problem_size;
    next_active = traversal_order.at(sweep_index);

}
/**
 * Get the next active spin in a sweeping Window scheme
*/
void MixedSA::get_next_block() {    
    block_index = (block_index + 1) % block_indices.size();
    activelist = block_indices[block_index];
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
double MixedSA::_M(){
    double _Mval = 0.0;
    for (int8_t si : state) {
        _Mval += si;
    }
    return _Mval;
}
std::ostream& operator<<(std::ostream& o, std::vector<int8_t> vec) {
    if (vec.size() > 0) {
        for (size_t i = 0; i < vec.size() - 1; i++) {
            o << static_cast<int>(vec[i]) << " ";
        }
        o << static_cast<int>(vec.back());
    }
    return o;
}
void MixedSA::dumplog(std::string outpath) {
    std::fstream outstream;
    outstream.open(outpath, std::fstream::out);
    outstream << "epoch,beta,M,ene" << std::endl;
    for (SALog entry: logdata) {
        outstream << entry.epoch << "," <<
                     entry.Beta << "," <<
                     entry.m    <<  ","<<
                     entry.ene    << std::endl;
    }
    outstream.close();
}

void MixedSA::dump_samples(std::string outpath) {
    std::fstream out;
    out.open(outpath, std::fstream::out);
    assert(out.is_open());
    out << "epoch,beta,ene,state" << std::endl;
    for (auto i: samples) {
        out << i.epoch << "," << 
               i.beta << "," << 
               i.ene << "," << 
               "\"" << i.spins << "\"" << std::endl;
    }
    out.close();
};
// anneal, return the energy found at the end
double MixedSA::anneal(size_t sample_count){
    std::uniform_int_distribution<size_t> index_sampler(0, active_size-1) ;
    size_t sample_epochs = epochs;
    if (sample_count > 0) {
        sample_epochs = (epochs*active_epochs - 1) / sample_count;
        samples.reserve(sample_count);
    }
    set_bias();
    ene = energy();
    M = _M();
    double best_ene = ene;
    std::vector<int8_t> best_state = state;
    flips = 0;
    if (active_size == problem_size) {
        active_epochs = 1;
    }
    size_t sample_counter = 0;
    for (size_t e = 0; e < epochs; e++) {
        for (size_t ae = 0; ae < active_epochs; ae++, sample_counter++) {
            for (size_t index : activelist){
                // size_t indval = index_sampler(random_gen);
                    
                    if (accept(index)) {
                        flips += 1;
                        ene += dE[index];
                        M -= 2*state[index];
        #ifndef NDEBUG
                        state[index] = -state[index];
                        std::cout << index << std::endl;
                        double ene_other = energy();
                        std::cout << ene << " " << dE[index] << " " <<ene_other << " "  << cut() << std::endl;
                        if (ene != ene_other) {
                            exit(-1);
                        }
                        assert(ene == ene_other);
                        state[index] = -state[index];
        #endif
                        flip(index);
                        if (ene < best_ene) {
                            best_ene = ene;
                            best_state = state;
                        }
                    }
                }
            if (sample_counter == sample_epochs) {
                samples.push_back({e+ae, Beta, ene, state});
                sample_counter = 0;
            }
            logdata.push_back({Beta, std::pow(M / problem_size, 2), ene, e + ae});
        }
        update_T();
        if (active_size != problem_size) {
            if (block) {
                get_next_block();
            } else {
                get_next_active();
            }
        }
        //exit(0);
    }
    state = best_state;
    std::cout << cut() << std::endl;
    return best_ene;
}; 