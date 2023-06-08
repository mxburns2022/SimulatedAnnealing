#include "sa.hh"
#include <stack>

void topological_sort_sparse(const std::vector<std::list<SparseEntry>>& graph, std::vector<size_t>& out, std::mt19937& rng) {
    size_t numnodes = graph.size();
    out.clear();
    out.reserve(numnodes);
    int8_t grey = 1;
    int8_t black = 2;
    int8_t white = 0;
    std::vector<int8_t> color(numnodes, 0); // color nodes to ensure that they aren't visited multiple times
    std::stack<size_t> visitation;
    std::vector<size_t> loop_order (numnodes);
    for (size_t i = 0; i < numnodes; i++) {
        loop_order[i] = i;
    }
    size_t numadded = 0;
    std::shuffle(loop_order.begin(), loop_order.end(), rng);
    std::unordered_set<size_t> added;
    for (size_t i : loop_order) {
        if (color[i] == white) {
            visitation.push(i);
            color[i] = grey;
        }
        while(!visitation.empty()) {
            size_t current = visitation.top();
            bool pushed = false;
            std::vector<SparseEntry> outdeg = std::vector<SparseEntry>(graph[current].begin(), graph[current].end());
            std::shuffle(outdeg.begin(), outdeg.end(), rng);
            for (SparseEntry target: outdeg) {
                if (color[target.u] == white) {
                    pushed = true;
                    visitation.push(target.u);
                    color[target.u] = grey;
                }
            }
            if (!pushed) {
                visitation.pop();
                if (color[current] != black) {
                    color[current] = black;
                    assert(added.find(current) == added.end());
                    out.push_back(current);
                    added.insert(current);
                    numadded++;
                }
                if (current == 11)
                    printf("11 Here %ld\n", i);
            }
        }
    }
    assert(numadded == numnodes);
}   
void topological_sort_dense(const std::vector<double>& graph, std::vector<size_t>& out, std::mt19937& rng) {
    size_t numnodes = static_cast<size_t>(std::sqrt(graph.size()));
    out.reserve(numnodes);
    int8_t black = 1;
    int8_t white = 0;
    std::vector<int8_t> color(numnodes, 0); // color nodes to ensure that they aren't visited multiple times
    std::stack<size_t> visitation;
    std::vector<size_t> loop_order (numnodes);
    for (size_t i = 0; i < numnodes; i++) {
        loop_order[i] = i;
    }

    std::vector<size_t> outdeg = std::vector<size_t>(loop_order.begin(), loop_order.end());
    std::shuffle(loop_order.begin(), loop_order.end(), rng);
    for (size_t i : loop_order) {
        if (color[i] == white) {
            visitation.push(i);
        }
        while(!visitation.empty()) {
            size_t current = visitation.top();
            bool pushed = false;
            std::shuffle(outdeg.begin(), outdeg.end(), rng);
            for (size_t target: outdeg) {
                if (std::abs(ACC2D(graph, current, target, numnodes)) > 0 && color[target] == white) {
                    pushed = true;
                    visitation.push(target);
                    color[target] = black;
                }
            }
            if (!pushed) {
                visitation.pop();
                out.push_back(current);
            }
        }
    }
    std::reverse(out.begin(), out.end());
}   