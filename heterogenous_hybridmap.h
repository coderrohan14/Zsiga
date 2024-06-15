#include <iostream>
#include <vector>
#include <unordered_map>
#include "heterogenous_split.h"

class HybridMap {
public:
    int num_nodes;
    int pid; // Processor ID
    HeterogenousSplitter* splitter;

    std::vector<int> local_nodes; // Indices of nodes that are local to this processor
    std::unordered_map<int, int> remote_nodes; // Maps remote nodes to their local indices

    HybridMap(int num_nodes, int pid, HeterogenousSplitter* splitter) : num_nodes(num_nodes), pid(pid), splitter(splitter) {
        distributeNodes();
    }

    void distributeNodes() {
        std::vector<int> counts(splitter->num_procs, 0);
        int total_nodes = num_nodes;

        // First, distribute nodes based on the capabilities proportionally
        for (int i = 0; i < splitter->num_procs; ++i) {
            float proportion = splitter->node_capabilities[i].weighted_sum(splitter->cpu_weight, splitter->memory_weight, splitter->bandwidth_weight) / splitter->total_power;
            counts[i] = static_cast<int>(proportion * total_nodes);
            total_nodes -= counts[i];
        }

        // If there are leftover nodes due to rounding, distribute them starting from the most capable processor
        int i = 0;
        while (total_nodes > 0) {
            counts[i % splitter->num_procs]++;
            total_nodes--;
            i++;
        }

        // Now assign nodes to local or remote depending on the pid and count
        int node_id = 0;
        for (int proc_id = 0; proc_id < splitter->num_procs; ++proc_id) {
            for (int j = 0; j < counts[proc_id]; ++j, ++node_id) {
                if (proc_id == pid) {
                    local_nodes.push_back(node_id);
                } else {
                    remote_nodes[node_id] = -1; // Placeholder for remote nodes
                }
            }
        }
    }

    int& operator[](int node) {
        if (splitter->get_pid_for_node(node) == pid) {
            return local_nodes[node];
        } else {
            return remote_nodes[node];
        }
    }
};