#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

struct NodeCapabilities {
    float cpu;      // CPU capability in MHz
    float memory;   // Memory capability in GB
    float bandwidth; // Bandwidth capability in Gbps

    float weighted_sum(float cpu_weight, float memory_weight, float bandwidth_weight) const {
        return cpu * cpu_weight + memory * memory_weight + bandwidth * bandwidth_weight;
    }
};

class HeterogenousSplitter {
    
public:
    int num_procs;
    std::vector<NodeCapabilities> node_capabilities;
    float total_power;
    std::vector<float> cumulative_ratios;
    float cpu_weight;
    float memory_weight;
    float bandwidth_weight;

    HeterogenousSplitter(const std::vector<NodeCapabilities>& node_capabilities,
             float cpu_weight = 1.0, float memory_weight = 1.0, float bandwidth_weight = 1.0)
    : num_procs(node_capabilities.size()), node_capabilities(node_capabilities),
      cpu_weight(cpu_weight), memory_weight(memory_weight), bandwidth_weight(bandwidth_weight) {
        calculateTotalPower();
        calculateCumulativeRatios();
    }

    void calculateTotalPower() {
        total_power = 0;
        for (const auto& nc : node_capabilities) {
            total_power += nc.weighted_sum(cpu_weight, memory_weight, bandwidth_weight);
        }
    }

    void calculateCumulativeRatios() {
        cumulative_ratios.resize(num_procs);
        float sum = 0.0f;
        for (int i = 0; i < num_procs; ++i) {
            sum += node_capabilities[i].weighted_sum(cpu_weight, memory_weight, bandwidth_weight) / total_power;
            cumulative_ratios[i] = sum;
        }
    }

    int get_pid_for_node(int node) {
        float position = (node % num_procs) / static_cast<float>(num_procs);
        auto it = std::lower_bound(cumulative_ratios.begin(), cumulative_ratios.end(), position);
        return std::distance(cumulative_ratios.begin(), it);
    }
};

// Usage:
// int main() {
//     std::vector<NodeCapabilities> node_capabilities = {
//         {2.0, 8.0, 10.0},  // Node 1 capabilities
//         {2.5, 16.0, 8.0},  // Node 2 capabilities
//         {3.0, 32.0, 12.0}, // Node 3 capabilities
//         {1.5, 4.0, 6.0},   // Node 4 capabilities
//         {4.0, 64.0, 14.0}, // Node 5 capabilities
//         {3.5, 48.0, 10.0}, // Node 6 capabilities
//         {2.8, 16.0, 7.0},  // Node 7 capabilities
//         {2.2, 8.0, 9.0},   // Node 8 capabilities
//         {1.8, 12.0, 5.0},  // Node 9 capabilities
//         {2.7, 24.0, 11.0}  // Node 10 capabilities
//     };
    
//     Splitter splitter(node_capabilities, 1.0, 1.0, 1.0); // Equal weights for CPU, memory, and bandwidth

//     for (int node = 0; node < 20; ++node) {
//         std::cout << "Node " << node << " is assigned to processor " << splitter.get_pid_for_node(node) << std::endl;
//     }

//     return 0;
// }
