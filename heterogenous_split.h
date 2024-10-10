#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

struct NodeCapabilities {
    float cpu;      // CPU capability in MHz
    float memory;   // Memory capability in GB
    float bandwidth; // Bandwidth capability in Gbps

    float weighted_sum(float cpu_weight, float memory_weight, float bandwidth_weight) const {
        return cpu * cpu_weight + memory * memory_weight + bandwidth * bandwidth_weight;
    }
};

class HeterogeneousSplitter {
    
public:
    int num_procs;
    int total_nodes;
    std::vector<NodeCapabilities> node_capabilities;
    float total_power;
    std::vector<float> cumulative_ratios;
    float cpu_weight;
    float memory_weight;
    float bandwidth_weight;

    HeterogeneousSplitter(const std::vector<NodeCapabilities>& node_capabilities,
                          int total_nodes,
                          float cpu_weight = 1.0, float memory_weight = 1.0, float bandwidth_weight = 1.0)
    : num_procs(node_capabilities.size()), total_nodes(total_nodes), node_capabilities(node_capabilities),
      cpu_weight(cpu_weight), memory_weight(memory_weight), bandwidth_weight(bandwidth_weight) {
        calculateTotalPower();
        calculateCumulativeRatios();
    }

    void calculateTotalPower() {
        total_power = 0.0f;
        for (const auto& nc : node_capabilities) {
            total_power += nc.weighted_sum(cpu_weight, memory_weight, bandwidth_weight);
        }
    }

    void calculateCumulativeRatios() {
        cumulative_ratios.resize(num_procs);
        float sum = 0.0f;
        for (int i = 0; i < num_procs; ++i) {
            sum += node_capabilities[i].weighted_sum(cpu_weight, memory_weight, bandwidth_weight);
            cumulative_ratios[i] = sum / total_power;
        }
    }

    int get_pid_for_node(int node) {
        float position = static_cast<float>(node) / total_nodes;
        auto it = std::lower_bound(cumulative_ratios.begin(), cumulative_ratios.end(), position);
        return std::distance(cumulative_ratios.begin(), it);
    }

    int getMaxNodesPerProcessor() {
        std::vector<int> distribution = this->getNodeDistribution();
        return *std::max_element(distribution.begin(), distribution.end());
    }

    std::vector<int> getNodeDistribution() {
        std::vector<int> distribution(num_procs, 0);
        for (int node = 0; node < total_nodes; ++node) {
            int pid = get_pid_for_node(node);
            distribution[pid]++;
        }
        return distribution;
    }

    std::vector<int> getNodesForProcessor(int processor_id) {
        std::vector<int> nodes;
        for (int node = 0; node < total_nodes; ++node) {
            if (get_pid_for_node(node) == processor_id) {
                nodes.push_back(node);
            }
        }
        return nodes;
    }
};