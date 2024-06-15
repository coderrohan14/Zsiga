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
        total_power = std::accumulate(node_capabilities.begin(), node_capabilities.end(), 0.0f,
                                      [this](float sum, const NodeCapabilities& nc) {
                                          return sum + nc.weighted_sum(cpu_weight, memory_weight, bandwidth_weight);
                                      });
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
        float position = static_cast<float>(node) / static_cast<float>(num_procs);
        auto it = std::lower_bound(cumulative_ratios.begin(), cumulative_ratios.end(), position);
        return std::distance(cumulative_ratios.begin(), it);
    }
};