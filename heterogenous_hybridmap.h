#include "ska/flat_hash_map.hpp"
#include "heterogenous_split.h" 
#include <vector>

typedef ska::flat_hash_map<int, int> hash_t;

class HybridMap {
public:
    int pid;
    int num_procs;
    int num_nodes;
    int num_local_nodes;  // Adding this to store the number of local nodes
    int* p_in;  // Array for storing "inside" vertices
    hash_t p_out;  // Hash table for storing "outside" vertices
    HeterogeneousSplitter* splitter;
    std::pair<int, int> node_range;

    HybridMap(int num_nodes, int pid, HeterogeneousSplitter* splitter)
    : num_nodes(num_nodes), pid(pid), num_procs(splitter->num_procs), splitter(splitter) {
        num_local_nodes = splitter->getNodeCountForProcessor(pid);  // Assigning the number of local nodes
        p_in = new int[num_local_nodes];
        node_range = splitter->getNodeRangeForProcessor(pid);
        int global_node_id = node_range.first;
        for(int i = 0; i < num_local_nodes; i++){
            p_in[i] = global_node_id;
            global_node_id++;
        }
    }

    ~HybridMap() {
        delete[] p_in;  // Free the allocated memory for inside vertices
    }

    int& operator[](int u) {
        // Use the splitter to determine the processor responsible for this vertex
        int responsible_pid = splitter->get_pid_for_node(u);

        // If this processor is responsible, determine local index for the vertex
        if (responsible_pid == pid) {
            int local_index = -1;
            int left = 0;
            int right = num_local_nodes - 1;
            while (left <= right) {
                int mid = left + (right - left) / 2;
                int node = node_range.first + mid;
                if (node == u) {
                    local_index = mid;
                    break;
                } else if (node < u) {
                    left = mid + 1;
                } else {
                    right = mid - 1;
                }
            }
            if (local_index == -1) {
                // Handle error: u not found in the range
                throw std::out_of_range("Node not found in local range");
            }
            return p_in[local_index];
        } else {
            // If not responsible, vertex is stored in the hash table
            if(p_out.find(u) == p_out.end()){
                p_out[u] = u;
            }
            return p_out[u];
        }
    }
};
