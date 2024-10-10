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

    HybridMap(int num_nodes, int pid, HeterogeneousSplitter* splitter)
    : num_nodes(num_nodes), pid(pid), num_procs(splitter->num_procs), splitter(splitter) {
        std::vector<int> distribution = splitter->getNodeDistribution();
        for(int i=0; i<num_procs; i++){
            std::cout << "Process " << i << " has " << distribution[i] << " nodes." << std::endl;
        }
        num_local_nodes = distribution[pid];  // Assigning the number of local nodes
        p_in = new int[num_local_nodes];
        int global_node_id = 0;
        for (int i = 0; i < num_local_nodes; i++) {
            while (splitter->get_pid_for_node(global_node_id) != pid) {
                global_node_id++;
            }
            p_in[i] = global_node_id;
            // std::cout<<"Process "<<pid<<": Node "<<global_node_id<<" is local"<<std::endl;
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
            //TODO: Optimize this
            int local_index = 0;
            for (int i = 0; i < u; i++) {
                if (splitter->get_pid_for_node(i) == pid) {
                    local_index++;
                }
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
