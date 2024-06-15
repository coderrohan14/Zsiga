#include "ska/flat_hash_map.hpp"
#include "split.h"
// #include <memory.h>
// using namespace std;

typedef ska::flat_hash_map<int, int> hash_t;

class HybridMap{
    public:
    int num_local_nodes;
    int pid;
    int num_procs;
    int num_slow_procs;
    int num_normal_procs;
    int num_nodes;
    float power_ratio;
    int mod;
    Splitter* splitter;

    int* p_in;
    hash_t p_out;

    HybridMap(int num_nodes, int pid, int num_procs, Splitter* splitter): pid(pid), num_nodes(num_nodes), num_procs(num_procs), splitter(splitter){
        
        num_slow_procs = splitter->num_slow_procs;
        num_normal_procs = splitter->num_normal_procs;
        power_ratio = splitter->power_ratio;
        mod = splitter->mod;
        
        int extra_nodes;
        if (pid<num_normal_procs) {
            num_local_nodes = floor(num_nodes/mod)*power_ratio;
            extra_nodes = num_nodes%mod;
            if (extra_nodes >= (pid*power_ratio)) {
                if (extra_nodes >= (pid+1)*power_ratio) {
                    num_local_nodes += power_ratio;
                }
                else {
                    num_local_nodes += extra_nodes-(pid*power_ratio);
                }
            }
        }
        else {
            num_local_nodes = floor(num_nodes/mod);
            extra_nodes = num_nodes%mod;
            if (extra_nodes > (num_normal_procs*power_ratio + (pid-num_normal_procs)))
                num_local_nodes++;
        }
        
        p_in = new int[num_local_nodes];
        
        if (pid<num_normal_procs) {
            int i = 0, power_iter_count = 0, power_ratio_iter=0;
            while (i<num_local_nodes) {
                p_in[i] = mod*power_iter_count + pid*power_ratio + power_ratio_iter;
                i++;
                power_ratio_iter = (power_ratio_iter+1)%int(power_ratio);
                if (power_ratio_iter==0) {
                    power_iter_count++;
                }
            }
        }
        else {
            for(int i=0; i<num_local_nodes; i++){
                p_in[i] = i*mod + num_normal_procs*power_ratio + (pid-num_normal_procs);
            }
        }
        // memset(p_in, -1, sizeof(int) * num_local_nodes);
    }

    int& operator[] (int u){  
        if(splitter->get_pid_for_node(u) == pid) {
            int local_u;
            if (pid<num_normal_procs) {
                local_u = floor(u/mod)*power_ratio + (u%mod)-(pid*power_ratio);
            }
            else {
                local_u = floor(u/mod);
            }
            return p_in[local_u];
        }
        else{
            if(p_out.find(u) == p_out.end()){
                p_out[u] = u;
            }
            return p_out[u];
        }
    }

};