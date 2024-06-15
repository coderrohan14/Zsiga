#include <math.h>

class Splitter {
    
    public:
    
    float power_ratio;
    int num_procs;
    int num_slow_procs;
    int num_normal_procs;
    int mod;

    Splitter(float power_ratio, int num_procs, int num_slow_procs): power_ratio(power_ratio), num_procs(num_procs), num_slow_procs(num_slow_procs) {
        num_normal_procs = num_procs - num_slow_procs;
        mod = power_ratio*num_normal_procs + num_slow_procs;
    }

    int get_pid_for_node(int node) {
        int unode = node%mod;
        if (unode < num_normal_procs*power_ratio)
        {
            unode = floor(unode/power_ratio);
        }
        else
        {
            unode = num_normal_procs + (unode - num_normal_procs*power_ratio);
        }
        return unode;
    }
};