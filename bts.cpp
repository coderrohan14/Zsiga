#include <mpi.h>
#include <cstdlib>
#include <stdlib.h>
#include <thread>
#include <iostream>
#include <chrono>
#include <sys/resource.h>
#include "rem_bts.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <fstream>
#include <cstring>
#include <string>
#include <array>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <map>
#include <vector>
#include <unistd.h>
#include <regex>
#include <iomanip>

#include "MessageManager.h"
#include "ska/flat_hash_map.hpp"
using namespace std;
using namespace std::chrono;
using namespace bts;

#define read_buffer_size 1024
#define recv_buffer_size 1024

Edge* read_buffer = new Edge[read_buffer_size];
Edge* recv_buffer = new Edge[recv_buffer_size];

typedef ska::flat_hash_map<int, int> hash_t;

MPI_Datatype dt_edge;
extern char **environ;


void receiver(int num_procs, MessageManager* msgmng, RemBTS* rembts, size_t* total_received_size){
    int remain = num_procs;

    MPI_Status status;
    int received_size;
    *total_received_size = 0;

    while(remain){
        MPI_Recv(recv_buffer, msgmng->capacity, dt_edge, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, dt_edge, &received_size);

        for(int i=0; i<received_size; i++){
            rembts->merge(recv_buffer[i].u, recv_buffer[i].v);
        }

        *total_received_size += received_size;
        
        if(status.MPI_TAG == EndMessage) remain--;
    }
}

void simply_emitting(int pid, int num_edges, int num_procs, MessageManager* msgmng, FILE* f){
    size_t num_read_edges = 0;
    for(int i=0; i<num_edges;){

        num_read_edges = num_edges - i < read_buffer_size ? num_edges - i : read_buffer_size;
        num_read_edges = fread(read_buffer, sizeof(Edge), num_read_edges, f);
        
        for(int j=0; j<num_read_edges; j++){
            Edge e = read_buffer[j];
            if(e.u < e.v){
                e = Edge{e.v, e.u};
            }
            int up = e.u % num_procs;
            int vp = e.v % num_procs;
            msgmng->emit(&e, up);
            if(up != vp){
                msgmng->emit(&e, vp);
            }
        }
        i += num_read_edges;
    }

    for (int step = 0; step < num_procs; step++) {
        int id_dst = (pid - step + num_procs) % num_procs;
        msgmng->flush(id_dst);
        msgmng->endAck(id_dst);
    }
}


// fix seg fault here
void initialization_step(int pid, int num_nodes, int num_edges, int num_procs, Edge* emit_buffer, int emit_buffer_size, MessageManager* msgmng, FILE* f, HeterogeneousSplitter* splitter){

    // std::cout << "Process " << pid << " Inside initialization  step." << std::endl;

    RemInit* reminit = new RemInit(num_nodes, num_procs, splitter);


    size_t num_read_edges = 0;
    for(int i=0; i<num_edges;){
        while(reminit->p.size() + read_buffer_size <= emit_buffer_size && i < num_edges){
            num_read_edges = num_edges - i < read_buffer_size ? num_edges - i : read_buffer_size;            
            num_read_edges = fread(read_buffer, sizeof(Edge), num_read_edges, f);
            for(int j=0; j<num_read_edges; j++){
                reminit->merge(read_buffer[j].u, read_buffer[j].v);
            }
            i += num_read_edges;
        }

        int size_to_emit = reminit->partition(emit_buffer);

        for(int j=0; j<size_to_emit; j++){
            int up = splitter->get_pid_for_node(emit_buffer[j].u);
            int vp = splitter->get_pid_for_node(emit_buffer[j].v);

            msgmng->emit(&emit_buffer[j], up);
            if(up != vp){
                msgmng->emit(&emit_buffer[j], vp);
            }
        }

        reminit->clear();
    }

    for (int step = 0; step < num_procs; step++) {
        int id_dst = (pid - step + num_procs) % num_procs;
        msgmng->flush(id_dst);
        msgmng->endAck(id_dst);
    }
    
    delete reminit;
}


FILE* open_and_seek(char* input, int pid, int num_procs, long& num_edges, HeterogeneousSplitter* splitter) {
    FILE* f = fopen(input, "rb");
    fseek(f, 0, SEEK_END);

    long num_total_edges = ftell(f) / sizeof(Edge);
    std::vector<long> edge_counts(num_procs);

    // Calculate the number of edges for each process based on their capabilities
    for (int i = 0; i < num_procs; ++i) {
        float ratio = splitter->node_capabilities[i].weighted_sum(splitter->cpu_weight, splitter->memory_weight, splitter->bandwidth_weight) / splitter->total_power;
        edge_counts[i] = std::ceil(num_total_edges * ratio);
    }

    // Adjust edge counts to ensure the total matches num_total_edges
    long total_assigned = std::accumulate(edge_counts.begin(), edge_counts.end(), 0L);
    int i = 0;
    while (total_assigned > num_total_edges) {
        if (edge_counts[i] > 0) {
            edge_counts[i]--;
            total_assigned--;
        }
        i = (i + 1) % num_procs;
    }

    // Calculate the starting position for this process
    long edge_start = 0;
    for (int i = 0; i < pid; ++i) {
        edge_start += edge_counts[i];
    }

    long pos_start = edge_start * sizeof(Edge);
    num_edges = edge_counts[pid];

    fseek(f, pos_start, SEEK_SET);

    fprintf(stderr, "(%d/%d) num_total_edges: %ld, pos_start: %ld, num_edges: %ld\n", pid, num_procs, num_total_edges, pos_start, num_edges);
    return f;
}

// Function to read a value from a file
std::string readFile(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file) {
        std::cerr << "Error: Cannot open file " << filePath << std::endl;
        return "";
    }
    std::string value;
    std::getline(file, value);
    return value;
}

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

std::string getBandwidthLimit() {
    std::string tc_output = exec("sudo tc class show dev $(ip -o -4 route show to default | awk '{print $5}')");
    std::regex rate_regex("rate ([0-9.]+[KMG]?bit)");
    std::smatch match;
    if (std::regex_search(tc_output, match, rate_regex)) {
        return match[1];
    }
    return "No limit set or unable to retrieve";
}

// Function to get CPU frequency
double getCPUFrequencyMHz() {
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line)) {
        if (line.substr(0, 7) == "cpu MHz") {
            std::istringstream iss(line.substr(line.find(":") + 1));
            double value;
            if (iss >> value) {
                return value;
            }
        }
    }
    return 0.0;  // Return 0 if unable to find CPU frequency
}

int main(int argc, char** argv){

    auto start = high_resolution_clock::now();

    int provided, num_procs, pid;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Type_contiguous(2, MPI_INT, &dt_edge);
    MPI_Type_commit(&dt_edge);

     // Paths to cgroup files
    std::string cpuMaxFile = "/sys/fs/cgroup/mpi_limits/cpu.max";
    std::string memMaxFile = "/sys/fs/cgroup/mpi_limits/memory.max";

    // Read CPU limit
    double limitedCPUFrequency = 0.0;
    std::string cpuMax = readFile(cpuMaxFile);
    if (!cpuMax.empty()) {
        std::istringstream cpuStream(cpuMax);
        long long cpuLimit, cpuPeriod;
        cpuStream >> cpuLimit >> cpuPeriod;
        double cpuPercentage = (double)cpuLimit / cpuPeriod;
        
        double totalCPUFrequency = getCPUFrequencyMHz();
        limitedCPUFrequency = totalCPUFrequency * cpuPercentage;
        
        std::cout << "Process " << pid << " - CPU Limit: " << std::fixed << std::setprecision(2) 
                  << limitedCPUFrequency << " MHz" << std::endl;
    }

    // Read Memory limit
    float memLimitMB = 0.0f;
    std::string memMax = readFile(memMaxFile);
    if (!memMax.empty()) {
        long long memLimitBytes = std::stoll(memMax);
        memLimitMB = static_cast<float>(memLimitBytes) / (1024 * 1024); // Convert to MB
        std::cout << "Process " << pid << " - Memory Limit (in MB): " << memLimitMB << std::endl;
    }

    // Get Bandwidth limit
    float bwLimitMBps = 0.0f;
    std::string bwLimit = getBandwidthLimit();
    if (bwLimit != "No limit set or unable to retrieve") {
        std::string numPart = bwLimit.substr(0, bwLimit.find("Mbit"));
        bwLimitMBps = std::stof(numPart) / 8.0f; // Convert Mbit/s to MByte/s
        std::cout << "Process " << pid << " - Bandwidth Limit (in MBps): " << bwLimitMBps << std::endl;
    }

    // Collect resource information from all nodes
    NodeCapabilities localCapabilities = {
        static_cast<float>(limitedCPUFrequency),
        memLimitMB,
        bwLimitMBps
    };

    // Define a buffer size that's large enough to hold the NodeCapabilities data
const int BUFFER_SIZE = sizeof(float) * 3; // Adjust if needed

// Create send and receive buffers
char sendBuffer[BUFFER_SIZE];
std::vector<char> recvBuffer(BUFFER_SIZE * num_procs);

// Pack the data
int position = 0;
MPI_Pack(&localCapabilities.cpu, 1, MPI_FLOAT, sendBuffer, BUFFER_SIZE, &position, MPI_COMM_WORLD);
MPI_Pack(&localCapabilities.memory, 1, MPI_FLOAT, sendBuffer, BUFFER_SIZE, &position, MPI_COMM_WORLD);
MPI_Pack(&localCapabilities.bandwidth, 1, MPI_FLOAT, sendBuffer, BUFFER_SIZE, &position, MPI_COMM_WORLD);

// Gather the data
MPI_Allgather(sendBuffer, BUFFER_SIZE, MPI_PACKED, recvBuffer.data(), BUFFER_SIZE, MPI_PACKED, MPI_COMM_WORLD);

// Unpack the data
std::vector<NodeCapabilities> allCapabilities(num_procs);
for (int i = 0; i < num_procs; ++i) {
    position = i * BUFFER_SIZE;
    MPI_Unpack(recvBuffer.data(), recvBuffer.size(), &position, &allCapabilities[i].cpu, 1, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(recvBuffer.data(), recvBuffer.size(), &position, &allCapabilities[i].memory, 1, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(recvBuffer.data(), recvBuffer.size(), &position, &allCapabilities[i].bandwidth, 1, MPI_FLOAT, MPI_COMM_WORLD);
}



    char* input = argv[1];
    int num_nodes = strtol(argv[2], NULL, 10);
    fprintf(stderr, "(%d/%d) num_nodes: %d\n", pid, num_procs, num_nodes);

    // Initialize HeterogeneousSplitter with the collected capabilities
    HeterogeneousSplitter* splitter = new HeterogeneousSplitter(allCapabilities, num_nodes);

    // if(pid == 0){
    //     std::cout << "Node Distribution: ";
    //     for(int i = 0; i < num_nodes; ++i) {
    //         std::cout << splitter->get_pid_for_node(i) << " ";
    //     }
    //     std::cout << std::endl;
    // }
    
    MessageManager* msgmng = new MessageManager(num_procs, dt_edge);
    RemBTS* rembts = new RemBTS(num_nodes, pid, splitter);
    size_t total_communication, total_received_size;

    // std::vector<int> node_distribution = splitter->getNodeDistribution();
    
    // Calculate the emit buffer size based on the distribution
    int emit_buffer_size = splitter->getNodeCountForProcessor(pid) * 10;  

    // Ensure the buffer size is at least as large as read_buffer_size
    if (emit_buffer_size < read_buffer_size) {
        emit_buffer_size = read_buffer_size;
    }

    // std::cout << "Emit Buffer Size: " << emit_buffer_size << std::endl;
    
    Edge* emit_buffer = new Edge[emit_buffer_size];

    long num_edges;

    // Before open_and_seek
// std::cout << "Process " << pid << " - About to call open_and_seek" << std::endl;

    FILE* f = open_and_seek(input, pid, num_procs, num_edges, splitter);

    // After open_and_seek
    // std::cout << "Process " << pid << " - open_and_seek completed" << std::endl;
    thread recv_thread(receiver, num_procs, msgmng, rembts, &total_received_size);
    initialization_step(pid, num_nodes, num_edges, num_procs, emit_buffer, emit_buffer_size, msgmng, f, splitter);
    // std::cout << "Process " << pid << "initialization  step completed" << std::endl;
    fclose(f);

    milliseconds elapsed_ms = duration_cast<milliseconds>(high_resolution_clock::now()-start);
    recv_thread.join();
    fprintf(stderr, "(Intialize, pid %d) recv %ld, %ld ms.\n", pid, total_received_size, elapsed_ms.count());

    MPI_Allreduce(&total_received_size, &total_communication, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    struct rlimit limit1;
    limit1.rlim_cur = RLIM_INFINITY;
    limit1.rlim_max = RLIM_INFINITY;
    setrlimit(RLIMIT_CORE, &limit1);

    int round = 1;
    size_t num_changes = 1;
    size_t num_total_changes = 1;
    while(num_total_changes){
        auto timer1 = high_resolution_clock::now();
        num_changes = 0;
        int toemit_size = rembts->getEdgesToEmit(emit_buffer, num_changes);
        int64_t t1 = duration_cast<milliseconds>(high_resolution_clock::now()-timer1).count();
        
        MPI_Barrier(MPI_COMM_WORLD);

        auto timer2 = high_resolution_clock::now();

        thread recv_thread(receiver, num_procs, msgmng, rembts, &total_received_size);

        for(int j=0; j<toemit_size; j++){
            int up = splitter->get_pid_for_node(emit_buffer[j].u);
            int vp = splitter->get_pid_for_node(emit_buffer[j].v);

            if(up != pid){
                msgmng->emit(&emit_buffer[j], up);
            }
            if(vp != pid && up != vp){
                msgmng->emit(&emit_buffer[j], vp);
            }
        }

        for (int step = 0; step < num_procs; step++) {
            int id_dst = (pid - step + num_procs) % num_procs;
            msgmng->flush(id_dst);
            msgmng->endAck(id_dst);
        }
        
        int64_t t2 = duration_cast<milliseconds>(high_resolution_clock::now()-timer2).count();

        recv_thread.join();
        
        MPI_Allreduce(&total_received_size, &total_communication, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&num_changes, &num_total_changes, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        
        int64_t t3 = duration_cast<milliseconds>(high_resolution_clock::now()-timer1).count();

        fprintf(stderr, "(Round %d, pid %d) recv %ld, changes %ld, refine: %ld ms, send: %ld ms, total: %ld ms.\n", round++, pid, total_received_size, num_total_changes, t1, t2, t3);
        
    }

    for(int i = pid; i < num_nodes; i += num_procs){
        rembts->find(i);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
}