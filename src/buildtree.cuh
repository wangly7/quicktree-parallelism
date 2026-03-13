#include <string>
#include <vector>
#include "tree.hpp"

void printGpuProperties();

struct RdcDst {
    uint32_t index_i;
    uint32_t index_j;
    double dst;
};


struct MergeInfo{
    uint32_t min_ij[2];
    double branch_len[2];
};

struct GpuTree {
    // CPU memory
    std::vector<double>& data;
    std::vector<TNode*> nodes;
    std::vector<std::string>& identifiers;
    uint32_t N;
    Tree* theTree;

    // GPU memory
    double *mat;  // input
    double *u_i; // input

    uint32_t *d_block_min_ij; // min ij from different blocks
    double* d_block_min_dst; // min distance from different blocks

    uint32_t *d_min_ij; // min i and min j packed into one
    double *d_merged_dst; // [0] distance from i to ij
                         // [1] distance from j to ij

    double* partial_sum;

    uint32_t *active;  // active nodes that need to merge

     GpuTree(std::vector<double>& d, std::vector<std::string>& ids, uint32_t n): 
        data(d), identifiers(ids), N(n), 
        theTree(nullptr), 
        mat(nullptr), 
        u_i(nullptr), 
        d_block_min_ij(nullptr),
        d_block_min_dst(nullptr),
        d_min_ij(nullptr),
        d_merged_dst(nullptr),
        active(nullptr) {}

    void allocateMem();
    void transferMat2Device();
    void build();
    void initNodesOnCPU();
    MergeInfo transferNode2Host();
    void transferMat2Host();
    void buildInternalNode(const MergeInfo& info);
    void handleLeftovers();
    void clearAndReset();
};