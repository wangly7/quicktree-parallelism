#include <vector>
#include <tuple>

struct GpuTree {
    // CPU memory
    std::vector<double>& data; 
    uint32_t h_size;
    uint32_t *h_nodes;

    // GPU memoory
    double *mat;  // input
    double *u_i; // input
    uint32_t d_size;  // input
    uint32_t *d_nodes;  // output
    uint32_t *d_leftovers; // output


    void allocateMem();
    void transferMat2Device();
    std::vector<double> buildTree (std::vector<double>);
    std::vector<double> transferTree2Host(std::vector<double>& data, uint32_t size);
    void clearAndReset();
};