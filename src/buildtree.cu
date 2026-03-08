#include "buildtree.cuh"
#include <cfloat>
#include <cstdint>
#include <limits>
#include <tuple>
#define BLOCKSIZE 256

/**
 * Building phylogenetic tree on GPU
 */
 __global__ void buildTreeOnCPU(uint32_t d_size, double *mat, uint32_t *d_nodes, double* u_i){
    // Kernel Configuration
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int bs = blockDim.x;
    int gs = gridDim.x;

    uint32_t u_len = d_size*(d_size-1)/2;
    
    __shared__ uint32_t min_i;
    __shared__ uint32_t min_j;
    __shared__ double min_dst;
    struct red_dst {
        uint32_t index_i;
        uint32_t index_j;
        double dst;
    };
    __shared__ rdc_dst reduction[BLOCKSIZE]; //shared memory limit per block

    // one grid for one distance matrix
    // 1. for each leaf compute u_i = summation from j!= i to N, (d_ij) / (N-2)
    for (uint32_t i = bs*bx+tx; i<N; i+=bs*gs) {
        double sum_ui = 0.0;
        for (uint32_t j=0; j<i; j++){
            uint32_t index = i*(i-1)/2+j;
            if (d_nodes[index] == 0) continue; // skip inactive nodes
            double d_ij = mat[index]; 
            sum_ui += d_ij;
        }
        u_i[i] = sum_ui / (N-2);
    }

    // initialize reduction
    for (uint32_t i=bs*bx+tx; i<BLOCKSIZE; i+=bs*gs) {
        reduction[i].index_i = -1;
        reduction[i].index_j = -1;
        reduction[i].dst = DBL_MAX;
    }

    // 2. find leaf pair i and j for which d_ij - u_i - u_j is minimum 
    for (uint32_t i=bs*bx+tx; i<N; i+=bs*gs){
        for (uint32_t j=0; j<i; j++) {
            uint32_t index = i*(i-1)/2+j;
            if (d_nodes[index] == 0) continue; // skip inactive nodes
            double leaf = mat[index] - u_i[i] - u_j[j];
            if (leaf < reduction[i].dst) { // Issue: Out of bounds bc N > BLOCKSIZE
                reduction[i].dst = leaf;
                reduction[i].index_i = i;
                reduction[i].index_j = j;
            }
        }
    }

    // parallel reduction
    for (int s=BLOCKSIZE/2; s>0; s>>=1) { 
        if (tx < s) {
            if(reduction[tx].dst > reduction[tx+s].dst) reduction[tx] = reduction[tx+s];
        }
        __syncthreads();
    }

    if (tx == 0) {
        min_dst = reduction[0].dst;
        min_i = reduction[0].index_i;
        min_j = reduction[0].index_j;
    }
    __syncthreads();

    // 3. Join i and j to form a common node (ij) such that
    // the distance of i from the (ij), d_i(ij) = 0.5 * (d_ij + (u_i-u_j)) and 
    // the distance of j from the (ij), d_j(ij) = d_ij - d_i(ij)
    double d_i_ij =  0.5 * (min_dst + (u_i[min_i] - u_i[min_j]));
    double d_j_ij = mid_dst - d_i_ij;

    // adjustment for negative branch length
    if (d_i_ij < 0.0) {
        d_i_ij = 0.0;
        d_j_ij = (min_dst < 0.0)? 0.0 : d_ij;
    }
    else if (d_j_ij < 0.0) {
        d_j_ij = 0.0;
        d_i_ij = (min_dst < 0.0)? 0.0 : d_ij;
    }

    // 4. Treat (ij) as a new tip, ignoring previous tips i and j
    // 5. Distance of (ij) from other tips k, d_k(ij) = (d_ik+d_jk-d_ij) * 0.5 
    uint32_t a = std::min(min_i, min_j);
    uint32_t b = std::max(min_i, min_j);

    // set new merged distance for (a, k) for all tips k 
    for (uint32_t k=bs*bx+tx; k<N; k+=gs*bs) {
        uint32_t ak_index = a*(a-1)/2+k;
        uint32_t bk_index = b*(b-1)/2+k;
        if (k != a && k != b && d_nodes[k] == 1){
             double new_d = (mat[ak_index] + mat[bk_index] - mat[min_dst]) * 0.5;
             mat[ak_index] = (new_d < 0.0)? 0.0: new_d; // adjust for negative branch length 
        }
        d_nodes[bk_index] = 0; // set tips b to inactive
    }

    // Store in d_nodes: index, branch length
    // Should have global min_i and min_j for merging

 }



/**
 * Main orchestration method.
 * 1. Allocates GPU memory
 * 2. Transfers data from CPU
 * 3. Launches Kernel
 * 4. Retrieves results and reconstructs trees
 */
 void GpuTree::build() {
    int numBlocks = 512;  // i.e. number of thread blocks on the GPU
    int blockSize = BLOCKSIZE; // i.e. number of GPU threads per thread block

    // 1. Allocate memory on Device
    allocateMem();

    // 2. Transfer sequence to device
    transferMat2Device();

    // 3. Perform building tree on GPU
    buildTreeOnCPU<<<numBlocks, blockSize>>>(d_size, mat, d_nodes, u_i);
    

    // 4. Transfer the final merged nodes from device
    transferTree2Host();

    // 5. Get the phylogenetic tree with merged nodes
    
 }

void GpuTree::getPhylogenneticTree () {
    
}

/**
 * Allocates GPU memory for distance matrix 1D array, nodes status, and merged nodes.
 */
void GpuTree::allocateMem() {
    cudaError_t err;

    // 1. allocate memory to distance matrix
    err = cudaMalloc(&mat, (h_size * (h_size - 1) / 2) * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 2. allocate memory to d_size
    err = cudaMalloc(&d_size, sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 3. allocate memory to the final merged nodes array 
    err = cudaMalloc(&d_nodes, h_size * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 4. allocate memory to the  
    err = cudaMalloc(&d_leftovers, 3 * sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
    
    // 5. allocate memory for u_i
    err = cudaMalloc(&u_i, h_size * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
}

/*
Transfer distance matrix data and size to device
*/
void GpuTree::transferMat2Device() {
    cudaError_t err;

    err = cudaMemcpy(mat, data, (h_size*(h_size-1)/2) * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    err = cudaMemcpy(d_size, h_size, sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
}

/**
 * Copies the computed merged nodes information from GPU back to Host.
 */
std::vector<uint32_t> GpuTree::transferTree2Host() {
    std::vector<uint32_t> nodes;

    cudaError_t err = cudaMemcpy(h_nodes, d_nodes, size, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
    return nodes;
}

void GpuTree::clearAndReset() {
    cudaFree(mat);
    cudaFree(u_i);
    cudaFree(d_nodes);
    cudaFree(d_leftovers);
}