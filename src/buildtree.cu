#include "buildtree.cuh"
#include <cfloat>
#include <cstdint>
#include <limits>
#include <tuple>
#define BLOCKSIZE 256
#define NUMBLOCKS 8


__global__ void initActive(uint32_t N, uint32_t* active) {
    uint32_t bx = blockIdx.x;
    uint32_t tx = threadIdx.x;
    uint32_t bs = blockDim.x;
    uint32_t gs = gridDim.x;
    for (uint32_t i=bs*bx+tx; i < N; i+=bs*gs) {
        active[i] = 1;
    }
}

/*
*  Compute u_i array values
*/
__global__ void computeU(
    uint32_t N,
    const double* mat,
    const uint32_t* active,
    uint32_t active_count,
    double* u_i
) {
    // Kernel Configuration
    uint32_t bx = blockIdx.x;
    uint32_t tx = threadIdx.x;
    uint32_t bs = blockDim.x;
    uint32_t gs = gridDim.x;

    for (uint32_t i = bs*bx+tx; i < N; i+=bs*gs) {
        if (active[i] == 0) continue;

        double sum_ui = 0.0;
        for (uint32_t j = 0; j < N; j++){
            if (j==i) continue;
            if (active[j] == 0) continue;

            uint32_t a = (i > j) ? i : j;
            uint32_t b = (i > j) ? j : i;
            uint32_t index = a * (a - 1) / 2 + b;
            double d_ij = mat[index];
            sum_ui += d_ij;
        }
        u_i[i] = sum_ui / double(active_count - 2);
    }
}

/*
*  Find minimum value in a block
*/
__global__ void findBlockMin(
    uint32_t N,
    const double* mat,
    const double* u_i,
    const uint32_t* active,
    uint32_t* d_block_min_ij,
    double* d_block_min_dst
) {
    uint32_t bx = blockIdx.x;
    uint32_t tx = threadIdx.x;
    uint32_t bs = blockDim.x;
    uint32_t gs = gridDim.x;

    __shared__ RdcDst reduction[BLOCKSIZE];
    
    reduction[tx].index_i = UINT32_MAX;
    reduction[tx].index_j = UINT32_MAX;
    reduction[tx].dst = DBL_MAX;
    __syncthreads();

    for (uint32_t i = bs * bx + tx; i < N; i += bs * gs) {
        if (active[i] == 0) continue;
        for (uint32_t j = 0; j < i; j++) {
            if (active[j] == 0) continue;
            uint32_t index = i * (i - 1) / 2 + j;
            double q = mat[index] - u_i[i] - u_i[j];
            if (q < reduction[tx].dst) { 
                reduction[tx].dst = q;
                reduction[tx].index_i = i;
                reduction[tx].index_j = j;
            }
        }
    }
    __syncthreads();

    for (uint32_t s = bs / 2; s>0; s>>=1) { 
        if (tx < s) {
            if(reduction[tx].dst > reduction[tx+s].dst) reduction[tx] = reduction[tx+s];
        }
        __syncthreads();
    }

    if (tx == 0) {
        d_block_min_ij[2*bx] = reduction[tx].index_i;
        d_block_min_ij[2*bx+1] = reduction[tx].index_j;
        d_block_min_dst[bx] = reduction[tx].dst;
    }
}

/**
 * Use 1 block to find global minimum
 */
__global__ void findGlobalMin (
    const uint32_t* d_block_min_ij,
    const double* d_block_min_dst,
    uint32_t* d_min_ij
) {
    uint32_t tx = threadIdx.x;
    uint32_t bs = blockDim.x;

    __shared__ RdcDst reduction[BLOCKSIZE];

    RdcDst best;
    best.index_i = UINT32_MAX;
    best.index_j = UINT32_MAX;
    best.dst = DBL_MAX;

    // in case of threads num less than blocks num
    for(uint32_t k = tx; k < NUMBLOCKS; k += bs){
        RdcDst current;
        current.index_i = d_block_min_ij[2*k];
        current.index_j = d_block_min_ij[2*k+1];
        current.dst = d_block_min_dst[k];
        
        if (current.dst < best.dst) best = current;
    }

    reduction[tx] = best;
    __syncthreads();

    for (uint32_t s = bs / 2; s > 0; s >>=1){
        if (tx < s) {
            if (reduction[tx+s].dst < reduction[tx].dst) {
                reduction[tx] = reduction[tx+s];
            }
        }
        __syncthreads();
    }
    
    if (tx == 0) {
        d_min_ij[0] = reduction[0].index_i;
        d_min_ij[1] = reduction[0].index_j;
    }
}

/*
*  Use 1 thread to update merge info
*/
__global__ void updateMergeInfo(
    double* mat,
    double* u_i,
    uint32_t* active,
    const uint32_t* d_min_ij,
    double* d_merged_dst
) {
    uint32_t min_i = d_min_ij[0];
    uint32_t min_j = d_min_ij[1];

    uint32_t a = (min_i > min_j) ? min_i : min_j;
    uint32_t b = (min_i > min_j) ? min_j : min_i;
    uint32_t ij_index = a * (a - 1) / 2 + b;

    double dij = mat[ij_index];

    double di = 0.5 * (dij + (u_i[min_i] - u_i[min_j]));
    double dj = dij - di;

    // simple negative branch correction
    if (di < 0.0) {
        di = 0.0;
        dj = (dij < 0.0) ? 0.0 : dij;
    } else if (dj < 0.0) {
        dj = 0.0;
        di = (dij < 0.0) ? 0.0 : dij;
    }

    d_merged_dst[0] = di;
    d_merged_dst[1] = dj;

    // keep min_i, drop min_j
    active[min_j] = 0;
}

__global__ void updateDistanceMatrix(
    uint32_t N,
    double* mat,
    const uint32_t* active,
    const uint32_t* d_min_ij
){
    uint32_t bx = blockIdx.x;
    uint32_t tx = threadIdx.x;
    uint32_t bs = blockDim.x;
    uint32_t gs = gridDim.x;

    uint32_t min_i = d_min_ij[0];
    uint32_t min_j = d_min_ij[1];

    uint32_t a = (min_i > min_j) ? min_i : min_j;
    uint32_t b = (min_i > min_j) ? min_j : min_i;
    uint32_t ij_index = a * (a - 1) / 2 + b;
    double dij = mat[ij_index];

    for (uint32_t k = bs*bx+tx; k < N; k += bs*gs) {
        if (active[k] == 0) continue;
        if (k == min_i || k == min_j) continue;
        
        uint32_t ik_a = (min_i > k) ? min_i : k;
        uint32_t ik_b = (min_i > k) ? k : min_i;
        uint32_t ik_index = ik_a * (ik_a - 1) / 2 + ik_b;

        uint32_t jk_a = (min_j > k) ? min_j : k;
        uint32_t jk_b = (min_j > k) ? k : min_j;
        uint32_t jk_index = jk_a * (jk_a - 1) / 2 + jk_b;

        double dik = mat[ik_index];
        double djk = mat[jk_index];

        double new_d = 0.5 * (dik + djk - dij);
        if (new_d < 0.0) new_d = 0.0;
        mat[ik_index] = new_d;
    }
}

/**
 * Main orchestration method.
 * 1. Allocates GPU memory
 * 2. Transfers data from CPU
 * 3. Launches Kernel
 * 4. Retrieves results and reconstructs trees
 */
 void GpuTree::build() {
    int numBlocks = NUMBLOCKS;  // i.e. number of thread blocks on the GPU
    int blockSize = BLOCKSIZE; // i.e. number of GPU threads per thread block

    // 1. Allocate memory on Device
    allocateMem();

    // 2. Transfer sequence to device
    transferMat2Device();

    // 3. initialize active nodes array
    initActive<<<numBlocks, blockSize>>>(N, active);
    cudaDeviceSynchronize();

    // 4. initialize nodes array on CPU
    initNodesOnCPU();

    // 5. Main NJ loop
    uint32_t active_count = N;

    while (active_count > 3) {
        computeU<<<numBlocks, blockSize>>>(N, mat, active, active_count, u_i);

        findBlockMin<<<numBlocks, blockSize>>>(
            N, mat, u_i, active, d_block_min_ij, d_block_min_dst
        );

        findGlobalMin<<<1, blockSize>>>(
            d_block_min_ij, d_block_min_dst, d_min_ij
        );

        updateMergeInfo<<<1, 1>>>(
            mat, u_i, active, d_min_ij, d_merged_dst
        );

        updateDistanceMatrix<<<numBlocks, blockSize>>>(
            N, mat, active, d_min_ij
        );

        cudaDeviceSynchronize();
        MergeInfo info = transferNode2Host();
        buildInternalNode(info);
        
        active_count--;
    }

    // instantiate tree
    theTree = new Tree();
    handleLeftovers(); 
    
 }

void GpuTree::initNodesOnCPU () {
    nodes.resize(N);

    // for(i = 0; i < numseqs; i++) {
    for(uint32_t i = 0; i < N; i++) {
        nodes[i] = new TNode();
        nodes[i]->nodenumber = i;
        nodes[i]->identifier = identifiers[i];
    }
}


void GpuTree::buildInternalNode(const MergeInfo& info) {
    uint32_t i = info.min_ij[0];
    uint32_t j = info.min_ij[1];

    TNode* parent = new TNode();
    parent->left = nodes[i];
    parent->right = nodes[j];
    parent->identifier = "";
    parent->distance = 0.0;

    nodes[i]->distance = info.branch_len[0];
    nodes[j]->distance = info.branch_len[1];

    nodes[i] = parent;
    nodes[j] = nullptr;
}

void GpuTree::handleLeftovers() {
    std::vector<uint32_t> leftovers(3);

    for (uint32_t k = 0, m = 0; k < N; k++) {
        if (nodes[k] != NULL){
            theTree->child[m] = nodes[k];
            leftovers[m++] = k;
            nodes[k] = NULL;  // going to merge these nodes
        };
    }

    uint32_t index_1_0 = leftovers[1] * (leftovers[1] - 1) + leftovers[0];
    uint32_t index_2_0 = leftovers[2] * (leftovers[2] - 1) + leftovers[0];
    uint32_t index_2_1 = leftovers[2] * (leftovers[2] - 1) + leftovers[1];
    double dist_i = theTree->child[0]->distance = (
        data[index_1_0]+ data[index_2_0] - data[index_2_1]
    ) * 0.5;
    double dist_j = data[index_1_0] - theTree->child[0]->distance;
    double dist_k = data[index_2_0] - theTree->child[0]->distance;

    if(dist_i < 0.0) {
        dist_i = 0.0;
        dist_j = data[index_1_0];
        dist_k = data[index_2_0];
        if (dist_j < 0.0) {
            dist_j = 0.0;
            dist_k = (
                data[index_2_0] + data[index_2_1]
            ) * 0.5;
            if (dist_k < 0.0) {
                dist_k = 0.0;
            }
        }else if (dist_k < 0.0) {
            dist_k = 0.0;
            dist_j = (
                data[index_1_0] + data[index_2_1]
            ) * 0.5;
            if (dist_j < 0.0) {
                dist_j = 0.0;
            }
        }
    }else if (dist_j < 0.0) {
        dist_j = 0.0;
        dist_i = data[index_1_0];
        dist_k = data[index_2_1];
        if (dist_i < 0.0) {
            dist_i = 0.0;
            dist_k = (
                data[index_2_0] + data[index_2_1]
            ) * 0.5;
            if (dist_k < 0.0) dist_k = 0.0;
        }else if (dist_k < 0.0) {
            dist_k = 0.0;
            dist_i = (
                data[index_1_0] + data[index_2_0]
            ) * 0.5;
            if (dist_i < 0.0) dist_i = 0.0;
        }
    }else if (dist_k < 0.0) {
        dist_k = 0.0;
        dist_i = data[index_2_0];
        dist_j = data[index_2_1];
        if (dist_i < 0.0) {
            dist_i = 0.0;
            dist_j = (
                data[index_1_0] + data[index_2_1]
            ) * 0.5;
            if (dist_j < 0.0) dist_j = 0.0;
        }else if (dist_j < 0.0) {
            dist_j = 0.0;
            dist_i = (
                data[index_1_0] + data[index_2_0]
            ) * 0.5;
            if (dist_i < 0.0) dist_i = 0.0;
        }
    }

    theTree->child[0]->distance = dist_i;
    theTree->child[1]->distance = dist_j;
    theTree->child[2]->distance = dist_k;

}

/**
 * Allocates GPU memory for distance matrix 1D array, nodes status, and merged nodes.
 */
void GpuTree::allocateMem() {
    cudaError_t err;

    // 1. allocate memory to distance matrix
    err = cudaMalloc(&mat, (N * (N - 1) / 2) * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 2. allocate memory to nodes status array 
    err = cudaMalloc(&active, N * sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
    
    // 3. allocate memory for u_i
    err = cudaMalloc(&u_i, N * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 4. allocate memory for min_i, min_j
    err = cudaMalloc(&d_min_ij, 2 * sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 5. allocate memory for distance of min i and min j from merged ij
    err = cudaMalloc(&d_merged_dst, 2 * sizeof(double));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // 6. allocate memeory for min i and min j from different blocks
    err = cudaMalloc(&d_block_min_ij, 2 * NUMBLOCKS * sizeof(uint32_t));
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    err = cudaMalloc(&d_block_min_dst, NUMBLOCKS * sizeof(double));
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

    err = cudaMemcpy(mat, data.data(), (N*(N-1)/2) * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
}

/**
 * Copies the computed merged nodes information from GPU back to Host.
 */
MergeInfo GpuTree::transferNode2Host() {
    MergeInfo info;

    cudaError_t err = cudaMemcpy(info.min_ij, d_min_ij, 2 * sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }

    // transfer merged nodes distance from device to host
    err = cudaMemcpy(info.branch_len, d_merged_dst, 2 * sizeof(double), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: %s (%s)\n", cudaGetErrorString(err), cudaGetErrorName(err));
        exit(1);
    }
    return info;
}

void GpuTree::clearAndReset() {
    cudaFree(mat);
    cudaFree(u_i);
    cudaFree(active);
    cudaFree(d_min_ij);
    cudaFree(d_merged_dst);
    cudaFree(d_block_min_ij);
    cudaFree(d_block_min_dst);
}