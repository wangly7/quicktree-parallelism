// NJ Algorithm
#include "matcal.hpp"
#include <cfloat>
#include <cstdint>
#include <limits>
#include <tuple>

std::tuple<uint32_t, uint32_t, double, double> compute_min_ij (const DistanceMatrix* mat) {
    uint32_t N = mat->size;
    double u_i[N];

    double min_leaf = DBL_MAX;
    uint32_t min_i = std::numeric_limits<uint32_t>::min();
    uint32_t min_j = std::numeric_limits<uint32_t>::min();

    if (N < 3) return std::make_tuple(min_i, min_j, DBL_MIN, DBL_MIN);

    // 1. for each leaf compute u_i = summation from j!= i to N, (d_ij) / (N-2)
    for (uint32_t i=0; i<N; i++) {
        double sum_ui = 0.0;
        for (uint32_t j=0; j<i; j++){
            double d_ij = mat->get(i, j);
            sum_ui += d_ij;
        }
        u_i[i] = sum_ui / (N-2);
    }

    // 2. find leaf pair i and j for which d_ij - u_i - u_j is minimum
    for (uint32_t i=1; i<N; i++) {
        for (uint32_t j=0; j<i; j++) {
            double ui = u_i[i];
            double uj = u_i[j];
            double d_ij = mat->get(i, j); 
            double leaf = d_ij - ui - uj;
            
            if (leaf < min_leaf) {
                min_leaf = leaf;
                min_i = i;
                min_j = j;
            }
        }
    }

    // 3. Join i and j to form a common node (ij) such that
    // the distance of i from the (ij), d_i(ij) = 0.5 * (d_ij + (u_i-u_j)) and 
    // the distance of j from the (ij), d_j(ij) = d_ij - d_i(ij)
    double d_ij = mat->get(min_i, min_j); 
    double ui = u_i[min_i];
    double uj = u_i[min_j];
    double d_i_ij =  0.5 * (d_ij + (ui - uj));
    double d_j_ij = d_ij - d_i_ij;

    // adjustment for negative branch length
    if (d_i_ij < 0.0) {
        d_i_ij = 0.0;
        d_j_ij = (d_ij < 0.0)? 0.0 : d_ij;
    }
    else if (d_j_ij < 0.0) {
        d_j_ij = 0.0;
        d_i_ij = (d_ij < 0.0)? 0.0 : d_ij;
    }

    return std::make_tuple(min_i, min_j, d_i_ij, d_j_ij);
}

DistanceMatrix update_matrix (const DistanceMatrix* mat) {

    
    uint32_t min_i, min_j;
    double d_i_ij, d_j_ij, new_d;
    
    std::tie(min_i, min_j, d_i_ij, d_j_ij) = compute_min_ij(mat);
    double d_ij = d_i_ij + d_j_ij;

    uint32_t N = mat->size;
    DistanceMatrix new_mat(N-1);

    // 4. Treat (ij) as a new tip, ignoring previous tips i and j
    // 5. Distance of (ij) from other tips k, d_k(ij) = (d_ik+d_jk-d_ij) * 0.5 

    // use minimum index for merged, and disable the other one
    uint32_t a = std::min(min_i, min_j);
    uint32_t b = std::max(min_i, min_j);

    // helper function
    auto map_idx = [b](uint32_t old_idx) -> uint32_t {
        return (old_idx < b) ? old_idx : old_idx - 1;
    };

    for (uint32_t i=0; i<N; i++) {
        if (i==b) continue;
        for (uint32_t j=0; j<i; j++) {
            if (j==b) continue;
            // shift indices above b down by 1
            uint32_t ni = map_idx(i);
            uint32_t nj = map_idx(j);

            // calculate new distance
            if (i!=a && j!=a) {
                new_d = mat->get(i, j);
            }
            else {
                double d_ik = (i==a)? mat->get(a, j): mat->get(i, b);
                double d_jk = (i==a)? mat->get(i, a): mat->get(b, j);
                new_d = (d_ik + d_jk - d_ij) * 0.5;
            }

            // adjust for negative branch length => check if we need this
            new_d = (new_d < 0.0)? 0.0: new_d;

            new_mat.set(ni, nj, new_d); 
        }
    }
    return new_mat;
}