#include "buildtree.hpp"


struct Tree* neighbour_joining_buildtree(DistanceMatrix* mat, const std::vector<std::string>& identifiers) {
    uint32_t numseqs, i, j;
    uint32_t k, m, nodecount;
    uint32_t row, column;
    uint32_t nextfreenode;
    uint32_t mini = 0, minj = 0;
    double fnumseqs, ri, minsofar, dist, dij, dist_i, dist_j, dist_k, dmj, dmi;
    TNode* newnode;
    std::vector<uint32_t> leftovers(3);

    numseqs = mat->size;
    nextfreenode = numseqs;

    fnumseqs = (double) numseqs;

    std::vector<TNode*> nodes(numseqs);

    for(i = 0; i < numseqs; i++) {
        nodes[i] = new TNode();
        nodes[i]->nodenumber = i;
        nodes[i]->identifier = identifiers[i];
    }

    Tree* theTree = new Tree();

    if (numseqs > 2) {
        std::vector<double> r(numseqs);
        
        // calculate ri to update distance
        for (i = 0; i < numseqs; i++) {
            ri = 0.0;
            for (k = 0; k < numseqs; k++) {
                if (i == k) continue;
                if (k > i) ri += mat->get(k, i);
                else ri += mat->get(i, k);
            }
            r[i] = ri / (fnumseqs - 2.0);
        }

        // // main loop: need to merge numseq - 3 times
        for (nodecount = 0; nodecount < numseqs - 3; nodecount++) {
            minsofar = __FLT_MAX__;

            for(i = 0; i < numseqs; i++) {
                if (nodes[i] == NULL) continue;  // node i already merged
                for (j = 0; j < i; j++) {
                    if (nodes[j] == NULL) continue; // node j already merged
                    dist = mat->get(i, j) - (r[i] + r[j]);
                    if(dist < minsofar) {
                        mini = i;
                        minj = j;
                        minsofar = dist;
                    }
                }
            }

            dij = mat->get(mini, minj);
            dist_i = (dij + r[mini] - r[minj]) * 0.5;
            dist_j = dij - dist_i;

            // adjustment: allow for negative branch length
            if (dist_i < 0.0) {
                dist_i = 0.0;
                dist_j = dij;
                if (dist_j < 0) dist_j = 0.0;
            }else if (dist_j < 0.0){
                dist_j = 0.0;
                dist_i = dij;
                if (dist_i < 0.0) dist_i = 0.0;
            }

            // merge two nodes with minimum distance
            nodes[mini]->distance = dist_i;
            nodes[minj]->distance = dist_j;
            
            newnode = new TNode();
            nextfreenode++;
            newnode->left = nodes[mini];
            newnode->right = nodes[minj];
            nodes[mini]->parent = newnode;
            nodes[minj]->parent = newnode;
            nodes[mini] = newnode;
            nodes[minj] = NULL;

            // update distance matrix after merging
            r[mini] = 0.0;
            for (m = 0; m < numseqs; m++) {
                if (nodes[m] == NULL) continue;
                
                if (m != mini) {
                    if (m > minj) dmj = mat->get(m, minj);
                    else dmj = mat->get(minj, m);

                    if (m > mini) {
                        row = m;
                        column = mini;
                    }else{
                        row = mini;
                        column = m;
                    }
                    dmi = mat->get(row, column);

                    mat->set(row, column, (dmi + dmj - dij) * 0.5);
                    r[m] = ((r[m] * (fnumseqs - 2.0)) - dmi - dmj + mat->get(row, column)) / (fnumseqs - 3.0);
                    r[mini] += mat->get(row, column);
                }
            }

            fnumseqs -= 1.0;
            r[mini] /= fnumseqs - 2.0;
        }
        // end of main loop
        
        // handle 3 nodes left. 
        for (k = 0, m = 0; k < numseqs; k++) {
            if (nodes[k] != NULL){
                theTree->child[m] = nodes[k];
                leftovers[m++] = k;
                nodes[k] = NULL;  // going to merge these nodes
            };
        }

        dist_i = theTree->child[0]->distance = (
            mat->get(leftovers[1], leftovers[0]) + 
            mat->get(leftovers[2], leftovers[0]) - 
            mat->get(leftovers[2], leftovers[1])
        ) * 0.5;
        dist_j = mat->get(leftovers[1], leftovers[0]) - theTree->child[0]->distance;
        dist_k = mat->get(leftovers[2], leftovers[0]) - theTree->child[0]->distance;

        if(dist_i < 0.0) {
            dist_i = 0.0;
            dist_j = mat->get(leftovers[1], leftovers[0]);
            dist_k = mat->get(leftovers[2], leftovers[0]);
            if (dist_j < 0.0) {
                dist_j = 0.0;
                dist_k = (
                    mat->get(leftovers[2], leftovers[0]) + mat->get(leftovers[2], leftovers[1])
                ) * 0.5;
                if (dist_k < 0.0) {
                    dist_k = 0.0;
                }
            }else if (dist_k < 0.0) {
                dist_k = 0.0;
                dist_j = (
                    mat->get(leftovers[1], leftovers[0]) + mat->get(leftovers[2], leftovers[1])
                ) * 0.5;
                if (dist_j < 0.0) {
                    dist_j = 0.0;
                }
            }
        }else if (dist_j < 0.0) {
            dist_j = 0.0;
            dist_i = mat->get(leftovers[1], leftovers[0]);
            dist_k = mat->get(leftovers[2], leftovers[1]);
            if (dist_i < 0.0) {
                dist_i = 0.0;
                dist_k = (
                    mat->get(leftovers[2], leftovers[0]) + mat->get(leftovers[2], leftovers[1])
                ) * 0.5;
                if (dist_k < 0.0) dist_k = 0.0;
            }else if (dist_k < 0.0) {
                dist_k = 0.0;
                dist_i = (
                    mat->get(leftovers[1], leftovers[0]) + mat->get(leftovers[2], leftovers[0])
                ) * 0.5;
                if (dist_i < 0.0) dist_i = 0.0;
            }
        }else if (dist_k < 0.0) {
            dist_k = 0.0;
            dist_i = mat->get(leftovers[2], leftovers[0]);
            dist_j = mat->get(leftovers[2], leftovers[1]);
            if (dist_i < 0.0) {
                dist_i = 0.0;
                dist_j = (
                    mat->get(leftovers[1], leftovers[0]) + mat->get(leftovers[2], leftovers[1])
                ) * 0.5;
                if (dist_j < 0.0) dist_j = 0.0;
            }else if (dist_j < 0.0) {
                dist_j = 0.0;
                dist_i = (
                    mat->get(leftovers[1], leftovers[0]) + mat->get(leftovers[2], leftovers[0])
                ) * 0.5;
                if (dist_i < 0.0) dist_i = 0.0;
            }
        }
    
        theTree->child[0]->distance = dist_i;
        theTree->child[1]->distance = dist_j;
        theTree->child[2]->distance = dist_k;
    }else{
        for(i=0; i < numseqs; i++) {
            theTree->child[i] = nodes[i];
            nodes[i] = NULL;  // less two nodes: merge directly
        }

        if (numseqs == 2) {
            theTree->child[0]->distance = mat->get(1, 0) * 0.5;
            theTree->child[1]->distance = mat->get(1, 0) * 0.5;
        }
    }

    theTree->numnodes = nextfreenode;
    return theTree;
    
}