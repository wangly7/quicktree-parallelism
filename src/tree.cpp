#include "tree.hpp"


void write_newhampshire_Tnode(FILE* out, TNode* node) {
    if (node != NULL) {
        if (node->left == NULL && node->right ==NULL) {
            /* leaf node */
            if (node->identifier.empty()) {
                fprintf(stderr, "Fatal ERROR: encountered a leaf node with no identifier info");
                exit(1);
            }else {
                fprintf(out, "%s:%.5f", node->identifier.c_str(), node->distance);
            }
        }else if (node->left != NULL && node->right != NULL) {
            /* internal node: node that been merged */
            fprintf(out, "(");
            write_newhampshire_Tnode(out, node->left);
            fprintf(out, ",");
            write_newhampshire_Tnode(out, node->right);
            fprintf(out, "):%.5f", node->distance);
        }else {
            fprintf(stderr, "ERROR: unvalid intrenale node\n");
            exit(1);
        }
    }
}

void write_newhampshire_Tree(FILE* out, Tree* tree) {
    if (tree != NULL) {
        if (tree->child[0] != NULL) {
            if (tree->child[1] == NULL) {
                TNode* root = tree->child[0];
                /* draw rooted tree */
                if(root->left != NULL && root->right != NULL) {
                    fprintf(out, "(");
                    write_newhampshire_Tnode(out, root->left);
                    fprintf(out, ",");
                    write_newhampshire_Tnode(out, root->right);
                    fprintf(out, ");");
                }else if (root->left == nullptr && root->right == nullptr){
                    fprintf(stderr, "ERROR: Could not build a tree with a single sequence");
                    exit(1);
                }
            }else {
                /* build neighbor joing unrooted tree */
                fprintf(out, "(");
                write_newhampshire_Tnode(out, tree->child[0]);
                fprintf(out, ",");
                write_newhampshire_Tnode(out, tree->child[1]);
                if (tree->child[2] != NULL) {
                    fprintf(out, ",");
                    write_newhampshire_Tnode(out, tree->child[2]);
                }
                fprintf(out, ");");
            }
        }
    }
}