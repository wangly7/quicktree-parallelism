#include <cstdint>
#include <iostream>


struct TNode {
  TNode* left;
  TNode* right;
  TNode* parent;  
  
  double distance = 0.0;
  uint32_t nodenumber;
  std::string identifier;

  TNode()
    : left(nullptr), right(nullptr), parent(nullptr),
      distance(0.0), nodenumber(0), identifier("") {}
};

struct Tree {
  TNode *child[3];
  uint32_t numnodes;

  Tree(): child{nullptr, nullptr, nullptr}, numnodes(0) {}
};


void write_newhampshire_Tnode(FILE* out, TNode* node);
void write_newhampshire_Tree(FILE* out, Tree* tree);
