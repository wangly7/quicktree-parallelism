#include "tree.hpp"
#include "distancemat.hpp"

struct Tree* neighbour_joining_buildtree(DistanceMatrix* mat, const std::vector<std::string>& identifiers);