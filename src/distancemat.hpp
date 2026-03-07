#pragma once


#include <vector>
#include <cstdint>
#include <iostream>
#include <sstream>
#define MAX_PHYLIP_NAME_LEN 100


struct DistanceMatrix {
  uint32_t size;
  std::vector<double> data;

  DistanceMatrix(uint32_t n) : size(n), data(n*(n-1)/2) {}

  inline uint32_t index(uint32_t i, uint32_t j) const {
    if (i < j) std::swap(i, j);
    return i * (i - 1) / 2 + j;
  }

  inline double get(uint32_t i, uint32_t j) const {
    return data[index(i, j)];
  }

  inline void set(uint32_t i, uint32_t j, double value) {
    data[index(i, j)] = value;
  }
};

void printDistanceMatrix(const DistanceMatrix* distanceMat, uint32_t row, uint32_t column);
DistanceMatrix readPhylipDistanceMatrix(std::istream& input, std::vector<std::string>& identifiers);


