#include "distancemat.hpp"

/**
 * Construct 1D distance matrix array from PHYLIP format distance matrix file
 */
DistanceMatrix readPhylipDistanceMatrix(std::istream& input, std::vector<std::string>& identifiers) {
    uint32_t size = 0;
    input >> size;
    DistanceMatrix mat(size);

    if (!size || size <= 0) fprintf(stderr, "ERROR: Invalid PHYLIP header");

    identifiers.clear();
    identifiers.reserve(size);
    std::string line;
    // start from third line
    std::getline(input, line);
    
    for(uint32_t i = 0 ; i < size; i++) {
        if (!std::getline(input, line)) fprintf(stderr, "ERROR: Unexpected EOF.");

        std::string identifier;
        std::istringstream value(line);
        // first column: specie's name
        value >> identifier;
        identifiers.push_back(identifier);
        for (uint32_t j = 0; j < i; j++) {
            double distance;
            if (!(value >> distance)) {
                // PHYLIP allow matrix format like:
                // T1 0.2 0.3 0.4
                //    0.5
                if(!(input >> distance)) fprintf(stderr, "ERROR: Missing distance at position (%d, %d).\n", i, j);
                std::getline(input, line);
            }
            mat.set(i, j, distance);
        }
    }
    // printDistanceMatrix(&mat, size, size);
    return mat;
}


/**
 * Prints the given distance matrix.
 */
void printDistanceMatrix(const DistanceMatrix* mat, uint32_t row, uint32_t column) {
  for(row=0; row < mat->size; row++) {
    fprintf(stdout, "%5d", row);
    for(column=0; column < row; column++)
	    fprintf(stdout, "%10.5f", mat->get(row, column));
        fprintf(stdout, "\n");
  }
}