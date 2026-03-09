#include <boost/program_options.hpp>
#include <fstream>
#include <vector>
#include "distancemat.hpp"
#include "buildtree.cuh"
#include "timer.hpp"
#include "tree.hpp"

namespace po = boost::program_options;


int main(int argc, char** argv) {
    Timer timer;
    std::string matrixFilename;
    std::string outputFilename;
    uint32_t maxTaxa;

    po::options_description desc{"Options"};
    desc.add_options()
    ("distanceMatrix,m", po::value<std::string>(&matrixFilename)->required(), "Input distance matrix in PHYLIP format [REQUIRED].")
    ("output,o", po::value<std::string>(&outputFilename), "Output file in the Newick/New-Hampshire format.")
    ("maxTaxa,N", po::value<uint32_t>(&maxTaxa)->default_value(1000), "Maximum number of taxa to read from the distance matrix file.")
    ("help,h", "Print help messages");

    po::options_description allOptions;
    allOptions.add(desc);
    
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(allOptions).run(), vm);
	    po::notify(vm);
    } catch(std::exception &e){
        std::cerr << "Argument error: " << e.what() << "\n\n";
        std::cerr << desc << std::endl;
        exit(1);	
    }

    std::ifstream fp(matrixFilename);
    if (!fp){
        fprintf(stderr, "ERROR: Cannot open file: %s\n", matrixFilename.c_str());
    }

    FILE* out = stdout;
    if (!outputFilename.empty()) {
        out = fopen(outputFilename.c_str(), "w");
        if (!out) {
            fprintf(stderr, "ERROR: cannot open file: %s\n", outputFilename.c_str());
            exit(1);
        }
    }

    // Print GPU information
    fprintf(stdout, "Printing GPU device properties.\n");
    printGpuProperties();
    
    // read distance matrix
    std::vector<std::string> identifiers;
    timer.Start();
    fprintf(stdout, "Reading distance matrix from file.\n");
    DistanceMatrix mat = readPhylipDistanceMatrix(fp, identifiers);
    fprintf(stdout, "Completed in %ld msec \n\n", timer.Stop());

    GpuTree TreeBuilder(mat.data, identifiers, mat.size);

    // launch tree building process
    timer.Start();
    TreeBuilder.build();
    TreeBuilder.clearAndReset();

    // write output tree
    write_newhampshire_Tree(out, TreeBuilder.theTree);

    fprintf(stdout, "\nProgram completed in %ld ms\n\n", timer.Stop());
    return 0;

}
