#include <boost/program_options.hpp>
#include <fstream>
#include <cstdio>
#include "matcal.hpp"
#include "distancemat.hpp"
#include "buildtree.hpp"
#include "timer.hpp"

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
    

    // DistanceMatrix mat = readPhylipDistanceMatrix(fp); 
    // DistanceMatrix computed_tree = compute_tree(mat);
    // printDistanceMatrix(mat);

    std::vector<std::string> identifiers;
    
    timer.Start();
    fprintf(stdout, "Reading distance matrix from file.\n");
    DistanceMatrix mat = readPhylipDistanceMatrix(fp, identifiers);
    fprintf(stdout, "Completed in %ld msec \n\n", timer.Stop());

    FILE* out = stdout;

    if (!outputFilename.empty()) {
        out = fopen(outputFilename.c_str(), "w");
        if (!out) {
            fprintf(stderr, "ERROR: cant open file: %s\n", outputFilename.c_str());
            exit(1);
        }
    }

    Tree* myTree = neighbour_joining_buildtree(&mat, identifiers);
    write_newhampshire_Tree(out, myTree);

    // test distance matrix update
    uint32_t row=0;
    uint32_t col=0;

    fprintf(stdout, "Before update matrix: \n");
    printDistanceMatrix(&mat, row, col);

    DistanceMatrix new_mat = update_matrix(&mat);
    fprintf(stdout, "Afte update matrix: \n");
    printDistanceMatrix(&new_mat, row, col);
    return 0;

}
