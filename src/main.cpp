#include <boost/program_options.hpp>
#include <fstream>
#include "distancemat.hpp"
#include "matcal.cpp"
#include "timer.hpp"

namespace po = boost::program_options;


int main(int argc, char** argv) {
    Timer timer;
    std::string matrixFilename;
    uint32_t maxTaxa;

    po::options_description desc{"Options"};
    desc.add_options()
    ("distanceMatrix,m", po::value<std::string>(&matrixFilename)->required(), "Input distance matrix in PHYLIP format [REQUIRED].")
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
    
    DistanceMatrix mat = readPhylipDistanceMatrix(fp); 
    DistanceMatrix computed_tree = compute_tree(mat);
    printDistanceMatrix(mat);

    std::vector<std::string> identifiers;
    
    timer.Start();
    fprintf(stdout, "Reading distance matrix from file.\n");
    DistanceMatrix mat = readPhylipDistanceMatrix(fp, identifiers);
    fprintf(stdout, "Completed in %ld msec \n\n", timer.Stop());
}
