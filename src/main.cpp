#include <iostream>
#include <boost/program_options.hpp>


namespace po = boost::program_options;


int main(int argc, char** argv) {

    std::string matrixFilename;
    uint32_t maxTaxa;

    po::options_description desc{"Options"};
    desc.add_options()
    ("distanceMatrix, m", po::value<std::string>(&matrixFilename)->required(), "Input distance matrix in PHYLIP format [REQUIRED].")
    ("maxTaxa, N", po::value<uint32_t>(&maxTaxa)->default_value(1000), "Maximum number of taxa to read from the distance matrix file.")
    ("help,h", "Print help messages");

    po::options_description allOptions;
    allOptions.add(desc);
    
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).options(allOptions).run(), vm);
	po::notify(vm);
    } catch(std::exception &e){
        std::cerr << desc << std::endl;
        exit(1);	
    }
}
