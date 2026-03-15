# Quicktree Parallelism
We parallelize the implementation of the Neighbor-Joning (NJ) algorithm based on QuickTree.

## Table of Contents
- [Quicktree Parallelism](#quicktree-parallelism)
  - [Table of Contents](#table-of-contents)
  - [Dataset](#dataset)
  - [ Requirements](#-requirements)
  - [Build Instructions](#build-instructions)
  - [Run the Program](#run-the-program)
  - [Execution Workflow on DSMLP Server](#execution-workflow-on-dsmlp-server)

## Dataset
Distance matrix with size 1K, 5K, and 10K are available within the data folder in the repository. Distance matrix in 20K size is available on [Google Drive](https://drive.google.com/file/d/1fgA6SpcN3XCDpbHKMmuSvPujOlT3G2lN/view?usp=sharing). Please download and move it within the data folder.

## <a name="req"></a> Requirements
The following dependencies are required:
- C++ compiler (GCC / Clang)
- CMake
- CUDA Toolkit


## Build Instructions
Create a build directory and compile the project:

```bash
mkdir -p build
cd build
cmake ..
make -j4
```

This will generate the executable:

```bash
quicktree
```

## Run the Program

Run the program using a distance matrix in PHILIP format as input (output will be printed to the console):

```bash
./quicktree -m ../data/matrix_10k.phy
```

Run the program using a distance matrix in PHILIP format as input (output will be written to the file):
```bash
./quicktree -m ../data/matrix_10k.phy -o tree.nwk
```

The program arguments are as follows:

| Argument | Description |
|----------|-------------|
| `-m` | input distance matrix (must in PHILIP format)|
| `-o` | output phylogenetic tree in Newick format |


## Execution Workflow on DSMLP Server 
1. Clone the project repository in HOME directory, and decompress the data files
   ```bash
    git clone https://github.com/wangly7/quicktree-parallelism.git
    cd quicktree-parallelism/data
    unzip matrix_10k.phy.zip
    cd ..
   ```

2. Run the code
     ```bash
    /opt/launch-sh/bin/launch.sh -v a30 -c 8 -g 1 -m 8 -i yatisht/ece213-wi26:latest -f ./quicktree-parallelism/run-commands.sh
    ```
