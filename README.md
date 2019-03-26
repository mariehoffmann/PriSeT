# PriSeT
Tool for assisting the search of primer sequences for metabarcoding experiments. Given a taxonomic identifier, a reference database, and a genetic target region, PriSeT identifies conserved sections suitable for PCR primers such that taxonomic coverage and separation are optimal.

## Requirements
  * OS: Linux, MacOS
  * genmap 

## Setup
In the background PriSet uses database queries for walking in the taxonomy tree and collecting relevant sequences. Therefore you first need to setup the NCBI Taxonomy database. NCBI provides all data files necessary to build locally a database according to the scheme below.

### Compilation
  1. Create a build directory and compile `priset` from there
  ```
  g++ priset.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -o priset
  ```
  2. Call binary with path to `genmap` binary, the source directory (location of fasta and taxonomy files) and a working directory for temporary output
  ```
 ./priset <path_genmap_bin> <src_dir> <work_dir>
  ```
