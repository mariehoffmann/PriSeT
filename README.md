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

### Instructions for Local Database
Assume `data_dir` your directory for downloading data files.
  1. In the terminal go to your data directory cd `data_dir`
  2. Download one of the archived `new_taxdump` files, e.g. the zipped one with the md5 checksum and README file:

  ```bash
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip.md5
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/taxdump_readme.txt
  ```
  3. Compute md5 hash to check if nothing went wrong (`md5sum` on Linux) and compare to provided one (output should be 1)

    ```bash
    expr `cat new_taxdump.zip.md5 | cut -d' ' -f1` = `md5 new_taxdump.zip | cut -d'=' -f2 | xargs`
    ```
  4. Download the accession number to taxonomic ID resolution file from GenBank:
  ```bash
  wget
  ```
  4. Download Postgresql Database Server
  5. Download your favorite Database Client
