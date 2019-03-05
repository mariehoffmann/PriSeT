# PriSeT
Tool for assisting the search of primer sequences for metabarcoding experiments. Given a set of higher order taxonomies, a reference database, and a region to search in, PriSeT tries to identify conserved sections suitable for PCR primers such that taxonomic coverage and separation are optimal.

## Requirements
  * OS: Linux, MacOS
  * Postgresql

## Setup
In the background PriSet uses database queries for walking in the taxonomy tree and collecting relevant sequences. Therefore you first need to setup the NCBI Taxonomy database. NCBI provides all data files necessary to build locally a database according to the scheme below.

### Compiling Unit Test
  1. Create a build directory, e.g. in your home directory and change into it
  ```
  mkdir -p ~/builds/priset
  cd ~/builds/priset
  ```
  2. Compile complete test suite
  ```
  cmake -DCMAKE_CXX_COMPILER=g++-mp-7 $PRISET_DIR/test/unit
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
