# PriSeT
Tool for assisting the search of primer sequences for metabarcoding experiments. Given a taxonomic identifier, a reference database, and a genetic target region, PriSeT identifies conserved sections suitable for PCR primers such that taxonomic coverage and separation are optimal.

## Requirements

| **Platform**                       | **Details**            | **Tested** |
|:---------------------------------: | :--------------------: | :-----------: |
| <img src="./.github/Linux.svg" width="100" height="100" /> | `Linux 64 bit` | - |
| <img src="./.github/MacOS.svg" width="100" height="100" /> | `Mac OS 64 bit` | High Sierra 10.12.6 |

### PostgreSQL
For setting up the database you need PostgreSQL (arbitrary version) and the server. Use e.g. `apt-get` or `port` to install the latest distributions.

On Linux:
```
sudo apt-get install postgresql postgresql-contrib
```

On MacOS
```
sudo port install postgresql11 postgresql11-server
```

## Setup
Clone PriSeT and its submodules
```
git clone --recurse-submodules https://github.com/mariehoffmann/PriSeT.git
```


### Compilation
  1. Create a build directory and compile `priset` from there
  ```
  g++ priset.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -o priset
  ```
  2. Call binary with path to `genmap` binary, the source directory (location of fasta and taxonomy files) and a working directory for temporary output
  ```
 ./priset <path_genmap_bin> <src_dir> <work_dir>
  ```
