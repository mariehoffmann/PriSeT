# PriSeT
Tool for assisting the search of primer sequences for metabarcoding experiments. Given a taxonomic identifier, a reference database, and a genetic target region, PriSeT identifies conserved sections suitable for PCR primers such that taxonomic coverage and separation are optimal.

## Requirements

| **Platform**                       | **Details**            | **Tested** |
|:---------------------------------: | :--------------------: | :-----------: |
| <img src="./.github/Linux.svg" width="100" height="100" /> | `Linux 64 bit` | - |
| <img src="./.github/MacOS.svg" width="100" height="100" /> | `Mac OS 64 bit` | High Sierra 10.12.6 |

### R
If not installed on your system yet, install `R` via your standard package manager. However, for MacOS I recommend not to install via `port`, but download the binaries from [CRAN](https://cran.r-project.org/bin/macosx), because I ran into installation errors when when trying to install `igraph` and others in an interactive `R` session. If you install the package from `cran.r-project.org`, open the `R.app`, go to the package installer (under `Packages & Data`), search for the below listed packages and click the install button.
If you use R in terminal, start an inter session, install the required `R` packages `shiny` and `DT` for table output, and `treemap` and `d3treeR` for an interactive tree map plot.
```shell
$ R
> install.packages("shiny")
> install.packages("DT")
> install.packages("igraph")
> install.packages("treemap")
```
To make the treemap interactive you need `d3treeR` hosted on github. To download it in an interactive session, follow these steps:
```R
> install.packages("devtools")
> library(devtools)
> install_github("d3treeR/d3treeR")
```

```shell
wget https://ftp.gnu.org/pub/gnu/libiconv/libiconv-1.15.tar.gz
tar -zxvf libiconv-1.15.tar.gz
cd libiconv-1.15
./configure --prefix=/usr/local
make
make install
```

### Taxonomic Tree and Library
In order to explore the potential primer sequences hierarchically w.r.t. an existing taxonomy, the taxonomic node identifier (taxid) needs to be related to reference sequences.
You can use the tool [tactac](https://github.com/mariehoffmann/tactac) to create
a subset of your reference library. Calling `python tactac.py --subtree <taxid>` (plus password)
will create three files:
  * Taxonomic subtree as tuples in csv format: `/subset/<taxid>/root_<taxid>.tax`
  * Taxonomic map of taxids assigned directly to accessions: `/subset/<taxid>/root_<taxid>.acc`
  * Library of all sequences under `<taxid>`: `/subset/<taxid>/root_<taxid>.fasta`

## Setup
Clone PriSeT and its submodules
```shell
git clone --recurse-submodules https://github.com/mariehoffmann/PriSeT.git
```


### Compilation
  1. Create a build directory and compile with `cmake`
  ```shell
  cmake -DCMAKE_BUILD_TYPE=Debug ../PriSeT
  ```
  2. Build
  ```shell
  make -j
  ```

  3. Call binary with the library directory containing:
    * Taxonomy `root_<taxid>.tax`
    * Taxid to accessions map `root_<taxid>.acc`
    * Library with reference sequences `root_<taxid>.fasta`

 Assume `lib_dir` and `work_dir` as your library and working directories. The latter one will store the FM index, mappings, the shiny app as R script and the input table for the shiny app.
```shell
 ./priset <lib_dir> <work_dir> [--skip-idx]
 ```

### Unit Tests

```shell
cd ../build
cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc -DCMAKE_CXX_COMPILER=/usr/local/bin/g++ ../PriSeT/tests/
make
```

## References
   [1] Pockrandt, C., Alzamel, M., Iliopoulos, C. S., Reinert, K.. GenMap: Fast and Exact Computation of Genome Mappability. bioRxiv, presented on RECOMB-Seq, 2019.

   [2] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178

   [3] PCR primer design constraints: www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
