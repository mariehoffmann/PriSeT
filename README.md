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
```shell
sudo apt-get install postgresql postgresql-contrib
```

On MacOS
```shell
sudo port install postgresql11 postgresql11-server
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
  cmake ../PriSeT/ -DGENMAP_NATIVE_BUILD=OFF
  ```
  2. Call binary with the library directory containing:
    * Taxonomy `root_<taxid>.tax`
    * Taxid to accessions map `root_<taxid>.acc`
    * Library with reference sequences `root_<taxid>.fasta`

 In the working directory the FM index and mappings are stored.
 ```shell
 ./priset <lib_dir> <work_dir>
 ```


 ## References

   [1] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178

   [2] PCR primer design constraints: www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
