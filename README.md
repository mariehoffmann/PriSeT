# PriSeT
Tool for assisting the search of primer sequences for metabarcoding experiments. Given a taxonomic identifier, a reference database, and a genetic target region, PriSeT identifies conserved sections suitable for PCR primers such that taxonomic coverage and separation are optimal.


| **Platform**                       | **Details**            | **Binaries** |
|:---------------------------------: | :--------------------: | :-----------: |
| <img src="./.github/MacOS.svg" width="100" height="100" /> | `Mac OS 64 bit` | [Link to Build](https://github.com/mariehoffmann/PriSeT/blob/master/binaries/priset_v10_darwin) |

[comment]: <> (| <img src="./.github/Linux.svg" width="100" height="100" /> | `Linux 64 bit` | Link to Build |)


## Usage 

### Usecase: Precompiled Binary

Download one of the above linked binaries (Ubuntu build will be provided soon) and call it from command line with a directory containing a FASTA file of the reference sequences. The binary is based on `apps/solver_fast.cpp`. In case you want to edit the app and change the primer properties via setter functions provied by `types/PrimerConfig.hpp`, follow the compilation instructions in the next section.
The FASTA file should not be larger than 500 MB. Index building and primer search are decoupled to allow for tweaking the constraints that affect the primer pair result set.

#### Build FM-Index
First an FM-index is built on the FASTA file and stored under the user-given directory. Note, that the FM-index consumes about four times more space than the input file. The index for a specific library has to be built only once. Create `work_dir` in advance.
```
./priset -i -l <path_to_fasta> -w <work_dir|idx_dir>
```

#### Discover Primer Pairs
Primer search is triggered via the `-s` flag. The results are output in csv format; the first 20 pairs are displayed in terminal.

[comment]: <> (Optionally, a couple of primer sequence parameters can also be set in an experimental configuration file. If omitted the default parameters defined in `PrimerConfig.hpp` are chosen.)

```
./priset -s -l <path_to_fasta> -w <work_dir>
```

## Usecase: Userwritten Apps

In case you want to modify the way how primer pairs are searched and processed, you can write and compile your own applications. 

### Requirements

  - SDSL-Lite, get it from here: [xxsds/sdsl-lite](https://github.com/xxsds/sdsl-lite)
  - Cmake 3.7 or higher
  - GNU C++ Compiler xx or higher

#### Setup
Clone PriSeT and its submodules
```shell
git clone --recurse-submodules https://github.com/mariehoffmann/PriSeT.git
```

#### Compilation

User-defined apps should go into the `apps` folder, and the therein contained `CMakeLists.txt` file modified. In the following we assume your application is called `solver_fast.cpp`. Add the following entry to the `apps/CMakeLists.txt`:
```
priset_app_macro(solver_fast.cpp)
```
The macro will create a target and link it against the required libraries.


  1. Create a build directory and compile with `cmake`
  ```shell
cmake ../PriSeT/apps -DCMAKE_BUILD_TYPE=Debug -B .
  ```
  2. Build
  ```shell
  make -j
  ```

  3. Run to build FM-index once 
  ```shell
./solver_fast -i -l <lib_dir> -w <work_dir>
  ```

  4. Discover primers 
  ```shell
./solver_fast -s -l <lib_dir> -w <work_dir>
  ```

### Unit Tests

Add new unit tests under test/unit and run cmake 
```shell
cmake ../PriSeT/test/unit/ -DGENMAP_NATIVE_BUILD=ON -DCMAKE_BUILD_TYPE=Debug -B .
```

## References
   [1] Pockrandt, C., Alzamel, M., Iliopoulos, C. S., Reinert, K.. GenMap: Fast and Exact Computation of Genome Mappability. bioRxiv, presented on RECOMB-Seq, 2019.

   [2] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178

   [3] PCR primer design constraints: www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
