# SpydrPick

## Get the package from Bioconda
```
conda -c bioconda install spydrpick
``
See detailed instructions [here](README.md/#the-SpydrPick-Conda-package


## Get the Docker image
```
docker pull quay.io/biocontainers/spydrpick
```


## Get the code
```
git clone --recursive https://github.com/santeripuranen/SpydrPick.git
```
See how to compile [here](README.md/#building-SpydrPick).


## About

Spydrpick is a command line tool for performing direct coupling analysis of aligned categorical datasets. It constructs a coupling graph based on thresholded node-pair (edge) scoring,
followed by pruning of edges that represent indirect couplings.


## The SpydrPick Conda package

An easy way to install SpydrPick is through conda, which is most easily accessed by first installing [miniconda](https://conda.io/miniconda.html). SpydrPick can then be installed by running:
```
conda install spydrpick
```

If the package cannot be found you will need to add the necessary channels:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```


## Building SpydrPick

In order to compile the `SpydrPick` binary, go to the `build` directory (create one if necessary; in-source builds are strongly discouraged) and give these commands:
```
cmake ..
make SpydrPick
```
This should set up the CMake project, compile the binary and place it into the `build/bin` directory. If not, then take a look at [Compile-time dependencies](README.md/#compile--time-dependencies).

The `SpydrPick` binary will by default be statically linked, except for [TBB](https://www.threadingbuildingblocks.org/) that can only be linked dynamically. Installing SpydrPick to another location is as easy as copying the binary (given that the [TBB](https://www.threadingbuildingblocks.org/) runtime library is properly installed on your system).


### Compile-time dependencies

`SpydrPick` is written in C++ and wrapped into a [CMake](https://cmake.org/) project. It relies on several common external libraries. Your build environment must satisfy the following requirements and compile-time dependencies:

* A C++14 compliant compiler (development was done using the 7-series [GNU C++ compiler](https://gcc.gnu.org/))
* [CMake](https://cmake.org/)
* [Boost](https://www.boost.org/)
* [Intel(R) Threading Building Blocks (TBB)](https://www.threadingbuildingblocks.org/)

You may need to set the [`CMAKE_MODULE_PATH`](https://cmake.org/cmake/help/latest/variable/CMAKE_MODULE_PATH.html) environment variable in order for `CMake` to find all relevant packages.


## Basic usage

Use `SpydrPick -h` or `SpydrPick --help` to get a list of available command line options.

To run SpydrPick with default settings use:
```
SpydrPick -v <name of input genome alignment file>
```
where the input alignment should be in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).

The default settings are designed to be overall reasonable and can usually be left as they are. Notable exceptions to this are the *linkage disequilibrium (LD) threshold distance* (`--ld-threshold`) and whether the alignment represents a *circular* (default) or *linear* chromosome/genome (`--linear-genome`). The level and reach of LD depends on factors such as the recombination characteristics of the organism under study. For bacterial genomes typical values of `--ld-threshold` are in the 500-20000 range. The *outlier* and *extreme outlier* thresholds calculated by `SpydrPick` (and reported in the console output) are affected by these options as is, consequently, the contents of the outlier listing; the MI and ARACNE steps are not affected *per se*.

A typical run of `SpydrPick` will perform the following steps: 1) parse the input alignment, 2) apply default SNP filtering rules, 3) determine sample weights designed to correct for population structure, 4) evaluate all pairwise MI values, 5) apply the ARACNE post-processing step to the resulting graph, and 6) finally output lists of the strongest edges (links) found.

Default filtering will extract positions with *more than* 1 allele (not counting gaps), *at least* 1% minor allele frequency and *at most* 15% gap frequency. These criteria are usually very effective at reducing the number of positions included in the subsequent computations, with little or no effect on SpydrPick's ability to find relevant associations, since columns with highly skewed allele distributions are non-informative in the information theoretical sense. The filtering criteria can be changed with the `--maf-threshold` and `--gap-threshold` command line options, or disabled completely with the `--no-filter-alignment` flag.

## Deciphering SpydrPick output

The main SpydrPick output file (`*.spydrpick_couplings.?-based.*edges`) contains a white space delimited, descending order list of MI values and pairs of position indices (using *1-based indexing* by default; control this with `--output-indexing-base`) numbered according to the columns (left to right) in the input alignment. The fields in the output are `[pos1 pos2 genome_distance ARACNE MI]`.

Position pairs with MI values above the *outlier threshold* and further apart than `--ld-threshold` can be found also in the file titled `*.outliers`, along with some additional data. Here the output fields are `[pos1, pos2, genome_distance, ARACNE, MI, MI_wo_gaps, gap_effect, extreme_outlier]`, where `gap_effect` is calculated as `(1-MI_wo_gaps/MI)*100` and `extreme_outlier` indicates values that surpass the *extreme_outlier_threshold*. A high value of `gap_effect` warrants a closer look at the allele compositions of the positions involved, as the high MI value assigned to the edge appears to be gap-driven, indicating a possible alignment-associated artifact.

