# SpydrPick

## Get SpydrPick

#### Bioconda
```
conda install -c bioconda spydrpick
```
See detailed instructions [here](README.md/#the-SpydrPick-Conda-package)

#### Docker
```
docker pull quay.io/biocontainers/spydrpick
```

#### C++ source code
```
git clone --recursive https://github.com/santeripuranen/SpydrPick.git
```
See how to compile [here](README.md/#building-SpydrPick-from-source).


## About

SpydrPick is a command line tool for performing direct coupling analysis of aligned categorical datasets. It constructs a coupling graph based on thresholded node-pair (edge) scoring,
followed by pruning of edges that represent indirect couplings.


## Using SpydrPick

### What does it do?

A typical run of `SpydrPick` will perform the following steps: 1) parse an input alignment, 2) apply default position filtering rules, 3) determine sample weights designed to correct for population structure, 4) estimate an MI threshold in order to save a pre-defined number of top ranking pairs, 5) evaluate all pairwise position-position MI values, 6) estimate outlier and extreme outlier thresholds (the threshold of likely significant signal) based on the data, 7) apply the ARACNE post-processing step to the MI graph, and finally 9) output lists of the strongest edges (links) found.


### Basic usage

Use `SpydrPick -h` or `SpydrPick --help` to get a list of available command line options.

To run SpydrPick with default settings use:
```
SpydrPick -v <name of input genome alignment file>
```
where the input alignment should be in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).

The default settings are designed to be overall reasonable and can usually be left as they are. Notable exceptions are the *linkage disequilibrium (LD) threshold distance* option `--ld-threshold=<int:distance>` and whether the alignment represents a *circular* (default) or a *linear* chromosome/genome (add `--linear-genome`). The level and reach of LD depends on factors such as the recombination characteristics of the organism under study. For bacterial genomes typical values of `--ld-threshold` are in the 500-20000 bp range. The *outlier* and *extreme outlier* thresholds calculated by SpydrPick (and reported in the console output) are affected by these options as is, consequently, the contents of the outlier listing; the MI and ARACNE steps are not affected *per se*. It is in general safer to err to longer, more conservative distance estimates here, as this will have less of an impact on *outlier* and *extreme outlier* thresholds than too short distance cutoffs.


### Advanced usage

Default filtering will extract positions with *more than* 1 allele (not counting gaps), where the less frequent will be of *at least* 1% frequency, and with *at most* 15% gap frequency. These criteria are usually very effective at reducing the number of positions included in the subsequent computations, with little or no effect on SpydrPick's ability to find relevant associations, since columns with highly skewed allele distributions are non-informative in the information theoretical sense. The filtering criteria can be changed with the `--maf-threshold=<float:[0,1]>` and `--gap-threshold=<float:[0,1]>` command line options, or disabled entirely with the `--no-filter-alignment` flag.

If the input alignment comes prefiltered, or the columns have non-unity stride for some other reason, the `--mappings-list=<string:filename>` option can be used to inform SpydrPick of the actual position indices. This will enable correct distances for LD threshold evaluations and will ensure that position indices are correct in the output. `--mappings-list` expects a file containing a white-space delimited list of indices the same size as the number of data columns in the input and matching each column in order from left to right. If the input genome is not linear and the full genome was not supplied as input then the correct genome size should be specified using the `--genome-size=<int:size>` option. Input and output indexing bases can be controlled separately using the `--input-indexing-base=<int:base>` and `--output-indexing-base=<int:base>` options. 

SpydrPick will by default attempt to correct for population structure by assigning weights across samples in the input alignment. This weighting can modified by specifying with `--sample-reweighting-threshold=<float:[0,1]>` the threshold sequence identity for when two samples/sequences are considered to be equal. Weighting can be turned off altogether with the `--no-sample-reweighting` option -- all samples will then have equal (=1) weight. SpyderPick can also use sample weights supplied by the user with the `--sample-weights=<string:filename>` option. `--sample-weights` expects a file containing a white-space delimited list of values the same size as the number of samples in the input and matching each sequence from top to bottom as they appear in the input alignment.

SpydrPick will by default return enough top ranking position pairs that these with certainty encompass all statistically significant interactions. However, an arbitrary number of top MI values can be requested using the `--mi-values=<int:#edges>` option.

## Deciphering SpydrPick output

The main SpydrPick output file (`*.spydrpick_couplings.?-based.*edges`) contains a white space delimited, descending order list of MI values and pairs of position indices (using *1-based indexing* by default; control this with `--output-indexing-base`) numbered according to the columns (left to right) in the input alignment. The fields in the output are `[pos1 pos2 genome_distance ARACNE MI]`.

Position pairs with MI values above the *outlier threshold* and further apart than `--ld-threshold` can be found also in the file titled `*.outliers`, along with some additional data. Here the output fields are `[pos1, pos2, genome_distance, ARACNE, MI, MI_wo_gaps, gap_effect, extreme_outlier]`, where `gap_effect` is calculated as `(1-MI_wo_gaps/MI)*100` and `extreme_outlier` indicates values that surpass the *extreme_outlier_threshold*. A high value of `gap_effect` warrants a closer look at the allele compositions of the positions involved, as the high MI value assigned to the edge appears to be gap-driven, indicating a possible alignment-associated artifact.


## Installation details

### The SpydrPick Conda package

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


### Building SpydrPick from source

In order to compile the `SpydrPick` binary, go to the `build` directory (create one if necessary; in-source builds are strongly discouraged) and give these commands:
```
cmake ..
make SpydrPick
```
This should set up the CMake project, compile the binary and place it into the `build/bin` directory. If not, then take a look at [Compile-time dependencies](README.md/#compile--time-dependencies).

The `SpydrPick` binary will by default be statically linked, except for [TBB](https://www.threadingbuildingblocks.org/) that can only be linked dynamically. Installing SpydrPick to another location is as easy as copying the binary. Note that the [TBB](https://www.threadingbuildingblocks.org/) runtime library will need to be present on the host system.


#### Compile-time dependencies

`SpydrPick` is written in C++ and wrapped into a [CMake](https://cmake.org/) project. It relies on several common external libraries. Your build environment must satisfy the following requirements and compile-time dependencies:

* A C++14 compliant compiler (development was done using the [GNU C++ compiler](https://gcc.gnu.org/))
* [CMake](https://cmake.org/)
* [Boost](https://www.boost.org/)
* [Intel(R) Threading Building Blocks (TBB)](https://www.threadingbuildingblocks.org/)

You may need to set the [`CMAKE_MODULE_PATH`](https://cmake.org/cmake/help/latest/variable/CMAKE_MODULE_PATH.html) environment variable in order for `CMake` to find all relevant packages.


## Cite

SpydrPick was developed as part of an academic project. Please cite:

* Johan Pensar, Santeri Puranen, Neil MacAlasdair, Juri Kuronen, Gerry Tonkin-Hill, Maiju Pesonen, Brian Arnold, Yingying Xu, Aleksi Sipola, Leonor Sanchez-Buso, John A Lees, Claire Chewapreecha, Stephen D Bentley, Simon R Harris, Julian Parkhill, Nicholas J Croucher, Jukka Corander. **Genome-wide epistasis and co-selection study using mutual information.** *bioRxiv* 2019 [**doi:** 10.1101/523407](https://doi.org/10.1101/523407)

* https://github.com/santeripuranen/SpydrPick
