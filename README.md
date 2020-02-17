## Get SpydrPick

#### Bioconda
[![Install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/spydrpick/README.html) ![Supported platforms](https://anaconda.org/bioconda/spydrpick/badges/platforms.svg) [![Anaconda-Server Badge](https://anaconda.org/bioconda/spydrpick/badges/version.svg)](https://anaconda.org/bioconda/spydrpick)
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

SpydrPick is a command line tool for performing 'direct coupling analysis' of aligned categorical datasets. It constructs a coupling graph based on thresholded node-pair [Mutual Information](https://en.wikipedia.org/wiki/Mutual_information) (MI) scoring, followed by pruning of edges that represent indirect couplings using the [ARACNE](https://doi.org/10.1186/1471-2105-7-S1-S7) algorithm.


## Using SpydrPick

### What does it do?

A typical run of `SpydrPick` will perform the following steps: 1) parse an input alignment, 2) apply default position filtering rules, 3) determine sample weights designed to correct for population structure, 4) estimate an MI threshold in order to save a pre-defined number of top ranking pairs, 5) evaluate all pairwise position-position MI values, 6) estimate outlier and extreme outlier thresholds (the threshold of likely significant signal) based on the data, 7) apply the ARACNE post-processing step to the MI graph, and finally 9) output lists of the strongest edges (links) found.


### Basic usage

Use `SpydrPick -h` or `SpydrPick --help` to get a list of available command line options.

To run SpydrPick with default settings use:
```
SpydrPick -v <name of input genome alignment file>
```
where the input alignment should be in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format). Five distinct categories are currently supported: nucleotide symbols `A`, `C`, `G` and `T` are mapped as distinct categories, whereas *all other symbols* map to a single gap category. The FASTA parser is currently case-insensitive, i.e. lower-case and upper-case letters map to the same categories.

The default settings are designed to be overall reasonable and can usually be left as they are. Notable exceptions are the *linkage disequilibrium (LD) threshold distance* option `--ld-threshold=<int:distance>` and whether the alignment represents a *circular* (default) or a *linear* chromosome/genome (add `--linear-genome`). The level and reach of LD depends on factors such as the recombination characteristics of the organism under study. For bacterial genomes typical values of `--ld-threshold` are in the 500-20000 bp range. The *outlier* and *extreme outlier* thresholds calculated by SpydrPick (and reported in the console output) are affected by these options as is, consequently, the contents of the outlier listing; the MI and ARACNE steps are not affected *per se*. It is in general safer to err to longer, more conservative distance estimates here, as this will have less of an impact on *outlier* and *extreme outlier* thresholds than too short distance cutoffs.


### Advanced usage

Default filtering will extract positions with *more than* 1 allele (not counting gaps), where the second most frequent will be of *at least* 1% frequency, and with *at most* 15% gap frequency. The primary aim of position filtering is to reduce the computational effort required in the subsequent analysis stage, by focusing analysis on the subset of positions where *detectable* signal is more likely to exist; columns with highly skewed allele distributions are non-informative in the information theoretical sense. The default criteria are usually very effective at reducing the number of positions included in the subsequent computations, with little or no effect on SpydrPick's ability to find relevant associations. Note though that empirical sensitivity power analysis using simulated genomes (see supplementary information for [10.1093/nar/gkz656](https://doi.org/10.1093/nar/gkz656)) suggests that the default minor allele frequency threshold is quite conservative, given the typical samples sizes and sample diversity of current (anno 2019) bacterial genomics data, and it *should* in most cases be perfectly safe with a moderately stricter threshold of up to 5% without risk of missing any otherwise detectable signal at the MI stage (*caveat utilitor* though). The filtering criteria can be changed with the `--maf-threshold=<float:[0,1]>` and `--gap-threshold=<float:[0,1]>` command line options, or disabled entirely with the `--no-filter-alignment` flag.

If the input alignment comes prefiltered, or the columns have non-unity stride for some other reason, the `--mappings-list=<string:filename>` option can be used to inform SpydrPick of the actual position indices. This will enable correct distances for LD threshold evaluations and will ensure that position indices are correct in the output. `--mappings-list` expects a file containing a white-space delimited list of indices the same size as the number of data columns in the input and matching each column in order from left to right. If the input genome is not linear and the full genome was not supplied as input then the correct genome size should be specified using the `--genome-size=<int:size>` option. Input and output indexing bases can be controlled separately using the `--input-indexing-base=<int:base>` and `--output-indexing-base=<int:base>` options. 

SpydrPick will by default attempt to correct for population structure by assigning weights across samples in the input alignment. This weighting can be modified by specifying with `--sample-reweighting-threshold=<float:[0,1]>` the threshold sequence identity for when two samples/sequences are considered to be equal. Weighting can be turned off altogether with the `--no-sample-reweighting` option -- all samples will then have equal (=1) weight. SpyderPick can also use sample weights supplied by the user with the `--sample-weights=<string:filename>` option. `--sample-weights` expects a file containing a white-space delimited list of values the same size as the number of samples in the input and matching each sequence from top to bottom as they appear in the input alignment.

SpydrPick will by default return enough top ranking position pairs that these with certainty encompass all statistically significant interactions. However, an arbitrary number of top MI values can be requested using the `--mi-values=<int:#edges>` option.


## Deciphering SpydrPick output

The main SpydrPick output file `*.spydrpick_couplings.?-based.*edges` contains a white space delimited, descending order list of MI values and pairs of position indices (using *1-based indexing* by default; control this with `--output-indexing-base`) numbered according to the columns (left to right) in the input alignment. The fields in the output are `[pos1 pos2 genome_distance ARACNE MI]`.

Position pairs with MI values above the *outlier threshold* and further apart than `--ld-threshold` can be found also in the file titled `*.outliers`, along with some additional data. Here the output fields are `[pos1, pos2, genome_distance, ARACNE, MI, MI_wo_gaps, gap_effect, extreme_outlier]`, where `gap_effect` is calculated as `(1-MI_wo_gaps/MI)*100` and `extreme_outlier` indicates values that surpass the *extreme_outlier_threshold*. A high value of `gap_effect` warrants a closer look at the allele compositions of the positions involved, as the high MI value assigned to the edge appears to be gap-driven, indicating a possible alignment-associated artifact.


## Plotting the results

The SpydrPick output can be visualized as a [Manhattan plot](https://en.wikipedia.org/wiki/Manhattan_plot). We provide an easy-to-use [R](https://www.r-project.org/) script called [`gwes_plot.r`](gwes_plot.r) to create such a plot.

To use `gwes_plot.r`, edit the `data_full_filepath` and (optionally) the `outliers_full_filepath` fields in the script file and run it in R.

Fill in the `ld_dist` field in order to mark the *ld_threshold* with a vertical line in the plot.

The plot file dimensions can be changed by modifying the values in section `Plot sizes`.


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
This should set up the SpydrPick [CMake](https://cmake.org/) project, compile the binary and place it into the `build/bin` directory. If not, then take a look at [Compile-time dependencies](README.md/#compile--time-dependencies) and [CMake configuration](README.md/#configuring-the-environment-for-CMake).

The `SpydrPick` binary will by default be statically linked, except for [TBB](https://www.threadingbuildingblocks.org/) that can only be linked dynamically. Installing SpydrPick to another location is as easy as copying the binary. Note that the [TBB](https://www.threadingbuildingblocks.org/) runtime library must be present on the host system.


#### Compile-time dependencies

`SpydrPick` is written in C++ and wrapped into a [CMake](https://cmake.org/) project. It relies on several common external libraries. Your build environment must satisfy the following requirements and compile-time dependencies:

* A C++14 compliant compiler (development was done using the [GNU C++ compiler](https://gcc.gnu.org/), but [Clang/LLVM](https://clang.llvm.org/) should work just as well)
* [CMake](https://cmake.org/) (version 3.3 or later)
* [Boost](https://www.boost.org/)
* [Intel(R) Threading Building Blocks (TBB)](https://www.threadingbuildingblocks.org/)

#### Configuring the environment for CMake

CMake locates and sets up the external dependencies on the local system using script files named `Find<package>.cmake`, which on linux systems are typically found in `/usr/share/cmake/Modules` or in the `Modules` subdirectory of custom CMake builds/installs; the [`CMAKE_MODULE_PATH`](https://cmake.org/cmake/help/latest/variable/CMAKE_MODULE_PATH.html) environment variable may need to be set in order for CMake to locate all relevant `Find*.cmake` scripts.

The CMake script for locating Boost usually comes bundled with the CMake installation itself; however, this is not the case for TBB. If your system does not contain the `FindTBB.cmake` file and you get an [error message like this](https://github.com/santeripuranen/SpydrPick/issues/2#issue-405332114) then you can borrow the `FindTBB.cmake` script found over at Kitware's [VTK project](https://github.com/Kitware/VTK/blob/master/CMake/FindTBB.cmake) (save the file and point `CMAKE_MODULE_PATH` to the parent directory).

In addition to `CMAKE_MODULE_PATH`, the following (shell) environment variables may need to be set for CMake project configuration to work properly: [`CMAKE_CXX_COMPILER`](), `BOOST_ROOT` and `TBB_ROOT`.

In our experience even setting all of the above does not always guarantee success on all systems. If all else fails, one may try to pass one or all of the above to CMake during project initialization like this `cmake .. -DTBB_ROOT=<path to tbb dir>`.


## Cite

SpydrPick was developed as part of an academic project, and builds on the same underlying library as [SuperDCA](https://github.com/santeripuranen/SuperDCA). Please cite:

* Johan Pensar, Santeri Puranen, Neil MacAlasdair, Juri Kuronen, Gerry Tonkin-Hill, Maiju Pesonen, Brian Arnold, Yingying Xu, Aleksi Sipola, Leonor Sanchez-Buso, John A Lees, Claire Chewapreecha, Stephen D Bentley, Simon R Harris, Julian Parkhill, Nicholas J Croucher, Jukka Corander. **Genome-wide epistasis and co-selection study using mutual information.** *Nucleic Acids Research* 2019, gkz656 [**doi:** 10.1093/nar/gkz656](https://doi.org/10.1093/nar/gkz656)

* Santeri Puranen, Maiju Pesonen, Johan Pensar, Ying Ying Xu, John A. Lees, Stephen D. Bentley, Nicholas J. Croucher and Jukka Corander. **SuperDCA for genome-wide epistasis analysis.** *Microbial Genomics* 2018;4, [**doi:** 10.1099/mgen.0.000184](https://doi.org/10.1099/mgen.0.000184)

* https://github.com/santeripuranen/SpydrPick
