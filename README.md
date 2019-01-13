# SpydrPick

## Get the code
```
git clone --recursive https://github.com/santeripuranen/SpydrPick.git
```
See how to compile [here](README.md/#building-SpydrPick).


## About

Spydrpick is a tool for performing direct coupling analysis of aligned categorical datasets. It constructs a coupling graph based on thresholded node-pair (edge) scoring,
followed by pruning of edges that represent indirect couplings.


## Building SpydrPick

Run `make` in the main directory to compile the prototype binaries.


## Running SpydrPick

Here are simple running instructions for the prototype version. You need a column-major binary data file that contains 32 bit integers, a binary weights file that contains doubles and a binary loci file that contains 32 bit integers.

First, run
```
./build/mi_gwes -f path_to_data -l path_to_loci -w path_to_weights -d number_of_variables -n number_of_samples -t number_of_threads -o output_file_prefix
```
to get two output files `output_file_prefix_count_edges` and `output_file_prefix_count_mi`. Then, run
```
./build/aracne -e output_file_prefix_count_edges -m output_file_prefix_count_mi -n count -t number_of_threads -o aracne_output_file
```
to get the output file `aracne_output_file` that contains a boolean list of direct/indirect edges.

You can then use plot.r with these three output files to draw a plot of the results. For now, user has to provide linkage disequilibrium distance `ld_dist` and genome length `genome_length`. The calculation for the two outlier thresholds `outlier_threshold` and `extreme_outlier_threshold` will be added soon.
