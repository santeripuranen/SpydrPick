/** @file aracne.cpp
    Stand-alone MI_GWES, prototype version.

    Copyright (c) 2018 Juri Kuronen and Santeri Puranen.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

    @author Juri Kuronen and Santeri Puranen
*/

#include <algorithm>
#include <cstdint>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

// Timekeeping macros before we move to apegrunt's stopwatch.
#define TIME_NOW (std::chrono::high_resolution_clock::now())
#define TIME_TAKEN(start, end) (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count())

/*
 * Neat printing format.
*/
std::string get_time(uint32_t time) {
    uint32_t hours = time / 1000 / 3600;
    uint32_t minutes = (time / 1000 / 60) % 60;
    uint32_t seconds = (time / 1000) % 60;
    return (hours ? std::to_string(hours) + " h " : "") + (minutes ? std::to_string(minutes) + " m " : "")
        + (!hours && ((minutes && seconds) || (!minutes && seconds >= 10)) ? std::to_string(seconds) + " s" : "")
        + (!hours && !minutes && seconds < 10 ? std::to_string(time) + " ms" : "");
}

/*
 * Struct for easier storage of edges.
*/
struct edge_struct {
    uint32_t var1;
    uint32_t var2;
    double mi;
    edge_struct(uint32_t v1, uint32_t v2, double m) : var1(v1), var2(v2), mi(m) {}
};

/*
 * Calculates the mutual information between variables 'var1' and 'var2'.
*/
double calculate_mi(const std::vector<uint32_t>& dummy_data, const std::vector<double>& weights, const std::vector<uint32_t>& p_indices, uint32_t var1, uint32_t var2, uint32_t n) {
    uint32_t var1_noc = p_indices[var1 + 1] - p_indices[var1];
    uint32_t var2_noc = p_indices[var2 + 1] - p_indices[var2];
    double mi = 0.0;
    std::vector<double> p;
    for (uint64_t col1 = p_indices[var1]; col1 < p_indices[var1 + 1]; ++col1) {
        for (uint64_t col2 = p_indices[var2]; col2 < p_indices[var2 + 1]; ++col2) {
            double p_tmp = 0.0;
            for (uint64_t row = 0; row < n; ++row) p_tmp += dummy_data[col1 * n + row] * weights[row] * dummy_data[col2 * n + row];
            p.push_back(p_tmp + 0.5);
        }
    }
    double sum = std::accumulate(p.begin(), p.end(), 0.0);
    for (double& x : p) x /= sum;
    
    // p12 * log(p12)
    for (double& x : p) mi += x * std::log(x);

    // p1 * log(p1)
    for (uint32_t i = 0; i < var1_noc; ++i) {
        double x = 0.0;
        for (uint32_t j = 0; j < var2_noc; ++j) x += p[i * var2_noc + j];
        mi -= x * std::log(x);
    }

    // p2 * log(p2)
    for (uint32_t j = 0; j < var2_noc; ++j) {
        double x = 0.0;
        for (uint32_t i = 0; i < var1_noc; ++i) x += p[i * var2_noc + j];
        mi -= x * std::log(x);
    }
    return mi;
}

/*
 * Calculate mutual information values for all pairs in the data.
 * Returns sorted list of edges by their MI values.
*/
std::vector<edge_struct> calculate_all_mi(const std::vector<uint32_t>& dummy_data, const std::vector<double>& weights, const std::vector<uint32_t>& p_indices,
        double threshold, uint32_t n_threads, uint32_t n, uint32_t d, uint32_t block_size, uint32_t verbose) {
    std::vector<edge_struct> edges;
    std::vector<std::vector<edge_struct>> store;
    std::vector<std::thread> threads(n_threads);
    auto lambda = [&dummy_data, &weights, &p_indices, &store, threshold, n_threads, n, d, block_size](uint32_t thr, uint32_t start, uint32_t end, uint32_t start2, uint32_t end2) {
        for (uint32_t var1 = start; var1 < end; var1 += n_threads) {
            for (uint32_t var2 = std::max<uint32_t>(var1 + 1, start2); var2 < end2; ++var2) {
                double mi = calculate_mi(dummy_data, weights, p_indices, var1, var2, n);
                if (mi > threshold) store[thr].push_back(edge_struct(var1, var2, mi));
            }
        }
    };
    auto start = TIME_NOW;
    for (uint32_t block1 = 0; block1 < d - 1; block1 += block_size) {
        auto start_block = TIME_NOW;
        if (verbose) std::printf("Calculating block (%d, %d)...", block1, std::min<uint32_t>(block1 + block_size, d - 1));
        for (uint32_t block2 = block1 + 1; block2 < d; block2 += block_size) {
            store = std::vector<std::vector<edge_struct>>(n_threads);
            for (uint32_t thr = 0; thr < n_threads; ++thr) threads[thr] = std::thread(lambda, thr, block1 + thr, std::min(block1 + block_size, d - 1), block2, std::min(block2 + block_size, d));
            for (auto& thr : threads) thr.join();
            for (uint32_t thr = 0; thr < n_threads; ++thr) {
                for (auto& e : store[thr]) edges.push_back(std::move(e));
            }
        }
        auto end_block = TIME_NOW;
        if (verbose) std::printf("\rCalculated block (%d, %d) in %s. Total time taken: %s.\n", block1, std::min<uint32_t>(block1 + block_size, d - 1), 
                get_time(TIME_TAKEN(start_block, end_block)).data(), get_time(TIME_TAKEN(start, end_block)).data());
    }
    std::sort(edges.begin(), edges.end(), [](const edge_struct& e1, const edge_struct& e2) {return e1.mi > e2.mi;});

    return edges;
}

// Can use boost/functional/hash.hpp.
class pair_hash {
public:
    std::size_t operator()(const std::pair<uint32_t, uint32_t> p) const {
        std::hash<uint64_t> hasher;
        std::size_t hash = hasher(p.first);
        hash ^= hasher(p.second) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        return hash;
    }
};

/*
 * Sample 'threshold_pairs' edge pairs.
*/
std::vector<std::pair<uint32_t, uint32_t>> sample_edges(uint32_t threshold_pairs, uint32_t d) {
    std::unordered_set<std::pair<uint32_t, uint32_t>, pair_hash> edges_map;
    std::mt19937 gen((uint32_t) TIME_NOW.time_since_epoch().count());
    std::uniform_int_distribution<uint32_t> dist(0, d - 1);
    while (edges_map.size() < threshold_pairs) {
        std::pair<uint32_t, uint32_t> edge = std::make_pair(dist(gen), dist(gen));
        while (edge.first == edge.second) edge = std::make_pair(dist(gen), dist(gen));
        if (edge.second < edge.first) std::swap(edge.first, edge.second);
        if (edges_map.count(edge) == 0) edges_map.insert(edge);
    }

    // Convert to vector.
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    for (auto& edge : edges_map) edges.push_back(std::move(edge));
    return edges;
}

/*
 * Calculate threshold for result storage.
*/
double estimate_saving_threshold(const std::vector<uint32_t>& dummy_data, const std::vector<double>& weights,
        const std::vector<uint32_t>& p_indices, uint32_t n_save, uint32_t n, uint32_t d, 
	    uint32_t threshold_iterations, uint32_t threshold_pairs, uint32_t n_threads) {
    uint32_t save_idx = (1.0 - (double) n_save / (d * (d - 1.0) / 2.0)) * threshold_pairs;
    std::vector<double> thresholds(threshold_iterations);
    auto lambda = [&dummy_data, &weights, &p_indices, &thresholds, n, d, threshold_iterations, threshold_pairs, save_idx, n_threads](uint32_t thr) {
        for (uint32_t iter = thr; iter < threshold_iterations; iter += n_threads) {
            std::vector<double> mi(threshold_pairs);
            std::vector<std::pair<uint32_t, uint32_t>> edges = sample_edges(threshold_pairs, d);
            for (uint32_t i = 0; i < threshold_pairs; ++i) mi[i] = calculate_mi(dummy_data, weights, p_indices, edges[i].first, edges[i].second, n);
            std::nth_element(mi.begin(), mi.begin() + save_idx, mi.end());
            thresholds[iter] = mi[save_idx];
        }
    };
    std::vector<std::thread> threads(n_threads);
    for (uint32_t thr = 0; thr < n_threads; ++thr) threads[thr] = std::thread(lambda, thr);
    for (auto& thr : threads) thr.join();
    std::sort(thresholds.begin(), thresholds.end());
    return (threshold_iterations % 2 ? thresholds[threshold_iterations / 2] : (thresholds[threshold_iterations / 2 - 1] + thresholds[threshold_iterations / 2]) / 2.0);
}

/*
 * Computes start/end column indices for variables in the dummy data.
*/
std::vector<uint32_t> compute_p_indices(const std::vector<uint32_t>& input, uint64_t n, uint64_t d) {
    std::vector<uint32_t> p_indices(d + 1);
    for (uint64_t col = 0; col < d; ++col) {
        std::set<uint32_t> outcomes;
        for (uint64_t row = 0; row < n; ++row) outcomes.insert(input[col * n + row]);
        p_indices[col + 1] = p_indices[col] + outcomes.size();
    }
    return p_indices;
}

/*
 * Converts input data into dummy data.
*/
std::vector<uint32_t> compute_dummy_data(const std::string& data_file, std::vector<uint32_t>& p_indices, uint64_t n, uint64_t d) {
    std::vector<uint32_t> input(n * d);
    std::ifstream(data_file, std::ios::in | std::ios::binary).read((char*) input.data(), n * d * sizeof(uint32_t));
    p_indices = compute_p_indices(input, n, d);
    std::vector<uint32_t> dummy_data(n * p_indices.back());
    for (uint64_t col = 0; col < d; ++col) {
        // Find unique outcomes in this column.
        std::set<uint32_t> outcomes;
        for (uint64_t row = 0; row < n; ++row) outcomes.insert(input[col * n + row]);
        // Map outcomes to {0, 1, ...}
        std::map<uint32_t, uint32_t> outcomes_map;
        uint32_t cnt = 0;
        for (uint32_t oc : outcomes) outcomes_map.emplace(oc, cnt++);
        // Create dummy columns.
        for (uint64_t dummy_col = p_indices[col]; dummy_col < p_indices[col + 1]; ++dummy_col) {
            uint32_t dummy_oc = dummy_col - p_indices[col];
            for (uint64_t row = 0; row < n; ++row) dummy_data[dummy_col * n + row] = outcomes_map[input[col * n + row]] == dummy_oc;
        }
    }
    return dummy_data;
}

// Todo: Use boost/program_options.
char* get_argument(char** begin, char** end, const std::string& opt, const std::string& alt) {
    auto itr = std::find(begin, end, opt);
    if (itr != end && ++itr != end) return *itr;
    itr = std::find(begin, end, alt);
    return (itr != end && ++itr != end ? *itr : nullptr);
}

/*
 * The main program.
 *
 * Compile with g++ -std=c++11 -pthread -O3 mi_gwes.cpp -o mi_gwes
 * Minimal run with ./mi_gwes -f ARG -l ARG -w ARG -d ARG -n ARG -o ARG 
 * Display help with '-h'.
 *
 * Command line arguments:
 *  -f, --data-file
 *  -l, --loci-file
 *  -w, --weights-file
 *  -d, --data-variables
 *  -n, --data-samples
 *  -o, --output-file-prefix
 * [-m, --mi-values]
 * [-i, --threshold-iterations]
 * [-p, --threshold-pairs
 * [-t, --threads]
 * [-b, --block-size]
 * [-d, --debug]
 * [-v, --verbose
 * [-h, --help]
*/
int main(int argc, char** argv) {
    // Read arguments
    if (argc == 1) {
        std::printf("mi_gwes - description goes here\n\nCopyright goes here\n\nUse '-h' or '--help' for a list of available options.\n");
        return 0;
    }
    char** argv_end = argv + argc;
    if (std::find(argv, argv_end, std::string("-h")) != argv_end || std::find(argv, argv_end, std::string("--help")) != argv_end) {
        std::printf("mi_gwes - description goes here\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n%-40s %-30s\n",
                "  -f [ --data-file ] ARG", "Path to file containing data.",
                "  -l [ --loci-file ] ARG", "Path to file containing loci identifiers.",
                "  -w [ --weights-file ] ARG", "Path to file containing weights.",
                "  -d [ --data-variables ] ARG", "Number of variables (columns) in data.",
                "  -n [ --data-samples ] ARG", "Number of data samples (rows) in data.",
                "  -o [ --output-file-prefix ] ARG", "Path (prefix) for output file.",
                " (-m [ --mi-values ] ARG)", "Approx. number of MI values to calculate from data (=100*d).",
                " (-i [ --threshold-iterations ] ARG)", "Number of iterations for estimating saving threshold.",
                " (-p [ --threshold-pairs ] ARG)", "Number of sampled pairs for estimating saving threshold.",
                " (-t [ --threads ] ARG)", "Number of threads (=1).",
                " (-b [ --block-size ] ARG)", "Block size for MI calculations (=512).",
                " (-d [ --debug ]", "Print debug information.",
                " (-v [ --verbose ]", "Be verbose.",
                " (-h [ --help ])", "Print this list.");
        return 0;
    }
    char* arg_reader = get_argument(argv, argv_end, "-f", "--data-file"); if (!arg_reader) {std::printf("Error: Missing data file (-f).\n"); return 0;}
    std::string data_file(arg_reader);

    arg_reader = get_argument(argv, argv_end, "-l", "--loci-file"); if (!arg_reader) {std::printf("Error: Missing loci identifiers file (-l).\n"); return 0;}
    std::string loci_file(arg_reader);

    arg_reader = get_argument(argv, argv_end, "-w", "--weights-file"); if (!arg_reader) {std::printf("Error: Missing weights file (-w).\n"); return 0;}
    std::string weights_file(arg_reader);

    arg_reader = get_argument(argv, argv_end, "-o", "--output-file"); if (!arg_reader) {std::printf("Error: Missing output file (-o).\n"); return 0;}
    std::string output_file(arg_reader);

    arg_reader = get_argument(argv, argv_end, "-d", "--data-variables"); if (!arg_reader) {std::printf("Error: Missing number of variables (columns) in data (-d).\n"); return 0;}
    uint64_t d(std::stoll(arg_reader));

    arg_reader = get_argument(argv, argv_end, "-n", "--data-samples"); if (!arg_reader) {std::printf("Error: Missing number of data samples (rows) in data (-n).\n"); return 0;}
    uint64_t n(std::stoll(arg_reader));

    arg_reader = get_argument(argv, argv_end, "-m", "--mi-values");
    uint64_t n_save = (arg_reader ? std::stoll(arg_reader) : 100LL * d);

    arg_reader = get_argument(argv, argv_end, "-p", "--threshold-pairs");
    uint64_t threshold_pairs = (arg_reader ? std::stoll(arg_reader) : 100000LL); // No safeguard in place at the moment.

    arg_reader = get_argument(argv, argv_end, "-i", "--threshold-iterations");
    uint32_t threshold_iterations = (arg_reader ? std::stoi(arg_reader) : 20);

    arg_reader = get_argument(argv, argv_end, "-t", "--threads");
    uint32_t n_threads = (arg_reader ? std::stoi(arg_reader) : 1);

    arg_reader = get_argument(argv, argv_end, "-b", "--block-size");
    uint32_t block_size = (arg_reader ? std::stoi(arg_reader) : 512);

    uint32_t debug = (std::find(argv, argv_end, std::string("-d")) != argv_end || std::find(argv, argv_end, std::string("--debug")) != argv_end ? 1 : 0);
    uint32_t verbose= (std::find(argv, argv_end, std::string("-v")) != argv_end || std::find(argv, argv_end, std::string("--verbose")) != argv_end ? 1 : 0);

    // Main program start
    auto start = TIME_NOW;
    // Read loci and weights files.
    std::vector<uint32_t> loci(d);
    std::ifstream(loci_file, std::ios::in | std::ios::binary).read((char*) loci.data(), d * sizeof(uint32_t));
    std::vector<double> weights(n);
    std::ifstream(weights_file, std::ios::in | std::ios::binary).read((char*) weights.data(), n * sizeof(double));
    if (debug) std::printf("(Debug) sum(weights)=%f\n", std::accumulate(weights.begin(), weights.end(), 0.0));
    std::vector<uint32_t> p_indices;
    // Read data file.
    std::vector<uint32_t> dummy_data = compute_dummy_data(data_file, p_indices, n, d);
    auto end = TIME_NOW;
    std::printf("Read and processed input files -- %s\n", get_time(TIME_TAKEN(start, end)).data());

    // Estimate saving threshold
    start = TIME_NOW;
    double threshold = estimate_saving_threshold(dummy_data, weights, p_indices, n_save, n, d, threshold_iterations, threshold_pairs, n_threads);
    end = TIME_NOW;
    std::printf("Using saving threshold th=%f -- %s\n", threshold, get_time(TIME_TAKEN(start, end)).data());

    // Calculate MI values.
    start = TIME_NOW;
    std::vector<edge_struct> edge_structs = calculate_all_mi(dummy_data, weights, p_indices, threshold, n_threads, n, d, block_size, verbose);
    end = TIME_NOW;
    std::printf("Calculated MI values -- %s\n", get_time(TIME_TAKEN(start, end)).data());

    // Write output files.
    start = TIME_NOW;
    uint32_t output_n_edges = edge_structs.size();
    std::vector<uint32_t> edges;
    for (auto& e : edge_structs) {edges.push_back(loci[e.var1]); edges.push_back(loci[e.var2]);}
    std::vector<double> mi;
    for (auto& e : edge_structs) mi.push_back(e.mi);
    std::string output_edges_file(output_file + "_" + std::to_string(output_n_edges) + "_edges");
    std::string output_mi_file(output_file + "_" + std::to_string(output_n_edges) + "_mi");
    std::ofstream(output_edges_file, std::ios::out | std::ios::binary).write((char*) edges.data(), edges.size() * sizeof(uint32_t));
    std::ofstream(output_mi_file, std::ios::out | std::ios::binary).write((char*) mi.data(), mi.size() * sizeof(double));
    end = TIME_NOW;
    std::printf("Wrote %d edges and MI values to output files %s and %s -- %s\n", output_n_edges, output_edges_file.data(), output_mi_file.data(), get_time(TIME_TAKEN(start, end)).data());
}
