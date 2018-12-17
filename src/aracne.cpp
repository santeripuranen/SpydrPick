/** @file aracne.cpp
	Stand-alone ARACNE, prototype version.

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
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

//#include "aracne.h"
//#include "aracne.hpp"

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
 * Stores edge data, where edge_id and edge vertices nodei and nodej are encoded as
 * 0, 1, ..., and uint32_t address space should be sufficient.
 * Todo: Can improve padding.
*/
struct edge {
    uint32_t edge_id;
    uint32_t nodei;
    uint32_t nodej;
    uint32_t marked_for_removal = false;
    double mi;
    edge() {edge(0, 0, 0, 0.0);}
    edge(uint32_t eid, uint32_t ni, uint32_t nj, double m) : edge_id(eid), nodei(ni), nodej(nj), mi(m) {}
    void update(uint32_t eid, uint32_t ni, uint32_t nj, double m) {edge_id = eid; nodei = ni; nodej = nj; mi = m;}
};

/*
 * Searches for key (=node_id) in the neighborhood vector.
 * The neighborhood vector stores each node's neighbors as (edge_id, node_id) pairs.
*/
int32_t binary_search(std::vector<std::pair<uint32_t, uint32_t>>& vector, uint32_t key) {
    for (int32_t a = 0, b = vector.size() - 1, mid = (a + b) / 2; a <= b; mid = (a + b) / 2) {
        uint32_t val = vector[mid].first;
        if (val == key) return vector[mid].second;
        if (val > key) b = mid - 1; else a = mid + 1;
    }
    return -1;
}

/*
 * Finds intersection of ne(node1) and ne(node2). Nodes in this intersection form 3-cliques/triangles with node1 and node2.
 * The nodes vector stores each node's neighbors as a vector of (edge_id, node_id) pairs.
*/
std::vector<std::pair<uint32_t, uint32_t>> intersection(std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& nodes, uint32_t nodei, uint32_t nodej) {
    std::vector<std::pair<uint32_t, uint32_t>> intersection_edges;
    for (auto& neighbor : nodes[nodei]) {
        int32_t res = binary_search(nodes[nodej], neighbor.first); // Encode as -1 if not found, edge_id if found.
        if (res >= 0) intersection_edges.emplace_back(res, neighbor.second);
    }
    return intersection_edges;
}

/*
 * This custom sort function attempts to move as much of data as possible before sorting a smaller range.
*/
void custom_sort(std::vector<std::pair<uint32_t, uint32_t>>& vector, uint32_t new_elements) {
    // Find min_val and max_val in new_elements.
    auto min_val = std::min_element(vector.end() - new_elements, vector.end());
    auto max_val = std::max_element(vector.end() - new_elements, vector.end());
    // Find range of elements in (min_val, max_val) in sorted part of the vector.
    auto lower_bound = std::lower_bound(vector.begin(), vector.end() - new_elements, *min_val);
    auto upper_bound = std::upper_bound(vector.begin(), vector.end() - new_elements, *max_val);
    // Save new values.
    std::vector<std::pair<uint32_t, uint32_t>> new_values(std::distance(vector.end() - new_elements, vector.end()));
    std::move(vector.end() - new_elements, vector.end(), new_values.begin());
    // Make space for new values.
    std::move_backward(upper_bound, vector.end() - new_elements, vector.end());
    std::move(new_values.begin(), new_values.end(), upper_bound);
    // Finally sort the smaller range.
    std::sort(lower_bound, upper_bound + new_elements);
}

/*
 * If node id was already mapped, returns mapped node id. Otherwise maps node id in edge list to {0, 1, ..., 2^32} and (possibly) creates a new node mutex.
*/
inline uint32_t get_node_mapping_and_do_allocations(std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& nodes, std::unordered_map<uint64_t, uint32_t>& node_mapping, uint64_t new_node, std::vector<std::shared_ptr<std::mutex>>& node_mtx, uint32_t node_mtx_grouping_size) {
    if (node_mapping.count(new_node) == 0) { // Map new node and create new mutex
        if (nodes.size() % node_mtx_grouping_size == 0) node_mtx.emplace_back(std::shared_ptr<std::mutex>(new std::mutex));
        node_mapping.emplace(new_node, (uint32_t) nodes.size());
        nodes.emplace_back();
    }
    return node_mapping[new_node];
}

/*
 * Thread-safe emplace_back.
*/
inline void safe_emplace_back(std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& nodes, std::shared_ptr<std::mutex> node_mtx, uint32_t nodei_id, uint32_t nodej_id, uint32_t edge_id) {
    node_mtx->lock();
    nodes[nodei_id].emplace_back(nodej_id, edge_id);
    node_mtx->unlock();
}

/*
 * Joins the maps created by different threads of add_edges to a single vector.
 * The vector contains (node_id, new_edges) pairs.
*/
std::vector<std::pair<uint32_t, uint32_t>> join_maps_to_vector(std::vector<std::map<uint32_t, uint32_t>>& processed_nodes_map, uint32_t n_threads) {
    std::vector<std::pair<uint32_t, uint32_t>> processed_nodes;
    for (uint32_t thr = 1; thr < n_threads; ++thr) {
        for (auto& pn : processed_nodes_map[thr]) processed_nodes_map[0][pn.first] += pn.second; // Zero-initializes if element doesn't exist
    }
    for (auto&& pn : processed_nodes_map[0]) processed_nodes.emplace_back(std::move(pn));
    return processed_nodes;
}

/*
 * Reads next block of edges into the edge vectors, maps new nodes and creates new node mutexes.
 * This function cannot be reasonably parallelized.
*/
void read_edges(const std::vector<double>& data, std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& nodes, std::unordered_map<uint64_t, uint32_t>& node_mapping,
        std::vector<edge>& edges, std::vector<std::shared_ptr<std::mutex>>& node_mtx, uint32_t node_mtx_grouping_size, uint32_t block_start, uint32_t block_end) {
    for (uint32_t edge_id = block_start; edge_id < block_end; ++edge_id) {
        // Maps node ids in edge list to {0, 1, ..., 2^32}
        uint32_t nodei_id = get_node_mapping_and_do_allocations(nodes, node_mapping, (uint64_t) data[edge_id * 3], node_mtx, node_mtx_grouping_size);
        uint32_t nodej_id = get_node_mapping_and_do_allocations(nodes, node_mapping, (uint64_t) data[edge_id * 3 + 1], node_mtx, node_mtx_grouping_size);
        double mi = data[edge_id * 3 + 2];
        edges[edge_id].update(edge_id, nodei_id, nodej_id, mi);
    }
}

/*
 * The ARACNE procedure.
*/
void aracne(const std::string& filename, uint32_t n_edges, std::vector<std::vector<std::pair<uint32_t, uint32_t>>>& nodes, std::unordered_map<uint64_t, uint32_t>& node_mapping, 
        std::vector<edge>& edges, double threshold, uint32_t n_threads, uint32_t block_size, uint32_t node_mtx_grouping_size) {
    // Start by reading all data to memory.
    // Todo: We don't have to, we can keep the ifstream open and read one block at a time.
    std::vector<double> data(n_edges * 3);
    std::ifstream(filename, std::ios::in | std::ios::binary).read((char*) data.data(), n_edges * 3 * sizeof(double));
    // Also: We need more file reading options (csv-file) and can use Santtu's data structures.

    std::vector<std::shared_ptr<std::mutex>> node_mtx;

    /*
     * Adds next block of edges into nodes' neighborhoods and keeps track of processed node neighborhoods for faster re-sorting.
     * Note: Trusts that user provided no duplicate edges.
    */
    auto add_edges = [&nodes, &edges, &node_mtx, node_mtx_grouping_size, n_threads]
        (std::vector<std::map<uint32_t, uint32_t>>* processed_nodes_map, uint32_t block_start, uint32_t block_end, uint32_t thr) {
        // The map contains (node_id, new_edges) pairs. One map per thread.
        std::vector<std::map<uint32_t, uint32_t>>& processed_nodes_map_deref = *processed_nodes_map;
        // Until TBB implementation with proper task sharing, all threads operate on the same block.
        for (uint32_t edge_id = block_start + thr; edge_id < block_end; edge_id += n_threads) {
            uint32_t nodei_id = edges[edge_id].nodei;
            uint32_t nodej_id = edges[edge_id].nodej;
            safe_emplace_back(nodes, node_mtx[nodei_id / node_mtx_grouping_size], nodei_id, nodej_id, edge_id);
            safe_emplace_back(nodes, node_mtx[nodej_id / node_mtx_grouping_size], nodej_id, nodei_id, edge_id);
            // Direct access operator zero-initializes if element doesn't exist.
            ++processed_nodes_map_deref[thr][nodei_id];
            ++processed_nodes_map_deref[thr][nodej_id];
        }
    };

    /*
     * Re-sorts node neighborhoods after new block of edges were added.
    */
    auto sort_nodes = [&nodes](std::vector<std::pair<uint32_t, uint32_t>>* processed_nodes, uint32_t thr, uint32_t thr_work) {
        std::vector<std::pair<uint32_t, uint32_t>>& processed_nodes_deref = *processed_nodes;
        for (uint32_t i = thr * thr_work; i < (thr + 1) * thr_work && i < (uint32_t) processed_nodes_deref.size(); ++i) {
            //std::sort(nodes[processed_nodes_deref[i].first].begin(), nodes[processed_nodes_deref[i].first].end()); // For debugging. Multithreaded std::sort.
            custom_sort(nodes[processed_nodes_deref[i].first], processed_nodes_deref[i].second);
        }
    };

    /*
     * Run the ARACNE procedure on the new block of edges.
    */
    auto process_block = [&nodes, &edges, threshold, n_threads](uint32_t block_start, uint32_t block_end, uint32_t thr) {
        // Until TBB implementation with proper task sharing, all threads operate on the same block.
        for (uint32_t i = block_start + thr; i < block_end; i += n_threads) {
                uint32_t eid1 = edges[i].edge_id;
                uint32_t nodei = edges[i].nodei;
                uint32_t nodej = edges[i].nodej;
                double mi1 = edges[i].mi;

                // Find nodes in intersection of ne(node1) and ne(node2). These nodes form 3-cliques/triangles.
                std::vector<std::pair<uint32_t, uint32_t>> intersection_edges = intersection(nodes, nodei, nodej);
                for (auto& edge_pair : intersection_edges) {
                    edge& edge2 = edges[edge_pair.first];
                    edge& edge3 = edges[edge_pair.second];
                    double mi2 = edge2.mi;
                    double mi3 = edge3.mi;
                    double min12 = std::min(mi1, mi2);
                    double min23 = std::min(mi2, mi3);
                    double min123 = std::min(min12, min23);
                    if (std::max(min12, min23) - min123 >= threshold) {
                        uint32_t min_eid = (mi1 == min123 ? eid1 : (mi2 == min123 ? edge2.edge_id : edge3.edge_id));
                        edges[min_eid].marked_for_removal = true; // Thread-safe because of atomic operation.
                    }
                }

        }
    };

    uint32_t time_reading_edges = 0;
    uint32_t time_adding_edges = 0;
    uint32_t time_sorting = 0;
    uint32_t time_processing = 0;
    std::vector<std::thread> threads(n_threads);
    for (uint32_t block = 0; block < n_edges; block += block_size) {
        // Read next block of edges.
        auto start_reading_edges = TIME_NOW;
        uint32_t block_end = std::min(block + block_size, n_edges);
        read_edges(data, nodes, node_mapping, edges, node_mtx, node_mtx_grouping_size, block, block_end); // Single-threaded
        auto end_reading_edges = TIME_NOW;
        time_reading_edges += TIME_TAKEN(start_reading_edges, end_reading_edges);

        // Add next block of edges into nodes' neighborhoods.
        auto start_adding_edges = TIME_NOW;
        // Store maps of nodes whose neighborhoods need to be re-sorted. The maps contain (node_id, new_edges) pairs. One map per thread.
        std::vector<std::map<uint32_t, uint32_t>> processed_nodes_map(n_threads);
        for (uint32_t thr = 0; thr < n_threads; ++thr) threads[thr] = std::thread(add_edges, &processed_nodes_map, block, block_end, thr);
        for (auto& thr : threads) thr.join();
        std::vector<std::pair<uint32_t, uint32_t>> processed_nodes = join_maps_to_vector(processed_nodes_map, n_threads);
        auto end_adding_edges = TIME_NOW;
        time_adding_edges += TIME_TAKEN(start_adding_edges, end_adding_edges);

        // Re-sort node neighborhoods. 
        auto start_sorting = TIME_NOW;
        //for (auto& pn : processed_nodes) std::sort(nodes[pn.first].begin(), nodes[pn.first].end()); // For debugging. Single-threaded std::sort.
        //for (auto& pn : processed_nodes) custom_sort(nodes[pn.first], pn.second); // For debugging. Single-threaded custom_sort.
        uint32_t thr_work = 1 + processed_nodes.size() / n_threads; 
        for (uint32_t thr = 0; thr < n_threads; ++thr) threads[thr] = std::thread(sort_nodes, &processed_nodes, thr, thr_work);
        for (auto& thr : threads) thr.join();
        auto end_sorting = TIME_NOW;
        time_sorting += TIME_TAKEN(start_sorting, end_sorting);

        // Run the ARACNE procedure on the new block of edges.
        auto start_processing = TIME_NOW;
        for (uint32_t thr = 0; thr < n_threads; ++thr) threads[thr] = std::thread(process_block, block, block_end, thr);
        for (auto& thr : threads) thr.join();
        auto end_processing = TIME_NOW;
        time_processing += TIME_TAKEN(start_processing, end_processing);
    }

    std::cout << "Time spent reading: " << get_time(time_reading_edges + time_adding_edges + time_sorting) << " (reading edges: " << get_time(time_reading_edges)
        << ", adding edges: " << get_time(time_adding_edges) << ", sorting: " << get_time(time_sorting) << ")\nTime spent processing: " << get_time(time_processing) << '\n';
}

/*
 * The main program.
*/
int main(int argc, char** argv) {
    if (argc < 6) {std::cerr << "Please call with filename, number of edges, threshold, block_size and threads!\n"; return 0;}

    // Read input and initialize
    auto start = TIME_NOW;
    std::string filename_in(argv[1]);
    uint64_t n_edges = std::stoll(argv[2]);
    if (n_edges > std::numeric_limits<uint32_t>::max()) {std::cerr << "n_edges = " << n_edges << " > 2^32\n"; return 0;} // Safeguard against more than 2^32 edges.
    double threshold = std::stod(argv[3]);
    if (threshold <= 0.0) threshold = std::numeric_limits<double>::epsilon(); // Use machine epsilon instead of 0 threshold, otherwise algorithm cannot work meaningfully.
    uint32_t block_size = std::stoi(argv[4]);
    uint32_t n_threads = std::max(1, std::stoi(argv[5]));
    uint32_t node_mtx_grouping_size = (argc > 6 ? std::max(1, std::stoi(argv[6])) : 16); // Not a required argument. Playing around with this doesn't change much.
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> nodes; // Vector of nodes stores each node's neighbors as a vector of (edge_id, node_id) pairs.
    std::unordered_map<uint64_t, uint32_t> node_mapping; // Nodes in user-provided edge list can range in {0, ..., 2^64} but will be mapped to {0, ..., 2^32}.
    std::vector<struct edge> edges(n_edges);
    // Start ARACNE procedure.
    aracne(filename_in, n_edges, nodes, node_mapping, edges, threshold, n_threads, block_size, node_mtx_grouping_size);
    auto end = TIME_NOW;
    std::cout << "Read input file, initialized and ran ARACNE procedure -- TOTAL " << get_time(TIME_TAKEN(start, end)) << '\n';

    // Write output as a boolean list of edges remaining after the ARACNE procedure, in same order as input.
    // Todo: Output is binary now, but we probably want a csv-format.
    start = TIME_NOW;
    std::string filename_out = filename_in + "_aracne";
    std::vector<double> data_out(n_edges);
    for (uint32_t i = 0; i < n_edges; ++i) data_out[i] = !edges[i].marked_for_removal;
    std::ofstream(filename_out, std::ios::out | std::ios::binary).write((char*) data_out.data(), n_edges * sizeof(double));
    end = TIME_NOW;
    std::cout << "Wrote output file -- " << get_time(TIME_TAKEN(start, end)) << '\n';

    // Debug: Print out count of edges remaining after the ARACNE procedure.
    std::cout << "(Debug) Count of edges remaining after the ARACNE procedure: " << std::accumulate(edges.begin(), edges.end(), 0, [](uint32_t sum, const edge& e) {return sum + !e.marked_for_removal;}) << '\n';

}
