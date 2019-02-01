/** @file aracne.hpp

	Copyright (c) 2018-2019 Juri Kuronen and Santeri Puranen.

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
	$Id: $
*/
#ifndef ARACNE_HPP
#define ARACNE_HPP

#include <algorithm>
#include <cstdint>
#include <chrono>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <unordered_map>
#include <vector>

#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
#pragma message("Compiling with TBB support")
//#include "tbb/tbb.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/blocked_range.h" // should be included by parallel_for.h
//#include "tbb/mutex.h"
#endif // SPYDRPICK_NO_TBB

#include "graph/Graph.hpp"

#include "ARACNE.h"

namespace aracne {

// Initializes keys for maps for later parallel access.
template< typename NodeIdT, typename EdgeIdT >
void initialize_node_mtx_and_neighborhoods( apegrunt::Graph_ptr network, std::unordered_map< NodeIdT, std::shared_ptr<std::mutex> >& node_mtx, std::unordered_map< NodeIdT, std::vector< std::pair<NodeIdT,EdgeIdT> > >& node_neighborhoods )
{
    using node_id_t = NodeIdT;
    using edge_id_t = EdgeIdT;
    const std::uint32_t node_mtx_grouping_size = aracne::ARACNE_options::node_grouping_size();
    for( const auto& edge : *network )
    {
        node_id_t node1 = edge.node1();
        node_id_t node2 = edge.node2();
        if( node_neighborhoods.count(node1) == 0)
        {
            node_neighborhoods[node1] = std::vector< std::pair<node_id_t,edge_id_t> >(); // Default-initialize.
        }
        if( node_neighborhoods.count(node2) == 0)
        {
            node_neighborhoods[node2] = std::vector< std::pair<node_id_t,edge_id_t> >(); // Default-initialize.
        }
        if( node_mtx.count(node1 / node_mtx_grouping_size) == 0)
        {
            node_mtx.emplace( node1 / node_mtx_grouping_size, new std::mutex );
        }
        if( node_mtx.count(node2 / node_mtx_grouping_size) == 0)
        {
            node_mtx.emplace( node2 / node_mtx_grouping_size, new std::mutex );
        }

    }
}

// Condenses information about processed edges into an easily accessible map.
template< typename NodeIdT >
std::unordered_map<NodeIdT,NodeIdT> create_processed_nodes_map( const std::vector< std::pair<NodeIdT,NodeIdT> >& processed_edges )
{
    using node_id_t = NodeIdT;
    std::unordered_map<node_id_t,node_id_t> processed_nodes;
    for( auto processed_edge : processed_edges )
    {
        ++processed_nodes[processed_edge.first];
        ++processed_nodes[processed_edge.second];
    }
    return processed_nodes;
}

template< typename NodeIdT, typename EdgeIdT >
class block_reader
{
public:
    using node_id_t = NodeIdT;
    using edge_id_t = EdgeIdT;

    block_reader( apegrunt::Graph_ptr network, std::unordered_map< node_id_t, std::vector< std::pair<node_id_t,edge_id_t> > >& node_neighborhoods, std::unordered_map< node_id_t, std::shared_ptr<std::mutex> >& node_mtx, std::vector< std::pair<node_id_t,node_id_t> >& processed_edges, std::size_t block_start )
    : m_network(network),
      m_node_neighborhoods(node_neighborhoods),
      m_node_mtx(node_mtx),
      m_processed_edges(processed_edges),
      m_block_start(block_start),
      m_node_mtx_grouping_size(aracne::ARACNE_options::node_grouping_size())
    { }

    // Read next block of edges into neighborhood vector. Trusts that network contains no duplicate edges.
    inline void operator()( tbb::blocked_range<std::size_t>& r ) const
    {
        for( auto idx = r.begin(); idx < r.end(); ++idx ) { read_edge( idx ); }
    }

    // Read next block of edges into neighborhood vector. Trusts that network contains no duplicate edges. Non-TBB implementation.
    inline void operator()( std::size_t block_start, std::size_t block_end ) const
    {
        for( auto idx = block_start; idx < block_end; ++idx ) { read_edge( idx ); }
    }

private:
    apegrunt::Graph_ptr m_network;
    std::unordered_map< node_id_t, std::shared_ptr<std::mutex> >& m_node_mtx;
    std::unordered_map< node_id_t, std::vector< std::pair<node_id_t,edge_id_t> > >& m_node_neighborhoods;
    std::vector< std::pair<node_id_t,node_id_t> >& m_processed_edges;
    std::size_t m_block_start;
    uint32_t m_node_mtx_grouping_size;

    // Read a single edge.
    inline void read_edge( std::size_t idx ) const
    {
        auto edge = (*m_network)[idx];
        edge_id_t edge_idx = idx;
        node_id_t node1 = edge.node1();
        node_id_t node2 = edge.node2();
        safe_emplace_back( node1, node2, edge_idx, m_node_mtx[node1 / m_node_mtx_grouping_size] );
        safe_emplace_back( node2, node1, edge_idx, m_node_mtx[node2 / m_node_mtx_grouping_size] );
        m_processed_edges[idx - m_block_start] = std::make_pair( node1, node2 );
    }

    // Thread-safe emplace_back.
    inline void safe_emplace_back( node_id_t node, node_id_t new_neighbor, edge_id_t edge_idx, std::shared_ptr<std::mutex> node_mtx ) const
    {
        node_mtx->lock();
        m_node_neighborhoods[node].emplace_back(new_neighbor, edge_idx);
        node_mtx->unlock();
    }

};

template< typename NodeIdT, typename EdgeIdT >
class block_sorter
{
public:
    using node_id_t = NodeIdT;
    using edge_id_t = EdgeIdT;
    using const_itr_t = typename std::unordered_map<node_id_t,node_id_t>::const_iterator;

    block_sorter( std::unordered_map< node_id_t, std::vector< std::pair<node_id_t,edge_id_t> > >& node_neighborhoods, const std::unordered_map<node_id_t,node_id_t>& processed_nodes )
    : m_node_neighborhoods(node_neighborhoods),
      m_processed_nodes(processed_nodes)
    { }

    // Re-sorts node neighborhood after new block of edges was added.
    inline void operator()( std::pair<node_id_t,node_id_t> processed_node ) const
    {
        node_id_t node = processed_node.first;
        node_id_t new_elements = processed_node.second;
        custom_sort( m_node_neighborhoods[node], new_elements );
    }

    inline const_itr_t begin() const { return std::begin(m_processed_nodes); }
    inline const_itr_t end() const { return std::end(m_processed_nodes); }

private:
    std::unordered_map< node_id_t, std::vector< std::pair<node_id_t,edge_id_t> > >& m_node_neighborhoods;
    const std::unordered_map<node_id_t,node_id_t>& m_processed_nodes;

    // Speed up sort by sorting only the relevant range in vector.
    inline void custom_sort( std::vector< std::pair<node_id_t,edge_id_t> >& vector, uint32_t new_elements ) const
    {
        // Find min_val and max_val in new elements.
        auto min_val = std::min_element( vector.end() - new_elements, vector.end() );
        auto max_val = std::max_element( vector.end() - new_elements, vector.end() );
        // Find relevant range in sorted part of the vector.
        auto lower_bound = std::lower_bound( vector.begin(), vector.end() - new_elements, *min_val );
        auto upper_bound = std::upper_bound( vector.begin(), vector.end() - new_elements, *max_val );
        // Save new values.
        std::vector< std::pair<node_id_t,edge_id_t> > new_values( vector.end() - new_elements, vector.end() );
        // Make space for new values.
        std::move_backward( upper_bound, vector.end() - new_elements, vector.end() );
        // Move new values to the end of the range.
        std::move( new_values.begin(), new_values.end(), upper_bound );
        // Finally sort the smaller range.
        std::sort( lower_bound, upper_bound + new_elements );
    }

};

template< typename NodeIdT, typename EdgeIdT, typename RealT >
class block_processor
{
public:
    using node_id_t = NodeIdT;
    using edge_id_t = EdgeIdT;
    using real_t = RealT;

    block_processor( apegrunt::Graph_ptr network, std::unordered_map< node_id_t, std::vector< std::pair<node_id_t,edge_id_t> > >& node_neighborhoods, double threshold )
    : m_network(network),
      m_node_neighborhoods(node_neighborhoods),
      m_threshold(threshold)
    { }

    // Run the ARACNE procedure on the new block of edges.
    inline void operator()( tbb::blocked_range<std::size_t>& r) const
    {
        for( auto edge_idx = r.begin(); edge_idx < r.end(); ++edge_idx ) { process_edge( edge_idx ); }
    }

    // Run the ARACNE procedure on the new block of edges. Non-TBB implementation.
    inline void operator()( std::size_t block_start, std::size_t block_end ) const
    {
        for( auto edge_idx = block_start; edge_idx < block_end; ++edge_idx ) { process_edge (edge_idx ); }
    }

private:
    const double m_threshold;
    apegrunt::Graph_ptr m_network;
    std::unordered_map< node_id_t, std::vector< std::pair<node_id_t,edge_id_t> > >& m_node_neighborhoods;

    // Process a single edge.
    inline void process_edge( std::size_t edge_idx ) const
    {
        auto edge = (*m_network)[edge_idx];
        node_id_t node1 = edge.node1();
        node_id_t node2 = edge.node2();
        real_t mi1 = edge.weight();

        std::vector< std::pair<edge_id_t,edge_id_t> > intersection_edges = intersection( m_node_neighborhoods, node1, node2 );
        for( auto& edge_pair : intersection_edges )
        {
            edge_id_t edge_idx2 = edge_pair.first;
            edge_id_t edge_idx3 = edge_pair.second;
            auto& edge2 = (*m_network)[edge_idx2];
            auto& edge3 = (*m_network)[edge_idx3];
            real_t mi2 = edge2.weight();
            real_t mi3 = edge3.weight();
            real_t min12 = std::min(mi1, mi2);
            real_t min23 = std::min(mi2, mi3);
            real_t min123 = std::min(min12, min23);
            if( std::max(min12, min23) - min123 >= m_threshold )
            {
                edge_id_t min_edge_idx = mi1 == min123 ? edge_idx : (mi2 == min123 ? edge_idx2 : edge_idx3);
                (*m_network)[min_edge_idx].set(); // Thread-safe becaues of atomic operation.
            }
        }
    }

    // Searches for key (=node) in the sorted neighborhood vector. Returns edge index if found, std::numeric_limits::max() otherwise.
    inline edge_id_t binary_search( const std::vector< std::pair<node_id_t,edge_id_t> >& vector, node_id_t key ) const
    {
        for( std::int32_t a = 0, b = vector.size() - 1, mid = (a + b) / 2; a <= b; mid = (a + b) / 2 )
        {
            node_id_t neighbor = vector[mid].first;
            if( neighbor == key )
            {
                return vector[mid].second;
            }
            if( neighbor > key )
            {
                b = mid - 1;
            }
            else
            {
                a = mid + 1;
            }
        }
        return std::numeric_limits<edge_id_t>::max();
    }

    /*
     * Finds intersection of ne(node1) and ne(node2). Nodes in this intersection form 3-cliques/triangles with node1 and node2.
     * Returns a vector of edge index pairs for edges (node1, node3) and (node2, node3) for each node3 in intersection.
    */
    inline std::vector< std::pair<edge_id_t,edge_id_t> > intersection( const std::unordered_map< node_id_t, std::vector< std::pair<node_id_t,edge_id_t> > >& node_neighborhoods, node_id_t node1, node_id_t node2 ) const
    {
        std::vector< std::pair<edge_id_t,edge_id_t> > intersection_edges;
        for( const auto& neighbor : node_neighborhoods.at(node1) )
        {
            edge_id_t res = binary_search( node_neighborhoods.at(node2), neighbor.first ); // Returns std::numeric_limits::max() if neighbor was not found.
            if( res < std::numeric_limits<edge_id_t>::max() )
            {
                intersection_edges.emplace_back( res, neighbor.second );
            };
        }
        return intersection_edges;
    }

};

// The ARACNE procedure.
//template< typename NodeIdT, typename EdgeIdT, typename RealT >
void aracne( const apegrunt::Graph_ptr input_graph )
{
    using node_id_t = uint32_t /*NodeIdT*/;
    using edge_id_t = uint64_t /*EdgeIdT*/;
    using real_t = double /*RealT*/;
	const uint32_t n_edges = input_graph->size();
	const double threshold = aracne::ARACNE_options::edge_threshold();
    const uint32_t block_size = aracne::ARACNE_options::block_size();
	bool debug = true; // Make into a (hidden?) command line argument or remove entirely.
    std::size_t verbose_output_interval = n_edges / 10;
    std::size_t verbose_previous_block_start = 0;

    // Setup timers.
    stopwatch::stopwatch steptimer( ARACNE_options::verbose() ? ARACNE_options::get_out_stream() : nullptr );
    stopwatch::stopwatch edgestimer( debug ? ARACNE_options::get_out_stream() : nullptr );
    stopwatch::stopwatch sorttimer( debug ? ARACNE_options::get_out_stream() : nullptr );
    stopwatch::stopwatch processtimer( debug ? ARACNE_options::get_out_stream() : nullptr );

    // Debug counters.
    uint64_t time_edges = 0;
    uint64_t time_sort = 0;
    uint64_t time_process = 0;
    uint64_t time_edges_global = 0;
    uint64_t time_sort_global = 0;
    uint64_t time_process_global = 0;

    std::unordered_map< node_id_t, std::shared_ptr<std::mutex> > node_mtx;
    std::unordered_map< node_id_t, std::vector< std::pair<node_id_t,edge_id_t> > > node_neighborhoods;

    steptimer.start();
    initialize_node_mtx_and_neighborhoods<node_id_t,edge_id_t>( input_graph, node_mtx, node_neighborhoods );
    steptimer.stop();
    if( aracne::ARACNE_options::verbose() ) { *aracne::ARACNE_options::get_out_stream() << "  initialization routine time=" << stopwatch::time_string(steptimer.elapsed_time()) << "\n"; }

    steptimer.start();
    for( std::size_t block_start = 0; block_start < n_edges; block_start += block_size )
    {
        std::size_t block_end = std::min( block_start + block_size, (std::size_t) n_edges );

        edgestimer.start(); // Debug timer.
        std::vector< std::pair<node_id_t,node_id_t> > processed_edges(block_end - block_start);
        // Read next block of edges into neighborhood vector. Trusts that network contains no duplicate edges.
        auto reader = block_reader<node_id_t,edge_id_t>( input_graph, node_neighborhoods, node_mtx, processed_edges, block_start );
        #ifndef SPYDRPICK_NO_TBB
        tbb::parallel_for( tbb::blocked_range<std::size_t>( block_start, block_end ), reader );
        #else
        reader()( block_start, block_end ); // Single-threaded.
        #endif // #ifndef SPYDRPICK_NO_TBB
        edgestimer.stop();
        time_edges += edgestimer.elapsed_time();

        sorttimer.start(); // Debug timer.
        // Re-sorts node neighborhoods after new block of edges was added.
        std::unordered_map<node_id_t,node_id_t> processed_nodes = create_processed_nodes_map<node_id_t>( processed_edges );
        auto sorter = block_sorter<node_id_t,edge_id_t>( node_neighborhoods, processed_nodes );
        #ifndef SPYDRPICK_NO_TBB
        tbb::parallel_for_each( sorter.begin(), sorter.end(), sorter );
        #else
        for( auto itr = sorter.begin(); itr < sorter.end(); ++itr ) { sorter()( *itr ); } // Single-threaded;
        #endif // #ifndef SPYDRPICK_NO_TBB
        sorttimer.stop();
        time_sort += sorttimer.elapsed_time();

        processtimer.start(); // Debug timer.
        auto processor = block_processor<node_id_t,edge_id_t,real_t>( input_graph, node_neighborhoods, threshold );

        // Safeguard for a particular threshold=0 case, where sequential edges with the same mi-value form a clique
        // that would not be processed correctly due to the sequence being cut by the block ending.
        if( threshold == 0.0 && block_start > 0 && (*input_graph)[block_start - 1].weight() == (*input_graph)[block_start].weight() )
        {
            // Move back to the start of the sequence and reprocess.
            std::size_t idx = block_start - 2;
            while( idx >= 0 && (*input_graph)[idx].weight() == (*input_graph)[block_start].weight() ) {--idx;}
            #ifndef SPYDRPICK_NO_TBB
            tbb::parallel_for( tbb::blocked_range<std::size_t>( idx, block_start ), processor );
            #else
            processor()( idx, block_start ); // Single-threaded.
            #endif
        }

        // Run the ARACNE procedure on the new block of edges.
        #ifndef SPYDRPICK_NO_TBB
        tbb::parallel_for( tbb::blocked_range<std::size_t>( block_start, block_end ), processor );
        #else
        processor()( block_start, block_end ); // Single-threaded.
        #endif
        processtimer.stop();
        time_process += processtimer.elapsed_time();

        if( aracne::ARACNE_options::verbose() && (block_start > verbose_previous_block_start + verbose_output_interval || block_end == n_edges ) )
        {
            std::ostringstream oss;
            std::size_t verbose_previous_block_end = verbose_previous_block_start + block_size;
            steptimer.stop();
            oss << "  processed blocks (" << verbose_previous_block_start << ", " << verbose_previous_block_end << "), ..., ("
                << block_start << ", " << block_end << ") time=" << stopwatch::time_string(steptimer.elapsed_time()) << "\n";
            if( debug )
            {
                oss << "  (Debug) reading blocks " << stopwatch::time_string(time_edges + time_sort)
                    << " (processing " << stopwatch::time_string(time_edges) 
                    << " + sorting " << stopwatch::time_string(time_sort) 
                    << "), processing blocks " << stopwatch::time_string(time_process) << '\n';
            }
            *aracne::ARACNE_options::get_out_stream() << oss.str();
            verbose_previous_block_start = block_start;
            steptimer.start();
            time_edges_global += time_edges;
            time_sort_global += time_sort;
            time_process_global += time_process;
            time_edges = 0;
            time_sort = 0;
            time_process = 0;
        }

    }

    if( debug ){
        *aracne::ARACNE_options::get_out_stream() << "  (DEBUG) TOTAL reading blocks " << stopwatch::time_string(time_edges_global + time_sort_global)
            << " (processing " << stopwatch::time_string(time_edges_global) 
            << " + sorting " << stopwatch::time_string(time_sort_global) 
            << "), processing blocks " << stopwatch::time_string(time_process_global) << '\n';
    }
}

void run_ARACNE( apegrunt::Graph_ptr network )
{
	const auto n_edges = network->size();

    // Start ARACNE procedure.
    aracne(network);

    // Write output as a boolean list of edges remaining after the ARACNE procedure, in same order as input, stored as uint32_t.
    // Todo: Input&output files are binary now, but we probably want a csv-format.
    std::vector<uint32_t> data_out(n_edges);

    auto& ntwrk = *network;
    for( std::size_t i = 0; i < n_edges; ++i)
    {
    	data_out[i] = !bool(ntwrk[i]);
    }

	// Ensure that we always get a unique output filename
	auto aracne_outfile = apegrunt::get_unique_ofstream( aracne::ARACNE_options::outfilename() );

	if( aracne_outfile->stream()->is_open() && aracne_outfile->stream()->good() )
	{
		aracne_outfile->stream()->write((char*) data_out.data(), n_edges * sizeof(uint32_t));
	}

	return;
}

} // namespace aracne

#endif // ARACNE_HPP
