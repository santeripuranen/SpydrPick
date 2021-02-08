/** @file SpydrPick.hpp

	Copyright (c) 2018-2019 Santeri Puranen, 2019 Juri Kuronen

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

	@author Santeri Puranen and Juri Kuronen
	$Id: $
*/
#ifndef SPYDRPICK_HPP
#define SPYDRPICK_HPP

#include <numeric> // for std::accumulate
#include <memory> // for std::shared_ptr and std::make_shared
#include <vector>
#include <set>
#include <random>
#include <unordered_set>

#include <boost/functional/hash/hash.hpp>

#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
#pragma message("Compiling with TBB support")
//#include "tbb/tbb.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h" // should be included by parallel_for.h
//#include "tbb/mutex.h"
#endif // SPYDRPICK_NO_TBB

#include "apegrunt/Apegrunt_IO_misc.hpp"
#include "apegrunt/Apegrunt_utility.hpp"
#include "apegrunt/Alignment.h"
#include "apegrunt/StateVector_utility.hpp"
#include "apegrunt/Alignment_utility.hpp"
#include "apegrunt/Loci.h"
#include "apegrunt/Distance.hpp"

#include "misc/Matrix_math.hpp" // apegrunt

#include "graph/Graph.hpp" // apegrunt
#include "graph/Graph_output_format.hpp" // apegrunt

#include "mi.hpp"

namespace spydrpick {

template< typename RealT >
struct MI_network
{
	using real_t = RealT;

	apegrunt::Graph_ptr network;
	apegrunt::Graph_ptr network_wo_gaps;
	real_t outlier_threshold;
	real_t extreme_outlier_threshold;
};

template< typename RealT, typename DistanceT=apegrunt::LinearDistance >
struct Outlier_Graph_formatter
{
	Outlier_Graph_formatter( MI_network<RealT>& graph, const apegrunt::Loci_ptr loci_translation, std::size_t n_original_positions )
	: m_graph(graph),
	  m_translation(loci_translation),
	  m_distance(n_original_positions)
	{
		m_graph.network_wo_gaps->sort( []( const auto& a, const auto& b ) { return a.id() < b.id(); } );
	}

	inline std::size_t distance( std::size_t pos1, std::size_t pos2 ) const
	{
		return m_distance(pos1,pos2);
	}
	MI_network<RealT> m_graph;
	const apegrunt::Loci_ptr m_translation;
	DistanceT m_distance;
};

template< typename RealT, typename Distance >
static std::ostream& operator<< ( std::ostream& os, const Outlier_Graph_formatter<RealT,Distance>& ogf )
{
	const auto& index_translation = *(ogf.m_translation);
	const std::size_t base_index = apegrunt::Apegrunt_options::get_output_indexing_base();
	const auto outlier_threshold = ogf.m_graph.outlier_threshold;
	const auto extreme_outlier_threshold = ogf.m_graph.extreme_outlier_threshold;
	const auto ld_threshold = SpydrPick_options::get_ld_threshold();

	for( const auto& edge: *(ogf.m_graph.network) )
	{
		if( edge.weight() < outlier_threshold )
		{
			break;
		}
		else
		{
			const auto outlier_edge = ogf.m_graph.network_wo_gaps->find(edge); // will return either an edge_t* or nullptr
			auto weight_wo_gaps = outlier_edge != ogf.m_graph.network_wo_gaps->end() ? outlier_edge->weight() : edge.weight();

			const auto index1 = index_translation[edge.node1()]+base_index;
			const auto index2 = index_translation[edge.node2()]+base_index;

			const auto distance = ogf.distance(index1,index2);

			if( distance > ld_threshold )
			{
				os << index1 << " " << index2
						<< " " << ogf.distance(index1,index2)
						<< " " << bool(edge)
						<< std::fixed << std::setprecision(6)
						<< " " << edge.weight()
						<< " " << weight_wo_gaps
						<< " " << std::setprecision(1) << ( 1.0 - ( weight_wo_gaps / edge.weight() ) ) * 100
						<< " " << bool(edge.weight() > extreme_outlier_threshold)
						<< "\n";
			}
		}
	}
	return os;
}

template< typename RealT, typename StateT >
MI_network<RealT> get_MI_network( std::vector< apegrunt::Alignment_ptr<StateT> >& alignments, RealT mi_threshold=0.0 )
{
	using state_t = StateT;
	using apegrunt::cbegin; using apegrunt::cend;

	const auto n_loci = alignments.front()->n_loci();
	auto block_indices = alignments.front()->get_block_indices();
	auto block_range = boost::make_iterator_range( cbegin(block_indices), cend(block_indices) );

	auto mi_solver = get_MI_solver( alignments, mi_threshold, SpydrPick_options::get_mi_pseudocount() );
	#ifndef SPYDRPICK_NO_TBB
	tbb::parallel_reduce( tbb::blocked_range<decltype(block_range.begin())>( block_range.begin(), block_range.end(), 1 ), mi_solver );
	#else
	spydrpick_ftor( block_range );
	#endif // #ifndef SPYDRPICK_NO_TBB

	const auto Q1 = mi_solver.get_quartiles(). template quartile<1>();
	const auto Q3 = mi_solver.get_quartiles(). template quartile<3>();

	const auto outlier_threshold = Q3 + 1.5*(Q3-Q1);
	const auto extreme_outlier_threshold = Q3 + 3.0*(Q3-Q1);

	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << std::setprecision(6)
			<< "SpydrPick: outlier threshold=" << outlier_threshold << "\n"
			<< "SpydrPick: extreme outlier threshold=" << extreme_outlier_threshold << "\n";
	}

	MI_network<RealT> network;
	network.network = mi_solver.get_graph();
	network.network_wo_gaps = mi_solver.get_graph_wo_gaps();
	network.outlier_threshold = outlier_threshold;
	network.extreme_outlier_threshold = extreme_outlier_threshold;

	return network;
}

template< typename NodeT >
std::vector< std::pair<NodeT,NodeT> > sample_pairs( std::size_t threshold_pairs, std::size_t max_range )
{
    using pair_t = std::pair<NodeT,NodeT>;

    std::unordered_set< pair_t, boost::hash<pair_t> > pairs_set;

    // Initialize random number generator.
    std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<NodeT> dist(0, max_range);

    // Sample unique pairs.
    while( pairs_set.size() < threshold_pairs )
    {
        pair_t pair = std::make_pair(dist(gen), dist(gen));

        while( pair.first == pair.second )
        {
            pair = std::make_pair(dist(gen), dist(gen));
        }

        if( pair.first > pair.second )
        {
            std::swap(pair.first, pair.second);
        }

        pairs_set.insert(pair);
    }

    // Convert to vector.
    std::vector< pair_t > pairs_vector;
    for( auto&& p : pairs_set )
    {
        pairs_vector.push_back(std::move(p));
    }

    return pairs_vector;
}

template< typename StateT, typename NodeT, typename RealT >
struct single_edge_MI_solver
{
	using state_t = StateT;
	using node_t = NodeT;
	using real_t = RealT;
	using my_type = single_edge_MI_solver<state_t,node_t,real_t>;

	single_edge_MI_solver( apegrunt::Alignment_ptr<StateT> alignment, const std::vector< std::pair<node_t,node_t> >& edges, std::vector<real_t>& storage )
	: m_alignment(alignment),
	  m_edges(edges),
	  m_storage(storage),
	  m_mi_solver( get_MI_solver( std::vector< apegrunt::Alignment_ptr<state_t> >{ alignment }, 0.0, SpydrPick_options::get_mi_pseudocount() ) )
	{ }

	single_edge_MI_solver( const my_type& other )
	: m_alignment( other.m_alignment ),
	  m_edges( other.m_edges ),
	  m_storage( other.m_storage ),
	  m_mi_solver( get_MI_solver( std::vector< apegrunt::Alignment_ptr<state_t> >{ other.m_alignment }, 0.0, SpydrPick_options::get_mi_pseudocount() ) )
	{ }

	~single_edge_MI_solver() = default;

	inline void operator()( tbb::blocked_range<node_t>& r ) const
	{
		for( auto pair_idx = r.begin(); pair_idx < r.end(); ++pair_idx )
		{
			(*this)( pair_idx );
		}
	}

	inline void operator()( node_t pair_idx ) const
	{
			const auto& pair = m_edges[pair_idx];
			//std::cout << "pair_idx=" << pair_idx << " pair=(" << pair.first << "," << pair.second << ")" << std::endl;
			m_storage[pair_idx] = m_mi_solver.single( pair.first, pair.second );
			//const auto& graph = *mi_solver.get_graph();
			//mi_values[pair_idx] = graph[ pair.second % apegrunt::StateBlock_size ].weight();
	}

	apegrunt::Alignment_ptr<state_t> m_alignment;
	const std::vector< std::pair<node_t,node_t> >& m_edges;
	std::vector<real_t>& m_storage;
	mutable MI_solver<real_t, state_t> m_mi_solver;

};

template< typename RealT  >
inline std::size_t determine_threshold_pairs( std::size_t threshold_pairs, std::size_t possible_pairs, RealT threshold_percentile )
{
    using real_t = RealT;

    if( threshold_pairs == 0 )
    {
        // Determine automatically.
        threshold_pairs = 100000;
        std::size_t desired_threshold_idx_from_end = 100;
        std::size_t desired_max_threshold_pairs = 500000; // If above this, 'n_values' is likely too small.
        while( threshold_pairs - threshold_percentile * threshold_pairs < desired_threshold_idx_from_end && threshold_pairs < desired_max_threshold_pairs )
        {
            threshold_pairs += 10000;
        }
    } 

    // Safeguard against small alignments.
    if( possible_pairs / 10 < threshold_pairs )
    {
        threshold_pairs = possible_pairs / 10;
    }

    return threshold_pairs;

}

template< typename RealT, typename StateT >
RealT determine_MI_threshold( apegrunt::Alignment_ptr<StateT> alignment, std::size_t n_values )
{
	using real_t = RealT;
	using state_t = StateT;
	using node_t = std::size_t;

    stopwatch::stopwatch cputimer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ); // for timing statistics

    const auto n_loci = alignment->n_loci();
    const std::size_t threshold_iterations = SpydrPick_options::get_mi_threshold_iterations();

    // Calculate threshold estimate index in the empirical CDF.
    std::size_t possible_pairs = n_loci * (n_loci - 1) / 2;
    real_t threshold_percentile = (1.0 - double(n_values) / possible_pairs);
    std::size_t threshold_pairs = determine_threshold_pairs<real_t>( SpydrPick_options::get_mi_threshold_pairs(), possible_pairs, threshold_percentile );
    std::size_t threshold_idx = threshold_percentile * threshold_pairs;

    if( SpydrPick_options::verbose() )
    {
        *SpydrPick_options::get_out_stream() << " (" << threshold_pairs << " pairs * " << SpydrPick_options::get_mi_threshold_iterations() << " iterations)\n";
        SpydrPick_options::get_out_stream()->flush();
    }

    std::vector<RealT> thresholds;

    for( std::size_t iter=0; iter<threshold_iterations; ++iter )
    {
    	if( SpydrPick_options::verbose() ) { *SpydrPick_options::get_out_stream() << "\rSpydrPick: " << std::setw(2) << std::setfill(' ') << iter+1 << "/" << SpydrPick_options::get_mi_threshold_iterations() << " | sample pairs.."; SpydrPick_options::get_out_stream()->flush(); }

    	cputimer.start();
    	std::vector<RealT> mi_values(threshold_pairs);
        const auto pairs = sample_pairs<node_t>( threshold_pairs, n_loci - 1 );

        if( SpydrPick_options::verbose() ) { *SpydrPick_options::get_out_stream() << " | evaluate MI.."; SpydrPick_options::get_out_stream()->flush(); }

        auto mi_solver = single_edge_MI_solver<state_t,node_t,real_t>( alignment, pairs, mi_values );

        #ifndef SPYDRPICK_NO_TBB
        tbb::parallel_for( tbb::blocked_range<std::size_t>(0, threshold_pairs), mi_solver );
        #else
        for( const auto edge: pairs ) { mi_solver(edge); }
        #endif // #ifndef SPYDRPICK_NO_TBB

        std::nth_element( mi_values.begin(), mi_values.begin() + threshold_idx, mi_values.end() );
        thresholds.push_back( mi_values[threshold_idx] );

        if( SpydrPick_options::verbose() ) { *SpydrPick_options::get_out_stream() << " " << cputimer; SpydrPick_options::get_out_stream()->flush(); }
    }
    if( SpydrPick_options::verbose() )
    {
        *SpydrPick_options::get_out_stream() << "\n"; SpydrPick_options::get_out_stream()->flush();
    }

    std::size_t median_idx = thresholds.size() / 2 - ( thresholds.size() % 2 ? 0 : 1 );
    auto nth = thresholds.begin() + median_idx;
    std::nth_element( thresholds.begin(), nth, thresholds.end() );
    return *nth; //thresholds[median_idx];
}

} // namespace spydrpick


#endif // SPYDRPICK_HPP
