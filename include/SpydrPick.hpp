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

#include "misc/Matrix_math.hpp" // apegrunt

#include "graph/Graph.hpp" // apegrunt
#include "graph/Graph_output_format.hpp" // apegrunt

#include "mi.hpp"

namespace spydrpick {

template< typename RealT, typename StateT >
apegrunt::Graph_ptr get_MI_network( std::vector< apegrunt::Alignment_ptr<StateT> >& alignments, RealT mi_threshold=0.0 )
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

	return mi_solver.get_graph();
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

template< typename RealT, typename StateT >
RealT determine_MI_threshold( apegrunt::Alignment_ptr<StateT> alignment, std::size_t n_values )
{
    using default_state_t = apegrunt::nucleic_acid_state_t;
    using alignment_default_storage_t = apegrunt::Alignment_impl_block_compressed_storage< apegrunt::StateVector_impl_block_compressed_alignment_storage<default_state_t> >;
    using node_t = std::size_t;

    const auto n_loci = alignment->n_loci();
    const std::size_t threshold_iterations = SpydrPick_options::get_mi_threshold_iterations();
    std::size_t threshold_pairs = SpydrPick_options::get_mi_threshold_pairs();

    // Calculate threshold estimate index in the empirical CDF.
    std::size_t possible_pairs = n_loci * (n_loci - 1) / 2;
    std::size_t threshold_idx = (1.0 - double(n_values) / possible_pairs) * threshold_pairs;
    std::vector<RealT> thresholds;

    // Safeguard against small alignments.
    if( possible_pairs / 10 < threshold_pairs )
    {
        threshold_pairs = possible_pairs / 10;
    }

    const auto verbose = SpydrPick_options::verbose();
    SpydrPick_options::set_verbose( false );

    for( std::size_t iter=0; iter<threshold_iterations; ++iter )
    {
    	std::vector<RealT> mi_values(threshold_pairs);
        const auto pairs = sample_pairs<node_t>( threshold_pairs, n_loci - 1 );
        #ifndef SPYDRPICK_NO_TBB
        tbb::parallel_for( tbb::blocked_range<std::size_t>(0, threshold_pairs), [alignment, &mi_values, &pairs](tbb::blocked_range<std::size_t>& r)
        {
            for( std::size_t pair_idx = r.begin(); pair_idx < r.end(); ++pair_idx)
            {
                const auto pair = pairs[pair_idx];
                auto mi_solver = get_MI_solver( std::vector< apegrunt::Alignment_ptr<StateT> >{alignment}, 0.0, SpydrPick_options::get_mi_pseudocount() );
                mi_solver( pair.first, apegrunt::get_block_index( pair.second), false );
                const auto& graph = *mi_solver.get_graph();
                mi_values[pair_idx] = graph[ pair.second % apegrunt::StateBlock_size ].weight();
            }
        });
        #else
        
        #endif // #ifndef SPYDRPICK_NO_TBB
        std::nth_element( mi_values.begin(), mi_values.begin() + threshold_idx, mi_values.end() );
        thresholds.push_back( mi_values[threshold_idx] );
    }

    SpydrPick_options::set_verbose( verbose );

    std::size_t median_idx = thresholds.size() / 2 - ( thresholds.size() % 2 ? 0 : 1 );
    auto nth = thresholds.begin() + median_idx;
    std::nth_element( thresholds.begin(), nth, thresholds.end() );
    return *nth; //thresholds[median_idx];
}

} // namespace spydrpick


#endif // SPYDRPICK_HPP
