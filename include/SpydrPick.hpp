/** @file spydrpick.hpp

	Copyright (c) 2017-2018 Santeri Puranen.

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

	@author Santeri Puranen
	$Id: $
*/
#ifndef SPYDRPICK_HPP
#define SPYDRPICK_HPP

#include <numeric> // for std::accumulate
#include <memory> // for std::shared_ptr and std::make_shared
#include <vector>
#include <set>

#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
#pragma message("Compiling with TBB support")
//#include "tbb/tbb.h"
#include "tbb/parallel_reduce.h"
//#include "tbb/parallel_for.h"
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

template< typename RealT, typename StateT >
RealT determine_MI_threshold( apegrunt::Alignment_ptr<StateT> alignment, std::size_t nvalues )
{
	using default_state_t = apegrunt::nucleic_acid_state_t;
	using alignment_default_storage_t = apegrunt::Alignment_impl_block_compressed_storage< apegrunt::StateVector_impl_block_compressed_alignment_storage<default_state_t> >;
/*
	auto include_list = apegrunt::make_Loci_list( );
	auto new_alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().include( alignment, include_list );
*/
	return 0.11;
}

} // namespace spydrpick


#endif // SPYDRPICK_HPP
