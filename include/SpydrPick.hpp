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
#include "apegrunt/Alignment_StateVector_weights.hpp"
#include "apegrunt/Loci.h"

#include "misc/Stopwatch.hpp" // apegrunt
#include "misc/Matrix_math.hpp" // apegrunt

#include "graph/Graph.hpp" // apegrunt
#include "graph/Graph_output_format.hpp" // apegrunt

#include "mi.hpp"

namespace spydrpick {

template< typename RealT, typename StateT >
bool run_SpydrPick( std::vector< apegrunt::Alignment_ptr<StateT> >& alignments /*, apegrunt::Loci_ptr loci_list*/ )
{
	using real_t = RealT;
	using state_t = StateT;
	using apegrunt::cbegin; using apegrunt::cend;

	if( alignments.size() == 0 )
	{
		*SpydrPick_options::get_err_stream() << "spydrpick error: no input alignment(s)\n";
		return false;
	}

	stopwatch::stopwatch cputimer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ); // for timing statistics
	stopwatch::stopwatch spydrpick_timer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ); // for timing statistics
	spydrpick_timer.start();

	if( SpydrPick_options::verbose() )
	{
		for( auto alignment: alignments )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: input alignment \"" << alignment->id_string() << "\" has " << alignment->size() << " sequences and " << alignment->n_loci() << " loci\n";
		}
		*SpydrPick_options::get_out_stream() << std::endl;
	}

	const auto n_loci = alignments.front()->n_loci();

	if( 0 == n_loci )
	{
		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: there are no loci to analyze (need at least 2)\n";
		}
		return false;
	}

	auto block_indices = alignments.front()->get_block_indices();
	auto block_range = boost::make_iterator_range( cbegin(block_indices), cend(block_indices) );

    {
		cputimer.start();

		auto network = apegrunt::make_Graph_ptr<apegrunt::Graph>(); // empty network, will store the output

		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: storage pool capacity is " << network->capacity() << " units\n";
		}
		cputimer.stop();
		if( SpydrPick_options::verbose() ) { cputimer.print_timing_stats(); *SpydrPick_options::get_out_stream() << "\n"; }

		// Evaluate pair-wise scores

		cputimer.start();

		// refresh block accounting
		//std::cout << "SpydrPick: alignments.size()=" << alignments.size() << std::endl;
/*
		for( auto& alignment: alignments )
		{
			std::cout << "SpydrPick: block accounting" << std::endl;
			alignment->get_block_accounting();
			std::cout << "SpydrPick: statepresence accounting" << std::endl;
			alignment->get_statepresence_blocks();
			std::cout << "SpydrPick: statecount accounting" << std::endl;
			alignment->get_statecount_blocks()->size();
			std::cout << "SpydrPick: refreshed" << std::endl;
		}
*/
		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: get MI solver\n"; SpydrPick_options::get_out_stream()->flush();
		}

		// The scoring stage -- this is where the magic happens
		auto spydrpick_ftor = get_MI_solver( alignments, network, SpydrPick_options::get_mi_threshold(), SpydrPick_options::get_mi_pseudocount() );
		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: evaluate MI\n"; SpydrPick_options::get_out_stream()->flush();
		}
	#ifndef SPYDRPICK_NO_TBB
		tbb::parallel_reduce( tbb::blocked_range<decltype(block_range.begin())>( block_range.begin(), block_range.end(), 1 ), spydrpick_ftor );
	#else
		spydrpick_ftor( block_range );
	#endif // #ifndef SPYDRPICK_NO_TBB
		cputimer.stop(); cputimer.print_timing_stats();

		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: network contains " << network->size() << " edges\n"; SpydrPick_options::get_out_stream()->flush();
		}

		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: sort edges\n"; SpydrPick_options::get_out_stream()->flush();
		}
		cputimer.start();
		network->sort();
		cputimer.stop(); cputimer.print_timing_stats();

	// /*
		// output final coupling scores
		std::ostringstream extension;
		extension << apegrunt::Apegrunt_options::get_output_indexing_base() << "-based"; // indicate base index

		{

			// Ensure that we always get a unique output filename
			auto couplings_file = apegrunt::get_unique_ofstream( alignments.front()->id_string()+"."+apegrunt::size_string(alignments.front())+(alignments.size() > 1 ? ".scan" : "")+".spydrpick_couplings."+extension.str()+".all" );

			if( couplings_file->stream()->is_open() && couplings_file->stream()->good() )
			{
				if( SpydrPick_options::verbose() )
				{
					*SpydrPick_options::get_out_stream() << "\nSpydrPick: writing coupling values (" << network->size() << ") to file \"" << couplings_file->name() << "\"\n";
				}
				cputimer.start();

#pragma message("Index translation not implemented!")
				*(couplings_file->stream()) << apegrunt::Graph_output_formatter<state_t>(network,alignments.front());

				// ignore index translations
				cputimer.stop(); cputimer.print_timing_stats();
			}
		}
    }

	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "\nSpydrPick: analysis completed\n";
	}
	spydrpick_timer.stop(); spydrpick_timer.print_timing_stats();

	return true;
}

} // namespace spydrpick


#endif // SPYDRPICK_HPP
