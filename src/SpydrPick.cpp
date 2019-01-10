/** @file spydrpick.cpp
	Genome-wide epistasis analysis with MI-ARACNE

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

#include <cstdlib>
#include <iostream>
#include <fstream>

#include "boost/program_options.hpp"
#include "boost/filesystem/operations.hpp" // includes boost/filesystem/path.hpp

#include "apegrunt/Apegrunt.h"
#include "apegrunt/Apegrunt_IO_misc.hpp"
#include "apegrunt/Loci_generators.hpp"
#include "apegrunt/ValueVector_parser.hpp"
#include "misc/Stopwatch.hpp"

#include "../include/SpydrPick.h"
#include "../include/SpydrPick.hpp"

/*
 * The main program
 */
int main(int argc, char **argv)
{
	namespace po = boost::program_options;
	namespace fs = boost::filesystem;

	#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
	tbb::task_scheduler_init tbb_task_scheduler(tbb::task_scheduler_init::deferred); // Threading task scheduler
	#endif // #ifndef SPYDRPICK_NO_TBB

	using namespace spydrpick;
	using std::exit;

	using real_t = double;
	using default_state_t = apegrunt::nucleic_acid_state_t;
	using alignment_default_storage_t = apegrunt::Alignment_impl_block_compressed_storage< apegrunt::StateVector_impl_block_compressed_alignment_storage<default_state_t> >;

	// Parse command line options

	// Check the command line options
	po::options_description	all_options;
	SpydrPick_options spydrpick_options( &std::cout, &std::cerr );
	spydrpick_options.AddOptions( &all_options ); // add options of the spydrpick routine

	apegrunt::Apegrunt_options apegrunt_options( &std::cout, &std::cerr );
	apegrunt_options.AddOptions( &all_options ); // add options of apegrunt

	po::variables_map options_map;

	// Catch the alignment file, even if not specified by a flag in the input
	po::positional_options_description popt;
	popt.add("alignmentfile", -1);

	std::ostringstream options_string;
	options_string << all_options;
	spydrpick_options.m_set_options_string( options_string.str() );

	try
	{
		po::store( po::command_line_parser(argc, argv).options(all_options).positional(popt).run(), options_map );
		po::notify(options_map);

		// spydrpick options
		if( !spydrpick_options.CheckOptions(&options_map) )
		{
			exit(EXIT_FAILURE);
		}

		// Apegrunt options
		apegrunt_options.set_verbose( SpydrPick_options::verbose() );
		if( !apegrunt_options.CheckOptions(&options_map) )
		{
			exit(EXIT_FAILURE);
		}
		apegrunt_options.set_threads( SpydrPick_options::threads() );

		std::cout << SpydrPick_options::s_get_version_string() << "\n"
				  << apegrunt::Apegrunt_options::s_get_version_string() << "\n\n"
				  << SpydrPick_options::s_get_copyright_notice_string() << "\n"
				  << std::endl;

		#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
		SpydrPick_options::threads() > 0 ? tbb_task_scheduler.initialize( SpydrPick_options::threads() ) : tbb_task_scheduler.initialize(); // Threading task scheduler
		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream()
				<< "SpydrPick: TBB interface version " << TBB_INTERFACE_VERSION << "\n"
				<< "SpydrPick: TBB runtime interface version " << tbb::TBB_runtime_interface_version() << "\n"
				<< "SpydrPick: TBB task scheduler is " << ( tbb_task_scheduler.is_active() ? "ACTIVE" : "INACTIVE" );
			if( tbb_task_scheduler.is_active() )
			{
				*SpydrPick_options::get_out_stream()
					<< ": using "
					<< ( SpydrPick_options::threads() > 0 ? SpydrPick_options::threads() : tbb_task_scheduler.default_num_threads() )
					<< " threads"
				;
			}
			*SpydrPick_options::get_out_stream() << "\n" << std::endl;
		}
		#endif // #ifndef SPYDRPICK_NO_TBB
	}

	catch( std::exception& e )
	{
		// probably an unknown option
		*SpydrPick_options::get_err_stream() << "SpydrPick error: " << e.what() << "\n\n";
		*SpydrPick_options::get_out_stream() << SpydrPick_options::s_get_usage_string() << std::endl;
		exit(EXIT_FAILURE);
	}
	catch(...)
	{
		*SpydrPick_options::get_err_stream() << "SpydrPick error: Exception of unknown type!\n\n";
		exit(EXIT_FAILURE);
	}

	// setup global timer
	stopwatch::stopwatch globaltimer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ); // for timing statistics
	globaltimer.start();

	apegrunt::Loci_ptr sample_list;

	stopwatch::stopwatch cputimer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ); // for timing statistics

	// get alignments
	auto alignments = apegrunt::get_alignments<default_state_t>( 2 );

	// Check if we have the compulsory input alignment
	if( alignments.empty()) { exit(EXIT_FAILURE); }

// /*
	stopwatch::stopwatch steptimer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ); // for timing statistics
	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "SpydrPick: pre-process " << alignments.size() << " alignment" << (alignments.size() > 1 ? "s" : "") << ":\n" << std::endl;
	}
	steptimer.start();
	for( auto& alignment: alignments )
	{
		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: alignment \"" << alignment->id_string() << "\":\n";
		}
		if( alignments.size() < 2 )
		{
			if( apegrunt::Apegrunt_options::has_includelist_filename() )
			{
				if( SpydrPick_options::verbose() )
				{
					*SpydrPick_options::get_out_stream() << "SpydrPick: get include list from file \"" << apegrunt::Apegrunt_options::get_includelist_filename() << "\"\n";
				}
				cputimer.start();
				auto include_list = apegrunt::parse_Loci_list( apegrunt::Apegrunt_options::get_includelist_filename(), apegrunt::Apegrunt_options::get_input_indexing_base() );
				cputimer.print_timing_stats();

				if( SpydrPick_options::verbose() )
				{
					*SpydrPick_options::get_out_stream() << "SpydrPick: include list has " << include_list->size() << " loci.\n";
					*SpydrPick_options::get_out_stream() << "SpydrPick: trim alignment based on include list \"" << include_list->id_string() << "\"\n";
				}
				cputimer.start();
				alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().include( alignments.front(), include_list );
				cputimer.stop();
				if( SpydrPick_options::verbose() ) { cputimer.print_timing_stats(); *SpydrPick_options::get_out_stream() << "\n"; }
			}

			if( apegrunt::Apegrunt_options::has_excludelist_filename() )
			{
				if( SpydrPick_options::verbose() )
				{
					*SpydrPick_options::get_out_stream() << "SpydrPick: get exclude list from file \"" << apegrunt::Apegrunt_options::get_excludelist_filename() << "\"\n";
				}
				cputimer.start();
				auto exclude_list = apegrunt::parse_Loci_list( apegrunt::Apegrunt_options::get_excludelist_filename(), apegrunt::Apegrunt_options::get_input_indexing_base() );
				cputimer.print_timing_stats();

				if( SpydrPick_options::verbose() )
				{
					*SpydrPick_options::get_out_stream() << "SpydrPick: exclude list has " << exclude_list->size() << " loci.\n";
					*SpydrPick_options::get_out_stream() << "SpydrPick: trim alignment based on exclude list \"" << exclude_list->id_string() << "\"\n";
				}
				cputimer.start();
				alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().exclude( alignments.front(), exclude_list );
				cputimer.stop();
				if( SpydrPick_options::verbose() ) { cputimer.print_timing_stats(); *SpydrPick_options::get_out_stream() << "\n"; }
			}
		}

		if( apegrunt::Apegrunt_options::filter_alignment() )
		{
	    	// apply filters as defined on the command line
			cputimer.start();
			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "SpydrPick: apply filter rules";
				SpydrPick_options::get_out_stream()->flush();
			}
			auto alignment_filter = apegrunt::Alignment_filter( apegrunt::Alignment_filter::ParameterPolicy::AQUIRE_GLOBAL );
			//auto alignment_filter = apegrunt::Alignment_filter( apegrunt::Alignment_filter::ParameterPolicy::FILTER_SNPS );
			alignment = alignment_filter.operator()<alignment_default_storage_t>( alignment );
			cputimer.stop();
			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "\n";
				alignment->statistics( SpydrPick_options::get_out_stream() );
				cputimer.print_timing_stats(); *SpydrPick_options::get_out_stream() << "\n";
			}

		}

		if( sample_list )
		{
			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "SpydrPick: trim alignment samples based on sample list \"" << sample_list->id_string() << "\"\n";
			}
			cputimer.start();
			alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().copy_selected( alignment, sample_list, sample_list->id_string() );
			cputimer.stop();
			if( SpydrPick_options::verbose() ) { cputimer.print_timing_stats(); }
 		}

		// parse from file or calculate sample weights
		apegrunt::cache_sample_weights( alignments.back() );

		// get state frequency profile and output to file
		apegrunt::output_state_frequencies( alignments.back() );

		// output the sample-sample Hamming distance matrix
		apegrunt::output_sample_distance_matrix( alignments.back() );

	}

	//auto loci_list = apegrunt::get_loci_list();

	steptimer.stop();
	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "SpydrPick: pre-processing completed\n";
		steptimer.print_timing_stats();
		*SpydrPick_options::get_out_stream() << "\n";
	}

	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "SpydrPick: run MI-ARACNE\n" << std::endl;
	}

	// run the inference
	bool spydrpick_success = spydrpick::run_SpydrPick<double>( alignments /*, loci_list */ );

	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "SpydrPick: analysis completed\n";
	}
	globaltimer.stop(); globaltimer.print_timing_stats();
	if(	!spydrpick_success )
	{
		*SpydrPick_options::get_err_stream() << "SpydrPick error: MI failed\n\n";
		exit(EXIT_FAILURE);
	}
// */
    exit(EXIT_SUCCESS);
}
