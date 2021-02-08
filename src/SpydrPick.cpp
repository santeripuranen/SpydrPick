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
#include <iomanip>
#include <fstream>

#include "boost/program_options.hpp"
#include "boost/filesystem/operations.hpp" // includes boost/filesystem/path.hpp

#include "apegrunt/Apegrunt.h"
#include "apegrunt/Apegrunt_IO_misc.hpp"
#include "apegrunt/Loci_generators.hpp"
#include "apegrunt/ValueVector_parser.hpp"
#include "graph/Graph_utility.hpp"
#include "misc/Stopwatch.hpp"

#include "SpydrPick.h"
#include "SpydrPick.hpp"
#include "ARACNE.h"

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

	aracne::ARACNE_options aracne_options( &std::cout, &std::cerr );
	aracne_options.AddOptions( &all_options ); // add options of apegrunt

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

		// SpydrPick options
		if( !spydrpick_options.CheckOptions(&options_map) )
		{
			exit(EXIT_FAILURE);
		}
		#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
		SpydrPick_options::set_threads( SpydrPick_options::threads() > 0 ? SpydrPick_options::threads() : tbb_task_scheduler.default_num_threads() );
		//SpydrPick_options::threads() > 0 ? tbb_task_scheduler.initialize( SpydrPick_options::threads() ) : tbb_task_scheduler.initialize() ); // Threading task scheduler
		#endif // #ifndef SPYDRPICK_NO_TBB

		// Apegrunt options
		apegrunt_options.set_verbose( SpydrPick_options::verbose() );
		if( !apegrunt_options.CheckOptions(&options_map) )
		{
			exit(EXIT_FAILURE);
		}
		apegrunt_options.set_threads( SpydrPick_options::threads() );

		// ARACNE options
		aracne_options.set_verbose( SpydrPick_options::verbose() );
		if( !aracne_options.CheckOptions(&options_map) )
		{
			exit(EXIT_FAILURE);
		}
		aracne_options.set_threads( SpydrPick_options::threads() );

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
					//<< ( SpydrPick_options::threads() > 0 ? SpydrPick_options::threads() : tbb_task_scheduler.default_num_threads() )
					<< SpydrPick_options::threads()
					<< " threads"
				;
			}
			*SpydrPick_options::get_out_stream() << "\n" << std::endl;
		}
		//aracne::ARACNE_options::set_threads( SpydrPick_options::threads() > 0 ? SpydrPick_options::threads() : tbb_task_scheduler.default_num_threads() );
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

	stopwatch::stopwatch cputimer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ); // for timing statistics

	// get alignments
	auto alignments = apegrunt::get_alignments<default_state_t>( 1 );

	// Check if we have the compulsory input alignment
	if( alignments.empty() ) { exit(EXIT_FAILURE); }

	// output alignments? I guess this is mostly useful for testing purposes
	if( apegrunt::Apegrunt_options::output_alignment() )
	{
		for( auto alignment: alignments )
		{
			apegrunt::output_alignment( alignment );
		}
	}

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
				alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().include( alignment, include_list );
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
				alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().exclude( alignment, exclude_list );
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
		if( apegrunt::Apegrunt_options::has_samplelist_filename() )
		{
			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "SpydrPick: sample include list from file \"" << apegrunt::Apegrunt_options::get_samplelist_filename() << "\"\n";
			}
			cputimer.start();
			auto sample_list = apegrunt::parse_Loci_list( apegrunt::Apegrunt_options::get_samplelist_filename(), apegrunt::Apegrunt_options::get_input_indexing_base() );
			cputimer.print_timing_stats();

			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "SpydrPick: include list has " << sample_list->size() << " samples.\n";
				*SpydrPick_options::get_out_stream() << "SpydrPick: trim alignment samples based on sample list \"" << sample_list->id_string() << "\"\n";
			}
			cputimer.start();
			alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().copy_selected( alignment, sample_list, sample_list->id_string() );
			cputimer.stop();
			if( SpydrPick_options::verbose() ) { cputimer.print_timing_stats(); }
 		}

		if( apegrunt::Apegrunt_options::output_filtered_alignment() )
		{
			// output alignment to file
			apegrunt::output_alignment( alignment );
		}

		// assign sample weights (parse from file or determine automatically)
		apegrunt::cache_sample_weights( alignment );

		// output sample weigths, if specified using the cmd line option
		output_sample_weights( alignment );

		// get state frequency profile and output to file
		apegrunt::output_state_frequencies( alignment );

		// automatically determine MI threshold
		if( SpydrPick_options::get_mi_threshold() < 0 )
		{
			cputimer.start();
			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "SpydrPick: determine MI save threshold";
			}
			const std::size_t top_pairs_to_save = SpydrPick_options::get_mi_values() != 0 ? SpydrPick_options::get_mi_values() : 100*alignment->n_loci();
			const auto mi_threshold = determine_MI_threshold<double>( alignment, top_pairs_to_save );
			SpydrPick_options::set_mi_threshold( mi_threshold );
			cputimer.stop()
			;if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "SpydrPick: MI save threshold = " << std::setprecision(6) << SpydrPick_options::get_mi_threshold() << " (save approx. " << top_pairs_to_save << " top pairs)\n";
				cputimer.print_timing_stats();	*SpydrPick_options::get_out_stream() << "\n";
			}
		}

		// output the sample-sample Hamming distance matrix
		apegrunt::output_sample_distance_matrix( alignment );

	}

	steptimer.stop();
	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "SpydrPick: alignment pre-processing completed\n";
		steptimer.print_timing_stats();	*SpydrPick_options::get_out_stream() << "\n";
	}

	steptimer.start();
	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "SpydrPick: evaluate MI\n"; SpydrPick_options::get_out_stream()->flush();
	}

	auto network = spydrpick::get_MI_network( alignments, SpydrPick_options::get_mi_threshold() );

	steptimer.stop();
	if( SpydrPick_options::verbose() )
	{
		steptimer.print_timing_stats();	*SpydrPick_options::get_out_stream() << "\n";
	}

	steptimer.start();
	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "SpydrPick: sort " << network.network->size() << " edges\n"; SpydrPick_options::get_out_stream()->flush();
	}

	network.network->sort();

	steptimer.stop();
	if( SpydrPick_options::verbose() )
	{
		steptimer.print_timing_stats();	*SpydrPick_options::get_out_stream() << "\n";
	}

	if( !SpydrPick_options::no_aracne() )
	{
		steptimer.start();
		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: run ARACNE\n"; SpydrPick_options::get_out_stream()->flush();
		}

		aracne::run_ARACNE( network.network );

		steptimer.stop();
		if( SpydrPick_options::verbose() )
		{
			steptimer.print_timing_stats();	*SpydrPick_options::get_out_stream() << "\n";
		}
	}

	// output final coupling scores
	{
		std::ostringstream extension;
		extension << apegrunt::Apegrunt_options::get_output_indexing_base() << "-based." << network.network->size() << "edges"; // indicate base index

		// Ensure that we always get a unique output filename
		auto couplings_file = apegrunt::get_unique_ofstream( alignments.front()->id_string()+"."+apegrunt::size_string(alignments.front())+(alignments.size() > 1 ? ".scan" : "")+".spydrpick_couplings."+extension.str() );

		if( couplings_file->stream()->is_open() && couplings_file->stream()->good() )
		{
			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "SpydrPick: write network (" << network.network->size() << " edges) to file \"" << couplings_file->name() << "\"\n";
			}
			cputimer.start();

			// produce output
			if( apegrunt::Apegrunt_options::linear_genome() )
			{
				*(couplings_file->stream()) << apegrunt::Graph_output_formatter<default_state_t,apegrunt::LinearDistance>(network.network,alignments.front());
			}
			else
			{
				*(couplings_file->stream()) << apegrunt::Graph_output_formatter<default_state_t,apegrunt::CircularDistance>(network.network,alignments.front());
			}

			cputimer.stop(); cputimer.print_timing_stats(); *SpydrPick_options::get_out_stream() << "\n";
		}
	}

	// output outliers coupling scores
	{
		std::ostringstream extension;
		extension << apegrunt::Apegrunt_options::get_output_indexing_base() << "-based." << "outliers"; // indicate base index

		// Ensure that we always get a unique output filename
		auto couplings_file = apegrunt::get_unique_ofstream( alignments.front()->id_string()+"."+apegrunt::size_string(alignments.front())+(alignments.size() > 1 ? ".scan" : "")+".spydrpick_couplings."+extension.str() );

		if( couplings_file->stream()->is_open() && couplings_file->stream()->good() )
		{
			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << "SpydrPick: write outlier network to file \"" << couplings_file->name() << "\"\n";
			}
			cputimer.start();

			// produce output
			if( apegrunt::Apegrunt_options::linear_genome() )
			{
				*(couplings_file->stream()) << Outlier_Graph_formatter<real_t,apegrunt::LinearDistance>(network,alignments.front()->get_loci_translation(),alignments.front()->n_original_positions());
			}
			else
			{
				*(couplings_file->stream()) << Outlier_Graph_formatter<real_t,apegrunt::CircularDistance>(network,alignments.front()->get_loci_translation(),alignments.front()->n_original_positions());
			}

			cputimer.stop(); cputimer.print_timing_stats(); *SpydrPick_options::get_out_stream() << "\n";
		}

		// output alignment positions that are involved in outlier edges
		if( SpydrPick_options::verbose() )
		{
			*SpydrPick_options::get_out_stream() << "SpydrPick: extract nodes involved in outlier edges:"; SpydrPick_options::get_out_stream()->flush();
		}
		cputimer.start();
		{ // extract outlier node indices
			const auto outlier_node_list = apegrunt::extract_node_indices( network.network, [=](const auto& edge){ return edge.weight() >= network.outlier_threshold; } );
			*SpydrPick_options::get_out_stream() << " found " << outlier_node_list->size() << " edges\n";
			SpydrPick_options::get_out_stream()->flush();

			if( outlier_node_list->size() < alignments[0]->n_loci() && outlier_node_list->size() != 0 )
			{
				auto outlier_node_alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().include( alignments[0], outlier_node_list );
				cputimer.stop(); cputimer.print_timing_stats(); SpydrPick_options::get_out_stream()->flush();

				if( outlier_node_alignment->n_loci() != 0 )
				{
					if( SpydrPick_options::verbose() ) { *SpydrPick_options::get_out_stream() << "\n"; }
					apegrunt::output_alignment( outlier_node_alignment );
				}
			}

			if( SpydrPick_options::verbose() )
			{
				*SpydrPick_options::get_out_stream() << std::endl;
			}
		}
	}

	if( SpydrPick_options::verbose() )
	{
		*SpydrPick_options::get_out_stream() << "SpydrPick: analysis completed\n";
		globaltimer.stop(); globaltimer.print_timing_stats();
		SpydrPick_options::get_out_stream()->flush();
	}

	exit(EXIT_SUCCESS);
}
