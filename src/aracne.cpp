/** @file aracne.cpp
	Stand-alone ARACNE

	Copyright (c) 2018 Santeri Puranen.

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
*/

#include "boost/program_options.hpp"
#include "boost/filesystem/operations.hpp" // includes boost/filesystem/path.hpp

#include "aracne.hpp"

/*
 * The main program
 */
int main(int argc, char **argv)
{
	namespace po = boost::program_options;
	namespace fs = boost::filesystem;

	#ifndef ARACNE_NO_TBB // Threading with Threading Building Blocks
	tbb::task_scheduler_init tbb_task_scheduler(tbb::task_scheduler_init::deferred); // Threading task scheduler
	#endif // #ifndef ARACNE_NO_TBB

	using std::exit;

	std::cout << SpydrPick_options::s_get_version_string() << "\n"
			  << apegrunt::Apegrunt_options::s_get_version_string() << "\n\n"
			  << SpydrPick_options::s_get_copyright_notice_string() << "\n"
			  << std::endl;

	// Parse command line options

	// Check the command line options
	po::options_description	all_options;
	ARACNE_options aracne_options( &std::cout, &std::cerr );
	aracne_options.AddOptions( &all_options ); // add options of the ARACNE routine

	po::variables_map options_map;

	// Catch the edgelist file, even if not specified by a flag in the input
	po::positional_options_description popt;
	popt.add("edgelistfile", -1);

	std::ostringstream options_string;
	options_string << all_options;
	aracne_options.m_set_options_string( options_string.str() );

	try
	{
		po::store( po::command_line_parser(argc, argv).options(all_options).positional(popt).run(), options_map );
		po::notify(options_map);

		// aracne options
		if( !aracne_options.CheckOptions(&options_map) )
		{
			exit(EXIT_FAILURE);
		}

		#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
		ARACNE_options::threads() > 0 ? tbb_task_scheduler.initialize( ARACNE_options::threads() ) : tbb_task_scheduler.initialize(); // Threading task scheduler
		if( ARACNE_options::verbose() )
		{
			*ARACNE_options::get_out_stream()
				<< "ARACNE: TBB interface version " << TBB_INTERFACE_VERSION << "\n"
				<< "ARACNE: TBB runtime interface version " << tbb::TBB_runtime_interface_version() << "\n"
				<< "ARACNE: TBB task scheduler is " << ( tbb_task_scheduler.is_active() ? "ACTIVE" : "INACTIVE" );
			if( tbb_task_scheduler.is_active() )
			{
				*ARACNE_options::get_out_stream()
					<< ": using "
					<< ( ARACNE_options::threads() > 0 ? ARACNE_options::threads() : tbb_task_scheduler.default_num_threads() )
					<< " threads"
				;
			}
			*ARACNE_options::get_out_stream() << "\n" << std::endl;
		}
		#endif // #ifndef SPYDRPICK_NO_TBB
	}

	catch( std::exception& e )
	{
		// probably an unknown option
		*ARACNE_options::get_err_stream() << "ARACNE error: " << e.what() << "\n\n";
		*ARACNE_options::get_out_stream() << ARACNE_options::s_get_usage_string() << std::endl;
		exit(EXIT_FAILURE);
	}
	catch(...)
	{
		*ARACNE_options::get_err_stream() << "ARACNE error: Exception of unknown type!\n\n";
		exit(EXIT_FAILURE);
	}

	// setup global timer
	stopwatch::stopwatch globaltimer( ARACNE_options::verbose() ? ARACNE_options::get_out_stream() : nullptr ); // for timing statistics
	globaltimer.start();

	stopwatch::stopwatch cputimer( ARACNE_options::verbose() ? ARACNE_options::get_out_stream() : nullptr ); // for timing statistics

#pragma warning("Kesken")

	// get network
	auto network = get_network();

	// Check if we have the compulsory input alignment
	if( network.empty()) { exit(EXIT_FAILURE); }

    exit(EXIT_SUCCESS);
}
