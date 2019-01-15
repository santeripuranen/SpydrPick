/** @file SpydrPick_options.cpp

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

#include "../include/SpydrPick_version.h"
#include "SpydrPick_options.h"

namespace spydrpick {

std::ostream* SpydrPick_options::s_out = nullptr;
std::ostream* SpydrPick_options::s_err = nullptr;
bool SpydrPick_options::s_verbose = false;

#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
int SpydrPick_options::s_threads = -1;
#else
int SpydrPick_options::s_threads = 1;
#endif // SPYDRPICK_NO_TBB

double SpydrPick_options::s_mi_threshold = -1.;
double SpydrPick_options::s_mi_pseudocount = 0.5;

std::string SpydrPick_options::s_options_string;

const std::string SpydrPick_options::s_title_string(
	  std::string("SpydrPick: Genome-wide epistasis analysis with MI-ARACNE.\n")
);

const std::string SpydrPick_options::s_usage_string(
	  std::string("Usage: SpydrPick") /*+ std::string(argv[0])*/ + " [options] <alignmentfile> [-o <outputfile>]\nOption '--help' will print a list of available options.\n"
);

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

const std::string SpydrPick_options::s_version_string(
	std::string("SpydrPick version ") + std::to_string(spydrpick_version::s_MajorVersion) + "." + std::to_string(spydrpick_version::s_MinorVersion) + "." + std::to_string(spydrpick_version::s_SubminorVersion)
	+ " | revision " + TOSTRING(SPYDRPICK_GIT_BRANCH) + "-" + TOSTRING(SPYDRPICK_GIT_COMMIT_HASH) + " | " +
#ifdef __AVX2__
	"AVX2"
#elif __AVX__
	"AVX"
#elif __SSE2__
	"SSE2"
#else
	"generic"
#endif
	+ " build " + std::string(__DATE__) + " " + std::string(__TIME__)
);

const std::string SpydrPick_options::s_copyright_notice(
	std::string("Copyright (c) 2018-2019 Santeri Puranen\nLicensed under the GNU Affero General Public License version 3.\n\n")
	+ "THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND."
);

const std::string SpydrPick_options::s_long_copyright_notice(
	std::string("Copyright (c) 2018-2019 Santeri Puranen\nLicensed under the GNU Affero General Public License version 3.\n\n")
	+ "This program is free software: you can redistribute it and/or modify\n"
    + "it under the terms of the GNU Affero General Public License as\n"
    + "published by the Free Software Foundation, either version 3 of the\n"
	+ "License, or (at your option) any later version.\n\n"
    + "This program is distributed in the hope that it will be useful,\n"
    + "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
    + "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
    + "GNU Affero General Public License for more details.\n\n"
	+ "You should have received a copy of the GNU Affero General Public License\n"
    + "along with this program. If not, see <http://www.gnu.org/licenses/>."
);

SpydrPick_options::SpydrPick_options() { this->m_init(); }
SpydrPick_options::~SpydrPick_options() { }

SpydrPick_options::SpydrPick_options( std::ostream *out, std::ostream *err )
{
	SpydrPick_options::s_out = out;
	SpydrPick_options::s_err = err;
	this->m_init();
}

const std::string& SpydrPick_options::s_get_copyright_notice_string() { return s_copyright_notice; }
const std::string& SpydrPick_options::s_get_usage_string() { return s_usage_string; }
const std::string& SpydrPick_options::s_get_version_string() { return s_version_string; }
const std::string& SpydrPick_options::s_get_title_string() { return s_title_string; }

void SpydrPick_options::set_out_stream( std::ostream *out ) { s_out = out->good() ? out : nullptr; }
void SpydrPick_options::set_err_stream( std::ostream *err ) { s_err = err->good() ? err : nullptr; }

std::ostream* SpydrPick_options::get_out_stream() { return s_out; }
std::ostream* SpydrPick_options::get_err_stream() { return s_err; }

const std::string& SpydrPick_options::s_get_options_string() { return s_options_string; }
void SpydrPick_options::m_set_options_string( const std::string& options_string ) { s_options_string = options_string; }

bool SpydrPick_options::verbose() { return ( s_verbose && s_out ); } // be verbose only if we have a valid outstream.

int SpydrPick_options::threads() { return s_threads; }

double SpydrPick_options::get_mi_threshold() { return s_mi_threshold; }
double SpydrPick_options::get_mi_pseudocount() { return s_mi_pseudocount; }

void SpydrPick_options::m_init()
{
	namespace po = boost::program_options;

	m_general_options.add_options()
		("help,h", "Print this help message.")
		("version", "Print version information.")
		("verbose,v", po::bool_switch( &SpydrPick_options::s_verbose )->default_value(SpydrPick_options::s_verbose)->notifier(SpydrPick_options::s_init_verbose), "Be verbose.")
		("mi-threshold", po::value< double >( &SpydrPick_options::s_mi_threshold )->default_value(SpydrPick_options::s_mi_threshold), "The MI threshold value. Experience suggests that a value of 0.11 is often reasonable. Zero indicates no threshold and negative values will trigger auto-define heuristics.")
		("mi-pseudocount", po::value< double >( &SpydrPick_options::s_mi_pseudocount )->default_value(SpydrPick_options::s_mi_pseudocount), "The MI pseudocount value.")
	;
	m_parallel_options.add_options()
#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
		("threads,t", po::value< int >( &SpydrPick_options::s_threads )->default_value(SpydrPick_options::s_threads)->notifier(SpydrPick_options::s_init_threads), "Number of threads per MPI/shared memory node (-1=use all hardware threads that the OS/environment exposes).")
#endif // SPYDRPICK_NO_TBB
	;
}

void SpydrPick_options::AddOptions( boost::program_options::options_description *opdesc )
{
	namespace po = boost::program_options;
	opdesc->add(m_general_options);
	opdesc->add(m_parallel_options);
}

/// Check options stored in varmap. Return false if a fatal inconsistency is detected.
bool SpydrPick_options::CheckOptions( boost::program_options::variables_map *varmap )
{
	try
	{
		if( varmap->count("version") && s_out )
		{
			*s_out << s_get_version_string() << std::endl;
			exit(EXIT_SUCCESS);
		}
		if( varmap->count("help") && s_out )
		{
			*s_out << s_title_string << "\n" << s_usage_string << s_options_string << std::endl;
			exit(EXIT_SUCCESS);;
		}
		if( !varmap->count("alignmentfile") && s_err )
		{
			*s_err << "spydrpick ERROR: No alignment file specified!" << std::endl;
			if( s_out )
			{
				*s_out << s_usage_string << std::endl;
			}
			return false;
		}

    }
	catch( std::exception& e)
	{
		if( s_err )
		{
			// probably an unknown option
			*s_err << "spydrpick error: " << e.what() << "\n\n";
			*s_err << s_usage_string << "\n\n";
		}
		return false;
	}
	catch(...)
	{
		if( s_err )
		{
			*s_err << "spydrpick error: Exception of unknown type!\n\n";
		}
		return false;
	}

	return true;
}

void SpydrPick_options::s_init_verbose( const bool& verbose )
{
	if( verbose && s_out )
	{
		*s_out << "spydrpick: being verbose.\n";
	}
}

#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
void SpydrPick_options::s_init_threads( const int& nthreads )
{
	if( s_verbose && s_out )
	{
		*s_out << "spydrpick: user requests " << nthreads << " compute threads.\n";
	}
}
#endif // SPYDRPICK_NO_TBB

} // namespace spydrpick

