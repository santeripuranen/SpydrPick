/** @file mi.hpp

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
	$Id: $
*/
#ifndef MI_HPP
#define MI_HPP

#include <memory> // for std::shared_ptr and std::make_shared
#include <functional> // for std::less, std::function

//#include <Eigen/Core>

#include "boost/accumulators/statistics/max.hpp"

#include "apegrunt/Alignment.h"
#include "apegrunt/Apegrunt_utility.hpp"
#include "apegrunt/aligned_allocator.hpp"
#include "apegrunt/Distance.hpp"

#include "misc/Crosstable.hpp"
#include "misc/Matrix_kernel_access_order.hpp"
#include "misc/Vector.h"

//#include "accumulators/distribution_std.hpp"
//#include "accumulators/distribution_bincount.hpp"
//#include "accumulators/distribution_ordered.hpp"
//#include "accumulators/distribution_cumulative.hpp"
//#include "accumulators/distribution_generator_svg.hpp"
//#include "accumulators/distribution_generator_csv.hpp"

#include "graph/Graph.hpp"

#include "mi_parameters.hpp"

namespace spydrpick {

template< typename StateT, typename RealT=double >
class mutual_information_block_kernel
{
public:
	using real_t = RealT;
	using state_t = StateT;
	using my_type = mutual_information_block_kernel<state_t, real_t>;

	using src_ptr_t = const real_t* const;
	using dest_ptr_t = real_t* const;

	using statepresence_t = typename apegrunt::Alignment<state_t>::statepresence_t;
	using statepresence_block_t = typename apegrunt::Alignment<state_t>::statepresence_block_t;

	using statecount_t = typename apegrunt::Alignment<state_t>::statecount_t;
	using statecount_block_t = typename apegrunt::Alignment<state_t>::statecount_block_t;

	enum { N=apegrunt::number_of_states<state_t>::value };
	using AccessOrder = apegrunt::MATRICES_AccessOrder_tag<N>;
	enum { BlockSize=apegrunt::StateBlock_size };

	using edge_t = apegrunt::EdgeID<true>;

	mutual_information_block_kernel( apegrunt::Alignment_ptr<state_t> alignment, const std::shared_ptr< std::vector<real_t> > weights, real_t mi_threshold, real_t pseudocount=0.5 )
	: m_alignment(alignment),
	  m_buffer( N*N*BlockSize*BlockSize ),
	  m_weights(weights),
	  m_kernel( apegrunt::Weighted_crosstable_2Dblock<state_t,real_t>( alignment, weights, 0.0 ) ),
	  m_cached_xtable_id(0,0),
	  m_mi_threshold(mi_threshold),
	  m_pseudocount(pseudocount),
	  m_nloci(alignment->n_loci()),
	  m_last_block_index( apegrunt::get_last_block_index(m_nloci) ),
	  m_last_block_size( apegrunt::get_last_block_size(m_nloci) ),
	  m_nblocks( apegrunt::get_number_of_blocks(m_nloci) ),
	  m_neff( alignment->effective_size() )
	{
	}

	mutual_information_block_kernel( const my_type& other )
	: m_alignment( other.m_alignment ),
	  m_buffer( N*N*BlockSize*BlockSize ),
	  m_weights( other.m_weights ),
	  m_kernel( other.m_kernel ),
	  m_cached_xtable_id(0,0),
	  m_mi_threshold( other.m_mi_threshold ),
	  m_pseudocount( other.m_pseudocount ),
	  m_nloci( other.m_nloci ),
	  m_last_block_index( other.m_last_block_index ),
	  m_last_block_size( other.m_last_block_size ),
	  m_nblocks( other.m_nblocks ),
	  m_neff( other.m_neff )
	{
	}

	inline void block( dest_ptr_t mi_solution, std::size_t iblock, std::size_t jblock, bool exclude_gaps=false )
    {
		const std::size_t iblock_size = iblock == m_last_block_index ? m_last_block_size : BlockSize;
		const std::size_t jblock_size = jblock == m_last_block_index ? m_last_block_size : BlockSize;

		// get statepresence tables
		const auto statepresence_blocks_ptr = ( exclude_gaps ? m_alignment->get_statepresence_blocks_wo_gaps() : m_alignment->get_statepresence_blocks() ); // keep shared_ptr alive
		const auto& statepresence_blocks = *statepresence_blocks_ptr;
		const auto& istatepresence = statepresence_blocks[iblock];
		const auto& jstatepresence = statepresence_blocks[jblock];

		// get the contigency/cross-table

		// the coincidence matrix is cached in m_buffer for short term reuse, e.g. for
		// recalculating MI without the gap contribution
		if( !bool(m_cached_xtable_id) || m_cached_xtable_id != edge_t(iblock,jblock) )
		{
			// get the crosstable/coincidence matrix
			m_kernel( m_buffer.data(), iblock, jblock );

			m_cached_xtable_id = edge_t(iblock,jblock,true);
		}
		for( std::size_t i=0; i<iblock_size; ++i )
		{
			// add pseudocounts, normalize and calculate mi
			this->normalize_and_get_mi_row( m_buffer.data()+N*N*jblock_size*i, mi_solution+i*jblock_size, istatepresence[i], jstatepresence, jblock_size );
		}
    }

	inline void normalize_and_get_mi_row( dest_ptr_t xtables, dest_ptr_t solution,
				statepresence_t ipresence, const statepresence_block_t& jpresence_block, std::size_t extent )
	{
		for( std::size_t k=0; k < extent; ++k ) // k runs over a row of crosstables
		{
			this->normalize_and_get_mi_single( xtables + AccessOrder::ptr_increment(0,k,extent), solution+k, ipresence, jpresence_block[k] );
		}
	}

	inline void normalize_and_get_mi_single(
			dest_ptr_t buffer, dest_ptr_t solution,
			statepresence_t ipresence, statepresence_t jpresence
		)
	{
		const apegrunt::Vector<real_t,N> masked_pseudocount( m_pseudocount, ipresence );

		real_t normconst(0);
		for( std::size_t j=0; j < N; ++j )
		{
			if( apegrunt::is_true( jpresence, j ) )
			{
				// add pseudocounts and collect row sums
				normconst += apegrunt::mask_sum( apegrunt::make_Vector_view<real_t,N>( buffer + j*N, 1 ) += masked_pseudocount, ipresence );
			}
		}

		real_t jointH(0), icondH(0), jcondH(0);
		apegrunt::Vector<real_t,N> jcondHvec(0.0);
		for( std::size_t j=0; j < N; ++j )
		{
			if( apegrunt::is_true( jpresence, j ) )
			{
				// normalize elements on row j and calculate sum( x * log(x) ) over elements indicated by ipresence mask
				// Note: we modify contents of buffer here
				auto row_view = apegrunt::make_Vector_view<real_t,N>( buffer + j*N, 1 );
				jointH += apegrunt::mask_sum_xlogx( row_view /= normconst, ipresence );
				icondH += apegrunt::xlogx( apegrunt::sum( row_view ) );
				jcondHvec += row_view; // here we do N-popcnt(ipresence) unnecessary adds
			}
		}

		jcondH = mask_sum_xlogx( jcondHvec, ipresence );

		*solution = jointH - icondH - jcondH;
	}

	inline void single( dest_ptr_t mi_solution, std::size_t ipos, std::size_t jpos, bool cached=true )
    {
		const std::size_t iblock = apegrunt::get_block_index(ipos);
		const std::size_t jblock = apegrunt::get_block_index(jpos);

		const std::size_t iblock_size = iblock == m_last_block_index ? m_last_block_size : BlockSize;
		const std::size_t jblock_size = jblock == m_last_block_index ? m_last_block_size : BlockSize;

		const std::size_t i = apegrunt::get_pos_in_block(ipos);
		const std::size_t j = apegrunt::get_pos_in_block(jpos);

		const std::size_t bufpos = N*N*jblock_size*i+N*N*j;

		// get statepresence tables
		const auto statepresence_blocks_ptr = m_alignment->get_statepresence_blocks(); // keep shared_ptr alive
		const auto& statepresence_blocks = *statepresence_blocks_ptr;
		const auto& istatepresence = statepresence_blocks[iblock];
		const auto& jstatepresence = statepresence_blocks[jblock];

		// get the contigency/cross-table
		if( cached )
		{
			if( !bool(m_cached_xtable_id) || m_cached_xtable_id != edge_t(iblock,jblock) )
			//if( m_cached_xtable_id != apegrunt::EdgeID(iblock,jblock) )
			{
				// get the crosstable/coincidence matrix
				m_kernel( m_buffer.data(), iblock, jblock );

				// update cache flags
				m_cached_xtable_id = edge_t(iblock,jblock,true);
			}
		}
		else
		{
			m_kernel.single( m_buffer.data()+bufpos, ipos, jpos );
		}

		// given the 1-by-1 coincidence matrix, add pseudocounts & normalize matrix,
		// calculate the joint entropy H(i,j) and the conditional entropies H(i|j) and H(j|i).
		// calculate MI as H(i,j) - H(i|j) - H(j|i)
		this->normalize_and_get_mi_single( m_buffer.data()+bufpos, mi_solution, istatepresence[i/*ipos%BlockSize*/], jstatepresence[j/*jpos%BlockSize*/] );
    }

private:
    apegrunt::Alignment_ptr<state_t> m_alignment;
    std::vector<real_t> m_buffer; // coincidence matrix buffer
    const std::shared_ptr< std::vector<real_t> > m_weights;
    apegrunt::Weighted_crosstable_2Dblock<state_t,real_t> m_kernel;
    edge_t m_cached_xtable_id;
    //apegrunt::weighted_coincidence_block<state_t,real_t> m_coincidence_block_kernel;
    const real_t m_mi_threshold;
    const real_t m_pseudocount;
    const std::size_t m_nloci;
    const std::size_t m_last_block_index;
    const std::size_t m_last_block_size;
    const std::size_t m_nblocks;
    const real_t m_neff;
};

namespace acc = boost::accumulators;

template< typename RealT >
class maxvaltracker
{
public:
	using real_t = RealT;
	using my_type = maxvaltracker<real_t>;

	maxvaltracker() = delete;
	explicit maxvaltracker( std::size_t npos ) : m_npos(npos) { m_bins.resize(m_npos); }

	inline void operator()( std::size_t pos, real_t val ) { m_bins[pos](val); }

	inline void join( const my_type& other )
	{
		if( m_npos == other.m_npos )
		{
			for( std::size_t i=0; i < m_npos; ++i )
			{
				m_bins[i]( acc::max(other.m_bins[i]) );
			}
		}
	}

	template< std::size_t Q >
	real_t quartile()
	{
		using std::begin; using std::end;
		if( m_vals.size() == 0 )
		{
			m_vals.reserve(m_npos);
			for( const auto& bin: m_bins )
			{
				m_vals.push_back( acc::max(bin) );
			}
			std::sort( begin(m_vals), end(m_vals), std::less<real_t>() );
		}
		return m_vals[m_npos/4*Q];
	}

	inline std::size_t size() const { return m_bins.size(); }

private:
	const std::size_t m_npos;
	std::vector< acc::accumulator_set<real_t, acc::stats<acc::tag::max> > > m_bins;
	std::vector<real_t> m_vals;

};

template< typename RealT, typename StateT >
class MI_solver
{
public:
	using real_t = RealT;
	using state_t = StateT;
	using mi_parameters_t = MI_Parameters<state_t,real_t>;

	MI_solver( std::vector< apegrunt::Alignment_ptr<state_t> > alignments,
			real_t mi_threshold,
			real_t pseudocount=0.5
	)
	: m_mi_parameters( alignments, mi_threshold, pseudocount ),
	  m_stored_edges( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  m_stored_edges_wo_gaps( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  m_mi_block_kernel( alignments.front(), m_mi_parameters.get_weights(), mi_threshold, pseudocount ),
	  m_cputimer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ),
	  m_solution( apegrunt::StateBlock_size*apegrunt::StateBlock_size, 0 ),
	  //m_mi_distribution( acc::tag::distribution::binwidth=0.0001 ),
	  m_colmax( alignments.front()->n_loci() )
	{
		if( apegrunt::Apegrunt_options::linear_genome() )
		{
			m_distance = apegrunt::GenomeDistance<apegrunt::LinearDistance>( m_mi_parameters.get_alignment() );
		}
		else
		{
			m_distance = apegrunt::GenomeDistance<apegrunt::CircularDistance>( m_mi_parameters.get_alignment() );
		}
	}

	MI_solver( MI_solver<real_t,state_t>&& other )
    : m_mi_parameters( other.m_mi_parameters ),
	  m_stored_edges( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  m_stored_edges_wo_gaps( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  m_mi_block_kernel( other.m_mi_block_kernel ),
	  m_cputimer( other.m_cputimer ),
	  m_solution( other.m_solution.size(), 0 ), // each instance has its own private solution buffer
	  //m_mi_distribution( acc::tag::distribution::binwidth=0.0001 ), // each instance has its own private accumulator
	  m_colmax( other.m_colmax.size() ), // each instance has its own private accumulator
	  m_distance( other.m_distance )
	  //m_loci_slice( std::move( other.m_loci_slice ) )
	{
	}
#ifndef SPYDRPICK_NO_TBB
    // TBB interface (required by tbb::parallel_reduce Body)
	template< typename TBBSplitT >
	MI_solver( MI_solver<real_t,state_t>& other, TBBSplitT s )
    : m_mi_parameters( other.m_mi_parameters ),
	  m_stored_edges( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  m_stored_edges_wo_gaps( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  m_mi_block_kernel( other.m_mi_block_kernel ),
	  m_cputimer( other.m_cputimer ),
	  m_solution( other.m_solution.size(), 0 ), // each instance has its own private solution vector
	  //m_mi_distribution( acc::tag::distribution::binwidth=0.0001 ), // each instance has its own private accumulator
	  m_colmax( other.m_colmax.size() ), // each instance has its own private accumulator
	  m_distance( other.m_distance )
	  //m_loci_slice( other.m_loci_slice )
	{
	}
#endif // #ifndef SPYDRPICK_NO_TBB

	// TBB interface (required by tbb::parallel_reduce Body)
	void join( MI_solver<real_t,state_t>& rhs )
	{
		m_stored_edges->join( *(rhs.m_stored_edges) );
		m_stored_edges_wo_gaps->join( *(rhs.m_stored_edges_wo_gaps) );
		// m_mi_distribution.join( rhs.m_mi_distribution ); // accumulator joining not implemented yet
		m_colmax.join( rhs.m_colmax );
	}

	inline apegrunt::Graph_ptr get_graph() const { return m_stored_edges; }
	inline apegrunt::Graph_ptr get_graph_wo_gaps() const { return m_stored_edges_wo_gaps; }

	// Calculate MIs for a pair of positions.
	inline real_t single( std::size_t ipos, std::size_t jpos )
	{
		real_t mi;
		m_mi_block_kernel.single( &mi, ipos, jpos, false ); // false == crosstable not cached
		return mi;
	}

	inline maxvaltracker<real_t>& get_quartiles() { return m_colmax; }

	// Calculate MIs for a range of blocks against all other blocks in an upper triangular matrix. Self interactions are excluded.
	template< typename RangeT >
    inline void operator()( const RangeT& block_index_range )
    {
		const std::size_t n_loci = m_mi_parameters.get_alignment()->n_loci(); // cache the number of loci
		const auto gappresence_blocks_ptr = m_mi_parameters.get_alignment()->get_gappresence_blocks(); // keep shared_ptr alive
		const auto& gappresence_blocks = *gappresence_blocks_ptr;

		auto& stored_edges = *m_stored_edges;
		std::size_t edges_added = stored_edges.size();

		std::vector< apegrunt::EdgeID<true> > gap_edges;
		gap_edges.reserve( apegrunt::StateBlock_size*apegrunt::StateBlock_size );

		for( const auto iblock: block_index_range )
		{
			m_cputimer.start();
			const auto iblock_size = iblock == m_mi_parameters.get_last_block_index() ? m_mi_parameters.get_last_block_size() : m_mi_parameters.get_n_loci_per_block();
			const auto threshold = m_mi_parameters.threshold();
			const auto& igappresence = gappresence_blocks[iblock];

			// symmetric matrix; compute upper triangular block matrices
			for( std::size_t jblock = iblock; jblock < apegrunt::get_number_of_blocks(n_loci); ++jblock )
			{
				// reset local solution buffer
				for( auto& e: m_solution ) { e = 0.0; }

				const auto& jgappresence = gappresence_blocks[jblock];

				// calculate mutual information for all iblock versus jblock data columns
				m_mi_block_kernel.block( m_solution.data(), iblock, jblock ); // calculate MI values in 16-by-16 blocks
				const auto jblock_size = jblock == m_mi_parameters.get_last_block_index() ? m_mi_parameters.get_last_block_size() : m_mi_parameters.get_n_loci_per_block();

				// lock shared storage for thread safety
				//const auto&& scoped_lock = m_storage->acquire_lock(); // lock will expire once out of scope
				stored_edges.lock(); // lock
				std::size_t delta = stored_edges.size();
				// store solutions
				if( iblock == jblock ) // blocks on the diagonal
				{
					for( std::size_t i = 0; i < iblock_size; ++i )
					{
						const auto ipos = iblock*apegrunt::StateBlock_size+i;
						for( std::size_t j = i+1; j < jblock_size; ++j ) // j = i would be on the diagonal, hence j = i+1
						{
							const auto mi = *(m_solution.data()+jblock_size*i+j);
							const auto jpos = jblock*apegrunt::StateBlock_size+j;
							if( m_distance(ipos,jpos) > SpydrPick_options::get_ld_threshold() )
							{
								m_colmax( ipos, mi );
								m_colmax( jpos, mi );
							}
							//m_mi_distribution(mi);
							// Store all solutions that exceed threshold
							if( threshold < mi )
							{
								stored_edges.add( ipos, jpos, mi );
								if( igappresence[i] || jgappresence[j] ) { gap_edges.emplace_back(i,j); }
							}
						}
					}
				}
				else
				{
					for( std::size_t i = 0; i < iblock_size; ++i )
					{
						const auto ipos = iblock*apegrunt::StateBlock_size+i;
						for( std::size_t j = 0; j < jblock_size; ++j )
						{
							const auto mi = *(m_solution.data()+jblock_size*i+j);
							const auto jpos = jblock*apegrunt::StateBlock_size+j;
							if( m_distance(ipos,jpos) > SpydrPick_options::get_ld_threshold() )
							{
								m_colmax( ipos, mi );
								m_colmax( jpos, mi );
							}
							//m_mi_distribution(mi);
							// Store all solutions that exceed threshold
							if( threshold < mi )
							{
								stored_edges.add( ipos, jpos, mi );
								if( igappresence[i] || jgappresence[j] ) { gap_edges.emplace_back(i,j); }
							}
						}
					}
				}
				delta = stored_edges.size() - delta;
				stored_edges.unlock(); // explicit unlock

				{
					// if we've stored MIs and the active node blocks contain gaps, then re-evaluate edges while excluding the gap contribution
					if( !gap_edges.empty() )
					{
						// reset local solution buffer
						for( auto& e: m_solution ) { e = 0.0; }
						// calculate a block of MI values, excluding gap contribution; this call should re-use the cached crosstable
						m_mi_block_kernel.block( m_solution.data(), iblock, jblock, true ); // true == exclude the gap contribution

						m_stored_edges_wo_gaps->lock();
						for( const auto& edge: gap_edges )
						{
							//const auto i = pair.first;
							//const auto j = pair.second;
							const auto i = edge.first();
							const auto j = edge.second();
							const auto ipos = iblock*apegrunt::StateBlock_size+i;
							const auto jpos = jblock*apegrunt::StateBlock_size+j;
							const auto mi = *(m_solution.data()+jblock_size*i+j);
							m_stored_edges_wo_gaps->add( ipos, jpos, mi );
						}

						m_stored_edges_wo_gaps->unlock();
						gap_edges.clear();
					}
				}
			}

			edges_added = stored_edges.size() - edges_added;

			m_cputimer.stop();

			if( SpydrPick_options::verbose() )
			{
				const auto col_begin = 1+iblock*m_mi_parameters.get_n_loci_per_block();
				const auto col_end = col_begin + iblock_size - 1;
				// buffer output in a ss before committing it to the ostream,
				// in order to keep output clean when run in multi-threaded mode.
				std::ostringstream oss;
				oss << "  " << col_begin << "-" << col_end << " / " << n_loci << " (" << edges_added << " new edges) time=" << m_cputimer << "\n";
				*SpydrPick_options::get_out_stream() << oss.str();
			}
		}
	}

private:
	mi_parameters_t m_mi_parameters;
	apegrunt::Graph_ptr m_stored_edges;
	apegrunt::Graph_ptr m_stored_edges_wo_gaps;

	mutual_information_block_kernel<state_t, real_t> m_mi_block_kernel;
	stopwatch::stopwatch m_cputimer; // for timing statistics

	using allocator_t = apegrunt::memory::AlignedAllocator<real_t>;
	std::vector<real_t,allocator_t> m_solution;
	maxvaltracker<real_t> m_colmax;
	std::function<std::size_t(std::size_t,std::size_t)> m_distance;

	//acc::accumulator_set<real_t, acc::stats<acc::tag::std(acc::from_distribution),acc::tag::distribution_bincount> > m_mi_distribution;
};

template< typename RealT, typename StateT >
MI_solver<RealT,StateT> get_MI_solver(
	std::vector< apegrunt::Alignment_ptr<StateT> > alignments,
	//apegrunt::Graph_ptr storage,
	RealT mi_threshold,
	RealT pseudocount
) { return MI_solver<RealT,StateT>( alignments, mi_threshold, pseudocount ); }

} // namespace spydrpick

#endif // MI_HPP
