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

//#include <Eigen/Core>

#include "apegrunt/Alignment.h"
#include "apegrunt/Apegrunt_utility.hpp"
#include "apegrunt/aligned_allocator.hpp"

#include "misc/Coincidence_matrix.hpp"
#include "misc/Matrix_kernel_access_order.hpp"
#include "misc/Vector.h"

#include "accumulators/distribution_std.hpp"
#include "accumulators/distribution_bincount.hpp"
#include "accumulators/distribution_ordered.hpp"
#include "accumulators/distribution_cumulative.hpp"
//#include "accumulators/distribution_generator_svg.hpp"
#include "accumulators/distribution_generator_csv.hpp"

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

	using statepresence_t = typename apegrunt::Alignment<state_t>::statepresence_t;
	using statepresence_block_t = typename apegrunt::Alignment<state_t>::statepresence_block_t;

	using statecount_t = typename apegrunt::Alignment<state_t>::statecount_t;
	using statecount_block_t = typename apegrunt::Alignment<state_t>::statecount_block_t;

	enum { N=apegrunt::number_of_states<state_t>::value };
	using AccessOrder = apegrunt::MATRICES_AccessOrder_tag<N>;
	enum { BlockSize=apegrunt::StateBlock_size };

	mutual_information_block_kernel( apegrunt::Alignment_ptr<state_t> alignment, const std::shared_ptr< std::vector<real_t> > weights, real_t mi_threshold, real_t pseudocount=0.5 )
	: m_alignment(alignment),
	  m_weights(weights),
	  //m_coincidence_block_kernel(alignment,weights,pseudocount),
	  m_buffer( N*N*BlockSize ),
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
	  m_weights( other.m_weights ),
	  //m_coincidence_block_kernel( other.m_coincidence_block_kernel ),
	  m_buffer( N*N*BlockSize ),
	  m_mi_threshold( other.m_mi_threshold ),
	  m_pseudocount( other.m_pseudocount ),
	  m_nloci( other.m_nloci ),
	  m_last_block_index( other.m_last_block_index ),
	  m_last_block_size( other.m_last_block_size ),
	  m_nblocks( other.m_nblocks ),
	  m_neff( other.m_neff )
	{
	}

	// Calculates the mutual information (MI) for all iblock versus jblock data columns, producing a BlockSize-by-BlockSize dense solution matrix.
	inline void operator()( RealT *mi_solution, std::size_t iblock, std::size_t jblock )
    {
		const std::size_t iblock_size = iblock == m_last_block_index ? m_last_block_size : BlockSize;
		const std::size_t jblock_size = jblock == m_last_block_index ? m_last_block_size : BlockSize;
		const std::size_t icol_begin = iblock*BlockSize;
		const std::size_t icol_end = icol_begin+iblock_size;

		std::vector<real_t> joint_H(iblock_size, 0.0);
		std::vector<real_t> icond_H(iblock_size, 0.0);
		std::vector<real_t> jcond_H(iblock_size, 0.0);

		// get statecount tables
		const auto statecount_blocks_ptr = m_alignment->get_statecount_blocks(); // keep shared_ptr alive
		const auto& statecount_blocks = *statecount_blocks_ptr;
		const auto& istatecounts = statecount_blocks[iblock];
		const auto& jstatecounts = statecount_blocks[jblock];

		// get statepresence tables
		const auto statepresence_blocks_ptr = m_alignment->get_statepresence_blocks(); // keep shared_ptr alive
		const auto& statepresence_blocks = *statepresence_blocks_ptr;
		const auto& istatepresence = statepresence_blocks[iblock];
		const auto& jstatepresence = statepresence_blocks[jblock];

		for( std::size_t icol=icol_begin, i=0; icol<icol_end; ++icol, ++i )
		{
			// get the crosstable/coincidence matrix
			auto kernel = apegrunt::weighted_coincidence_block<state_t,real_t>( m_alignment, m_weights, 0.0 );

			// the coincidence kernel collects a tall 1-by-iblock_size coincidence matrix (icol versus all columns in jblock).
			kernel( m_buffer.data(), icol, jblock );

			//std::cout << "mi: normalize icol=" << icol << std::endl;
			// normalize matrices, accounting for pseudocounts

			this->normalize_and_get_jointH_icondH( joint_H.data(), icond_H.data(), istatepresence[icol%BlockSize], jstatepresence, istatecounts[icol%BlockSize], jstatecounts, jblock_size );

			//std::cout << "mi: sum_mat icol=" << icol << std::endl;
			// given the 1-by-iblock_size coincidence matrix, calculate the joint entropy of i and j, H(i,j)
			//this->sum_mat( joint_H.data(), istatepresence[icol%BlockSize], jstatepresence, iblock_size );

			//std::cout << "mi: sum_rows icol=" << icol << std::endl;
			// given the 1-by-iblock_size coincidence matrix, calculate the conditional entropies H(i|j)
			//this->sum_rows( icond_H.data(), istatepresence[icol%BlockSize], jstatepresence, iblock_size );

			//std::cout << "mi: sum_cols icol=" << icol << std::endl;
			// given the 1-by-iblock_size coincidence matrix, calculate the conditional entropies H(j|i)
			this->get_jCondH( jcond_H.data(), istatepresence[icol%BlockSize], jstatepresence, jblock_size );

			//std::cout << "mi: mi_solution icol=" << icol << std::endl;
			// given the 1-by-iblock_size conditional and joint entropies, calculate MI as H(i,j) - H(i|j) - H(j|i)
			this->mi( mi_solution+BlockSize*i, joint_H.data(), icond_H.data(), jcond_H.data(), jblock_size );
			//this->mi( mi_solution+i*iblock_size, joint_H.data(), icond_H.data(), jcond_H.data(), iblock_size );
		}
    }

	inline void single( RealT *mi_solution, std::size_t icol, std::size_t jblock )
    {
		const std::size_t iblock = apegrunt::get_block_index(icol);
		const std::size_t jblock_size = jblock == m_last_block_index ? m_last_block_size : BlockSize;
		const std::size_t jcol_begin = jblock*BlockSize;
		const std::size_t jcol_end = jcol_begin+jblock_size;

		std::vector<real_t> joint_H(jblock_size, 0.0);
		std::vector<real_t> icond_H(jblock_size, 0.0);
		std::vector<real_t> jcond_H(jblock_size, 0.0);

		// get statecount tables
		const auto statecount_blocks_ptr = m_alignment->get_statecount_blocks(); // keep shared_ptr alive
		const auto& statecount_blocks = *statecount_blocks_ptr;
		const auto& istatecounts = statecount_blocks[iblock];
		const auto& jstatecounts = statecount_blocks[jblock];

		// get statepresence tables
		const auto statepresence_blocks_ptr = m_alignment->get_statepresence_blocks(); // keep shared_ptr alive
		const auto& statepresence_blocks = *statepresence_blocks_ptr;
		const auto& istatepresence = statepresence_blocks[iblock];
		const auto& jstatepresence = statepresence_blocks[jblock];

		// get the crosstable/coincidence matrix
		auto kernel = apegrunt::weighted_coincidence_block<state_t,real_t>( m_alignment, m_weights, 0.0 );

		// the coincidence kernel collects a tall 1-by-iblock_size coincidence matrix (icol versus all columns in jblock).
		kernel( m_buffer.data(), icol, jblock );

		this->normalize_and_get_jointH_icondH( joint_H.data(), icond_H.data(), istatepresence[icol%BlockSize], jstatepresence, istatecounts[icol%BlockSize], jstatecounts, jblock_size );

		// given the 1-by-iblock_size coincidence matrix, calculate the conditional entropies H(j|i)
		this->get_jCondH( jcond_H.data(), istatepresence[icol%BlockSize], jstatepresence, jblock_size );

		// given the 1-by-iblock_size conditional and joint entropies, calculate MI as H(i,j) - H(i|j) - H(j|i)
		this->mi( mi_solution, joint_H.data(), icond_H.data(), jcond_H.data(), jblock_size );
    }

	inline void normalize_and_get_jointH_icondH( real_t* const jointH, real_t* const icondH, statepresence_t ipresence, const statepresence_block_t& jpresence_block, statecount_t istatecount, const statecount_block_t& jstatecount_block, std::size_t extent )
	{
		// As the coincidence matrix is always N-by-N, the following function needs to use masking
		// in order to add pseudocounts only to elements present in the crosstable/coincidence matrix.

		const apegrunt::Vector<real_t,N> masked_pseudocount( m_pseudocount, ipresence );
		for( std::size_t i=0; i < extent; ++i )
		{
			const auto jmask = jpresence_block[i];
			//const auto jmask = std::bitset<std::numeric_limits<statepresence_t>::digits>(jpresence_block[i]);

			real_t normconst(0);
			//const real_t normconst( real_t(istatecount*jstatecount_block[i])*m_pseudocount + m_neff ); // didn't work, normconst != sum(elements); apparently due to coincidence matrix round-off errors
			for( std::size_t j=0; j < N; ++j )
			{
				if( jmask & (1<<j) )
				//if( jmask[j] )
				{
					// add pseudocounts and collect row sums
					normconst += apegrunt::mask_sum( apegrunt::make_Vector_view<real_t,N>( m_buffer.data() + AccessOrder::ptr_increment(j,i,extent), 1 ) += masked_pseudocount, ipresence );
				}
			}

			real_t jointHsum(0), icondHsum(0);
			for( std::size_t j=0; j < N; ++j )
			{
				if( jmask & (1<<j) )
				//if( jmask[j] )
				{
					// normalize elements on row j and calculate sum( x * log(x) ) over elements indicated by ipresence mask
					auto row_view = apegrunt::make_Vector_view<real_t,N>( m_buffer.data() + AccessOrder::ptr_increment(j,i,extent), 1 );
					jointHsum += apegrunt::mask_sum_xlogx( row_view /= normconst, ipresence );
					icondHsum += apegrunt::xlogx( apegrunt::sum( row_view ) );
				}
			}
			*(jointH+i) = jointHsum;
			*(icondH+i) = icondHsum;
		}
	}
/*
	inline void sum_mat( real_t* const destination, statepresence_t ipresence, const statepresence_block_t& jpresence_block, std::size_t extent )
	{
		// As the coincidence matrix is always N-by-N, the following function needs to use masking
		// in order apply log() only to elements present in the crosstable/coincidence matrix
		for( std::size_t i=0; i < extent; ++i )
		{
			const auto jmask = std::bitset<std::numeric_limits<statepresence_t>::digits>(jpresence_block[i]);
			apegrunt::Vector<real_t,N> row_sums;
			for( std::size_t j=0; j < N; ++j )
			{
				if( jmask[j] )
				{
					//row_sums[j] = apegrunt::sum( apegrunt::make_Vector_view<real_t,N>( m_buffer.data() + AccessOrder::ptr_increment(j,i,extent), 1 ) );
					row_sums[j] = apegrunt::sum( apegrunt::mask_xlogx( apegrunt::make_Vector_view<real_t,N>( m_buffer.data() + AccessOrder::ptr_increment(j,i,extent), 1 ), ipresence ) );
				}
			}
			*(destination+i) = apegrunt::sum( row_sums );
			//std::cout << "i=" << i << " sum=" << *(destination+i) << std::endl;
		}
	}
*/
/*
	inline void sum_rows( real_t* const destination, statepresence_t ipresence, const statepresence_block_t& jpresence_block, std::size_t extent )
	{
		// As the coincidence matrix is always N-by-N, the following function needs to use masking
		// in order apply log() only to elements present in the crosstable/coincidence matrix

		for( std::size_t i=0; i < extent; ++i )
		{
			const auto jmask = std::bitset<std::numeric_limits<statepresence_t>::digits>(jpresence_block[i]);
			//statepresence_t rowmask(0);
			apegrunt::Vector<real_t,N> row_sums;
			for( std::size_t j=0; j < N; ++j )
			{
				//std::cout << apegrunt::make_Vector_view<real_t,N>( m_buffer.data() + AccessOrder::ptr_increment(j,i,extent), 1 ) << " " << jmask[j] << std::endl;
				if( jmask[j] )
				{
					//rowmask = rowmask | (1 << j);
					row_sums[j] = apegrunt::sum( apegrunt::make_Vector_view<real_t,N>( m_buffer.data() + AccessOrder::ptr_increment(j,i,extent), 1 ) );
				}
			}
			//std::cout << "i=" << i << " jmask=" << jmask << " rowmask=" << std::bitset<std::numeric_limits<statepresence_t>::digits>(rowmask) << " v=" << row_sums << std::endl;
			*(destination+i) = apegrunt::mask_sum_xlogx( row_sums, jpresence_block[i] );
		}
	}
*/
	inline void get_jCondH( real_t* const jcondH, statepresence_t ipresence, const statepresence_block_t& jpresence_block, std::size_t extent )
	{
		// As the coincidence matrix is always N-by-N, the following function needs to use masking
		// in order apply log() only to elements present in the crosstable/coincidence matrix

		//const auto jmask = std::bitset<std::numeric_limits<statepresence_t>::digits>(ipresence);
		for( std::size_t i=0; i < extent; ++i )
		{
			real_t jcondHsum(0);
			for( std::size_t j=0; j < N; ++j )
			{
				if( ipresence & (1<<j) )
				//if( jmask[j] )
				{
					jcondHsum += apegrunt::xlogx( apegrunt::sum( apegrunt::make_Vector_view<real_t,N>( m_buffer.data() + AccessOrder::ptr_increment(0,i,extent)+j, AccessOrder::column_stride(extent) ) ) );
				}
			}
			*(jcondH+i) = jcondHsum;
		}
	}

	inline void mi( real_t* const destination, const real_t* const joint_H, const real_t* const icond_H, const real_t* const jcond_H, std::size_t extent )
	{
		for( std::size_t i=0; i < extent; ++i )
		{
			*(destination+i) = *(joint_H+i) - *(icond_H+i) - *(jcond_H+i);
		}
	}

private:
    apegrunt::Alignment_ptr<state_t> m_alignment;
    std::vector<real_t> m_buffer; // coincidence matrix buffer
    const std::shared_ptr< std::vector<real_t> > m_weights;
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

template< typename RealT, typename StateT > //, typename OptimizerT >
class MI_solver
{
public:
	using real_t = RealT;
	using state_t = StateT;
	using mi_parameters_t = MI_Parameters<state_t,real_t>;

	MI_solver( std::vector< apegrunt::Alignment_ptr<state_t> > alignments,
			//apegrunt::Graph_ptr storage,
			real_t mi_threshold,
			real_t pseudocount=0.5
	)
	: m_mi_parameters( alignments, mi_threshold, pseudocount ),
	  m_storage( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  m_mi_block_kernel( alignments.front(), m_mi_parameters.get_weights(), mi_threshold, pseudocount ),
	  m_cputimer( SpydrPick_options::verbose() ? SpydrPick_options::get_out_stream() : nullptr ),
	  m_solution( apegrunt::StateBlock_size*apegrunt::StateBlock_size, 0 ),
	  m_mi_distribution( acc::tag::distribution::binwidth=0.0001 )
	  //m_loci_slice( loci_slice )
	{
	}

	MI_solver( MI_solver<real_t,state_t>&& other )
    : m_mi_parameters( other.m_mi_parameters ),
	  m_storage( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  //m_storage( other.m_storage ),
	  m_mi_block_kernel( other.m_mi_block_kernel ),
	  m_cputimer( other.m_cputimer ),
	  m_solution( other.m_solution.size(), 0 ), // each instance has its own private solution buffer
	  m_mi_distribution( acc::tag::distribution::binwidth=0.0001 ) // each instance has its own private accumulator
	  //m_loci_slice( std::move( other.m_loci_slice ) )
	{
	}
#ifndef SPYDRPICK_NO_TBB
    // TBB interface (required by tbb::parallel_reduce Body)
	template< typename TBBSplitT >
	MI_solver( MI_solver<real_t,state_t>& other, TBBSplitT s )
    : m_mi_parameters( other.m_mi_parameters ),
	  m_storage( apegrunt::make_Graph_ptr<apegrunt::Graph>() ), // empty network; private for each MI_solver instance
	  //m_storage( other.m_storage ),
	  m_mi_block_kernel( other.m_mi_block_kernel ),
	  m_cputimer( other.m_cputimer ),
	  m_solution( other.m_solution.size(), 0 ), // each instance has its own private solution vector
	  m_mi_distribution( acc::tag::distribution::binwidth=0.0001 ) // each instance has its own private accumulator
	  //m_loci_slice( other.m_loci_slice )
	{
	}
#endif // #ifndef SPYDRPICK_NO_TBB

	// TBB interface (required by tbb::parallel_reduce Body)
	void join( MI_solver<real_t,state_t>& rhs )
	{
		m_storage->join( *(rhs.m_storage) );
		// m_mi_distribution.join( rhs.m_mi_distribution ); // accumulator joining not implemented yet
	}

	inline apegrunt::Graph_ptr get_graph() const { return m_storage; }

	// Calculate MIs for a single position against a block of positions. Self interactions are included.
	inline void operator()( std::size_t icol, std::size_t jblock, bool exclude_selfinteraction=true )
	{
		const std::size_t iblock = apegrunt::get_block_index(icol);
		const std::size_t n_loci = m_mi_parameters.get_alignment()->n_loci(); // cache the number of loci
		auto& storage = *m_storage;
		std::size_t delta = storage.size();

		m_cputimer.start();
		const auto jblock_size = jblock == m_mi_parameters.get_last_block_index() ? m_mi_parameters.get_last_block_size() : m_mi_parameters.get_n_loci_per_block();
		const auto threshold = m_mi_parameters.threshold();

		// reset local solution buffer
		for( auto& e: m_solution ) { e = 0.0; }

		// calculate mutual information for all iblock versus jblock data columns
		m_mi_block_kernel.single( m_solution.data(), icol, jblock ); // calculate MI values in 1-by-16 blocks

		// lock shared storage for thread safety
		//const auto&& scoped_lock = m_storage->acquire_lock(); // lock will expire once out of scope
		storage.lock();
		// store solutions
		if( iblock == jblock && exclude_selfinteraction )
		{
			const std::size_t exclude = icol%apegrunt::StateBlock_size;
			for( std::size_t i = 0; i < jblock_size; ++i )
			{
				if( i != exclude )
				{
					const auto mi = *(m_solution.data()+i);
					//std::cout << " " << mi;
					// Store all solutions that exceed threshold
					if( threshold < mi )
					{
						storage.add( icol, jblock*apegrunt::StateBlock_size+i, mi );
					}
				}
			}
		}
		else
		{
			for( std::size_t i = 0; i < jblock_size; ++i )
			{
				const auto mi = *(m_solution.data()+i);
				// Store all solutions that exceed threshold
				//std::cout << " " << mi;
				if( threshold < mi )
				{
					storage.add( icol, jblock*apegrunt::StateBlock_size+i, mi );
				}
			}
		}
		storage.unlock();

		delta = storage.size() - delta;
		//std::cout << " delta=" << delta << "\n" << std::endl;

		m_cputimer.stop();

		if( SpydrPick_options::verbose() )
		{
			const auto col_begin = 1+jblock*m_mi_parameters.get_n_loci_per_block();
			const auto col_end = col_begin + jblock_size - 1;
			// buffer output in a ss before committing it to the ostream,
			// in order to keep output clean when run in multi-threaded mode.
			std::ostringstream oss;
			oss << "  " << icol << " x " << col_begin << "-" << col_end << " / " << n_loci << " (" << delta << " new edges) time=" << m_cputimer << "\n";
			*SpydrPick_options::get_out_stream() << oss.str();
		}
	}

	// Calculate MIs for a range of blocks against all other blocks in an upper triangular matrix. Self interactions are excluded.
	template< typename RangeT >
    inline void operator()( const RangeT& block_index_range )
    {
		const std::size_t n_loci = m_mi_parameters.get_alignment()->n_loci(); // cache the number of loci
		const auto& index_translation = *(m_mi_parameters.get_alignment()->get_loci_translation());
		//const auto ld_threshold = m_mi_parameters.get_alignment()->get_mi_ld_threshold();

		auto& storage = *m_storage;
		std::size_t delta = storage.size();

		for( const auto iblock: block_index_range )
		{
			m_cputimer.start();
			const auto iblock_size = iblock == m_mi_parameters.get_last_block_index() ? m_mi_parameters.get_last_block_size() : m_mi_parameters.get_n_loci_per_block();
			const auto threshold = m_mi_parameters.threshold();

			// symmetric matrix; compute upper triangular block matrices
			for( std::size_t jblock = iblock; jblock < apegrunt::get_number_of_blocks(n_loci); ++jblock )
			{
				// reset local solution buffer
				for( auto& e: m_solution ) { e = 0.0; }
				// calculate mutual information for all iblock versus jblock data columns
				m_mi_block_kernel( m_solution.data(), iblock, jblock ); // calculate MI values in 16-by-16 blocks
				const auto jblock_size = jblock == m_mi_parameters.get_last_block_index() ? m_mi_parameters.get_last_block_size() : m_mi_parameters.get_n_loci_per_block();

				// lock shared storage for thread safety
				//const auto&& scoped_lock = m_storage->acquire_lock(); // lock will expire once out of scope
				storage.lock();
				// store solutions
				if( iblock == jblock ) // blocks on the diagonal
				{
					for( std::size_t i = 0; i < iblock_size; ++i )
					{
						for( std::size_t j = i+1; j < jblock_size; ++j ) // j = i would be on the diagonal, hence j = i+1
						{
							const auto mi = *(m_solution.data()+iblock_size*i+j);
							//m_mi_distribution(mi);
							// Store all solutions that exceed threshold
							if( threshold < mi )
							{
								storage.add( iblock*apegrunt::StateBlock_size+i, jblock*apegrunt::StateBlock_size+j, mi );
							}
						}
					}
				}
				else
				{
					for( std::size_t i = 0; i < iblock_size; ++i )
					{
						for( std::size_t j = 0; j < jblock_size; ++j )
						{
							const auto mi = *(m_solution.data()+iblock_size*i+j);
							//m_mi_distribution(mi);
							// Store all solutions that exceed threshold
							if( threshold < mi )
							{
								storage.add( iblock*apegrunt::StateBlock_size+i, jblock*apegrunt::StateBlock_size+j, mi );
							}
						}
					}
				}
				storage.unlock();
			}

			delta = storage.size() - delta;

			m_cputimer.stop();

			if( SpydrPick_options::verbose() )
			{
				const auto col_begin = 1+iblock*m_mi_parameters.get_n_loci_per_block();
				const auto col_end = col_begin + iblock_size - 1;
				// buffer output in a ss before committing it to the ostream,
				// in order to keep output clean when run in multi-threaded mode.
				std::ostringstream oss;
				oss << "  " << col_begin << "-" << col_end << " / " << n_loci << " (" << delta << " new edges) time=" << m_cputimer << "\n";
				*SpydrPick_options::get_out_stream() << oss.str();
			}
		}
	}

private:
	mi_parameters_t m_mi_parameters;
	apegrunt::Graph_ptr m_storage;
	mutual_information_block_kernel<state_t, real_t> m_mi_block_kernel;
	stopwatch::stopwatch m_cputimer; // for timing statistics

	using allocator_t = apegrunt::memory::AlignedAllocator<real_t>;
	std::vector<real_t,allocator_t> m_solution;

	acc::accumulator_set<real_t, acc::stats<acc::tag::std(acc::from_distribution),acc::tag::distribution_bincount> > m_mi_distribution;
};

template< typename RealT, typename StateT >
MI_solver<RealT,StateT> get_MI_solver(
	std::vector< apegrunt::Alignment_ptr<StateT> > alignments,
	//apegrunt::Graph_ptr storage,
	RealT mi_threshold,
	RealT pseudocount
) { return MI_solver<RealT,StateT>( alignments, /*storage,*/ mi_threshold, pseudocount ); }

} // namespace spydrpick

#endif // MI_HPP
