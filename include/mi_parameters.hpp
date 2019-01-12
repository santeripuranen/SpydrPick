/** @file mi_parameters.hpp

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
#ifndef MI_PARAMETERS_HPP
#define MI_PARAMETERS_HPP

#include <memory> // for std::shared_ptr and std::make_shared
#include <vector>

namespace spydrpick {

template< typename StateT, typename RealT >
class MI_Parameters
{
public:
	using real_t = RealT;
	using state_t = StateT;
	using my_type = MI_Parameters<state_t,real_t>;

	using alignment_t = apegrunt::Alignment_ptr<state_t>;
	using alignments_t = std::vector< alignment_t >;

	using weights_t = std::vector<real_t>;
	using weights_ptr = std::shared_ptr< weights_t >;

	MI_Parameters( alignments_t alignments, real_t mi_threshold, real_t mi_pseudocount=0.5 ) :
		m_alignments(alignments),
		m_mi_threshold(mi_threshold),
		m_mi_pseudocount(mi_pseudocount),
		m_r(0),
		m_n_loci(alignments.front()->n_loci()),
		m_n_loci_per_block(apegrunt::StateBlock_size),
		m_last_block_size(apegrunt::get_last_block_size(m_n_loci)),
		m_last_block_index(apegrunt::get_last_block_index(m_n_loci))
	{
		m_weights = std::make_shared< weights_t >();
		auto& weights = *m_weights.get();
		weights.reserve( this->get_alignment()->size() );
		for( auto& seq: this->get_alignment() )
		{
			weights.push_back( seq->weight() );
		}
	}

	//~MI_Parameters() = default;

	MI_Parameters( const my_type& other ) :
		m_alignments( other.m_alignments ),
		m_weights( other.m_weights ),
		m_mi_threshold( other.m_mi_threshold ),
		m_mi_pseudocount( other.m_mi_pseudocount ),
		m_r( other.m_r ),
		m_n_loci( other.m_n_loci ),
		m_n_loci_per_block( other.m_n_loci_per_block ),
		m_last_block_size( other.m_last_block_size ),
		m_last_block_index( other.m_last_block_index )
	{
	}

	alignment_t get_alignment() { return m_alignments.front(); }
	alignments_t get_alignments() { return m_alignments; }

	const weights_ptr get_weights() const { return m_weights; }

	std::size_t number_of_alignments() const { return m_alignments.size(); }

	real_t threshold() const { return m_mi_threshold; }
	real_t pseudocount() const { return m_mi_pseudocount; }

	void set_target_column( std::size_t r ) { m_r = r; }
	std::size_t get_target_column() const { return m_r; }

	std::size_t get_n_loci() const { return m_n_loci; }
	std::size_t get_n_loci_per_block() const { return m_n_loci_per_block; }
	std::size_t get_last_block_size() const { return m_last_block_size; }
	std::size_t get_last_block_index() const { return m_last_block_index; }

private:
	alignments_t m_alignments;
	weights_ptr m_weights;
	real_t m_mi_threshold;
	real_t m_mi_pseudocount;
	std::size_t m_r; // target column

    const std::size_t m_n_loci; // number of columns in the alignment
    const std::size_t m_n_loci_per_block;
    const std::size_t m_last_block_size;
    const std::size_t m_last_block_index;
};

} // namespace spydrpick

#endif // MI_PARAMETERS_HPP


