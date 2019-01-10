/** @file spydrpick.h
	Top-level include file for spydrpick.

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
#ifndef SPYDRPICK_H
#define SPYDRPICK_H

#ifndef SPYDRPICK_NO_TBB // Threading with Threading Building Blocks
#include "tbb/tbb_stddef.h"
#include "tbb/task_scheduler_init.h"
#endif // #ifndef SPYDRPICK_NO_TBB

#include "apegrunt/Alignment_forward.h"
#include "apegrunt/StateVector_forward.h"
#include "apegrunt/StateVector_state_types.hpp"
#include "apegrunt/StateVector_utility.hpp"
#include "apegrunt/Alignment_parsers.hpp"
#include "apegrunt/Alignment_impl_block_compressed_storage.hpp"
#include "apegrunt/StateVector_impl_block_compressed_alignment_storage.hpp"

#include "SpydrPick_options.h"

#endif // SPYDRPICK_H
