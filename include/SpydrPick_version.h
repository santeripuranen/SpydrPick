/** @file spydrpick_version.h

	Copyright (c) 2017-2018 Santeri Puranen.

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License version 3 as
	published by the Free Software Foundation.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Affero General Public License for more details.

	You should have received a copy of the GNU Affero General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.

	@author Santeri Puranen
	$Id: $
*/
#ifndef SPYDRPICK_VERSION_H
#define SPYDRPICK_VERSION_H

namespace spydrpick {

struct spydrpick_version
{
	static const int s_MajorVersion = 0; // substantial rewrite
	static const int s_MinorVersion = 1; // feature change
	static const int s_SubminorVersion = 0; // bugfix, small enhancement
};

} // namespace spydrpick

#endif // SPYDRPICK_VERSION_H
