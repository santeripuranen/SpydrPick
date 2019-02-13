/** @file spydrpick_version.h

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
#ifndef SPYDRPICK_VERSION_H
#define SPYDRPICK_VERSION_H

namespace spydrpick {

struct spydrpick_version
{
	static const int s_MajorVersion = 1; // substantial rewrite
	static const int s_MinorVersion = 1; // feature change
	static const int s_SubminorVersion = 1; // bugfix, small enhancement
};

} // namespace spydrpick

#endif // SPYDRPICK_VERSION_H
