/** @file aracne_version.h

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
#ifndef ARACNE_VERSION_H
#define ARACNE_VERSION_H

namespace aracne {

struct ARACNE_version
{
	static const int s_MajorVersion = 0; // substantial rewrite
	static const int s_MinorVersion = 1; // feature change
	static const int s_SubminorVersion = 0; // bugfix, small enhancement
};

} // namespace aracne

#endif // ARACNE_VERSION_H
