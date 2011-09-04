/*
 * Copyright (C) 2008-2010  Daniel F. Schwarz
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

/*
 * Includes
 */
#include <string>

/*
 * Class declaration
 */

class Exception //: public std::exception
{
public:

	Exception(std::string str);
	Exception(int status, std::string str) :
		str(str), status(status) {
	}
	;

	virtual ~Exception() {
	}
	;

	virtual const char* what() const throw ();

	std::string str;
	int status;

};

#endif /*EXCEPTION_H_*/
