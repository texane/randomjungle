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

#ifndef INIT_H_
#define INIT_H_

#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <streambuf>
#include <cstring>

// Had to be additionally defined,
// because "char" otherwise will not be printed as a number.
//overloaded input operator
template<class T>
std::istream& operator>>(std::istream& is, std::vector<std::vector<T> >& tVec) {
	std::istringstream linestream;
	std::istringstream tokenstream;
	std::string line;
	std::string token;
	size_t pos, pcount, i;
	std::vector<T> result;

	if ((char) is.get() != '(')
		exit(1);
	std::getline(is, line);
	if (*line.rbegin() != ')')
		exit(1);
	linestream.str(line);
	tVec.clear();

	// if empty return to sender
	if ((line.length() > 0) && (line[0] == ')'))
		return is;

	while (!linestream.eof()) {
		std::getline(linestream, token);

		pos = pcount = 0;
		for (i = 0; i < token.size(); ++i) {
			if (token[i] == '(')
				++pcount;
			if ((pcount == 0) && ((token[i] == ')') || (token[i] == ','))) {
				pos = i;
				break;
			}
			if (token[i] == ')')
				--pcount;
		}
		if (pos == 0) {
			std::cerr << "Error in file init.h, in function:" << std::endl;
			std::cerr
					<< "std::istream& operator>>(std::istream& is, std::vector<std::vector<T> >& tVec)"
					<< std::endl;
			std::cout << token << std::endl;
			throw;
		}

		tokenstream.clear();
		tokenstream.str(token.substr(0, pos));
		tokenstream >> result;
		tVec.push_back(result);

		linestream.clear();
		if (pos + 1 != token.size()) {
			linestream.str(token.substr(pos + 1, token.size() - pos - 1));
		} else {
			break;
		}

	}

	return is;
}

// Had to be additionally defined,
// because "char" otherwise will not be printed as a number.
//overloaded input operator
template<class T>
std::istream& operator>>(std::istream& is, std::vector<T>& tVec) {
	std::istringstream linestream;
	std::istringstream tokenstream;
	std::string line;
	std::string token;
	size_t pos, pcount, i;
	T result;
	double pre;

	// get first char
	if ((char) is.get() != '(')
		exit(1);

	// get line (without first char)
	std::getline(is, line);
	if (*line.rbegin() != ')')
		exit(1);
	linestream.str(line);
	tVec.clear();

	// if empty return to sender
	if ((line.length() > 0) && (line[0] == ')'))
		return is;

	while (!linestream.eof()) {
		std::getline(linestream, token);

		pos = pcount = 0;
		for (i = 0; i < token.size(); ++i) {
			if (token[i] == '(')
				++pcount;
			if ((pcount == 0) && ((token[i] == ')') || (token[i] == ','))) {
				pos = i;
				break;
			}
			if (token[i] == ')')
				--pcount;
		}
		if (pos == 0) {
			std::cerr << "Error in file init.h, in function:" << std::endl;
			std::cerr
					<< "std::istream& operator>>(std::istream& is, std::vector<T>& tVec)"
					<< std::endl;
			throw;
		}

		tokenstream.clear();
		tokenstream.str(token.substr(0, pos));
		if (tokenstream.str() == "MAX") {
			tVec.push_back(std::numeric_limits<T>::max());
		} else {
			tokenstream >> pre;
			result = (T) pre;
			tVec.push_back(result);
		}

		linestream.clear();
		if (pos + 1 != token.size()) {
			linestream.str(token.substr(pos + 1, token.size() - pos - 1));
		} else {
			break;
		}

	}

	return is;
}

// Had to be additionally defined,
// because "char" otherwise will not be printed as a number.
template<class T>
std::ostream & operator<<(std::ostream &os, std::vector<std::vector<T> > &vec) {
	typename std::vector<std::vector<T> >::iterator it;

	os << "(";

	it = vec.begin();
	while (it != vec.end()) {
		os << *it;
		++it;
		if (it != vec.end())
			os << ",";
	}

	os << ")";

	return os;
}

// Had to be additionally defined,
// because "char" otherwise will not be printed as a number.
template<class T>
std::ostream & operator<<(std::ostream &os, std::vector<T> &vec) {
	typename std::vector<T>::iterator it;

	os << "(";

	it = vec.begin();
	while (it != vec.end()) {
		if (*it == std::numeric_limits<T>::max())
			os << "MAX";
		else
			os << (double) *it;
		++it;
		if (it != vec.end())
			os << ",";
	}

	os << ")";

	return os;
}

namespace rjungle {
// this def.s make magic access fast
#ifdef __SPARSE_DATA__
static const int magicAtMsk[4] = {192,48,12,3};
static const int magicAtOfs[4] = {6,4,2,0};
inline char magicAt(char* vec, size_t j) {
	return (vec[j / 4] & rjungle::magicAtMsk[j % 4]) >> rjungle::magicAtOfs[j % 4];
}
#else
template<class T>
inline T magicAt(T* vec, uli_t j) {
	return vec[j];
}
#endif

static const char strLom[4][8] = { "undef", "nominal", "ordinal", "numeric" };
}

#endif /* INIT_H_ */
