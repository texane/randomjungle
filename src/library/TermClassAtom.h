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

#ifndef TERMCLASSATOM_H_
#define TERMCLASSATOM_H_


#include <vector>
#include <limits>
#include "ClassAtom.h"
#include "treedefs.h"

template <class T, class C >
class TermClassAtom : public ClassAtom<T, C > {
public:
	TermClassAtom(): val(T()) { };
	TermClassAtom(const TermClassAtom<T, C > &termCA) : val(termCA.val) { };
	TermClassAtom(T val) : val(val) { };

	virtual ~TermClassAtom(void) {
	};

	virtual inline T getResult() const {
		return this->val;
	}

	virtual inline bool isConsistent() const {
		return true;
	}

	virtual std::ostream& print(std::ostream& os) const {
		os << "term:" << (double)this->val;
		return os;
	}

	virtual void printXml(std::ostream& os) const {
    os
    << "<classifier id=\"0\" size=\"1\" type=\"term\">"
    << "<value id=\"0\" type=\"eq\">" << (double)this->val << "</value>"
		<< "</classifier>" << std::endl;
	}

	virtual inline uli_t getType() const { return ca_TERM; }

	T val;

};

#endif /*TERMCLASSATOM_H_*/
