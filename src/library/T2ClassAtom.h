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

#ifndef T2CLASSATOM_H_
#define T2CLASSATOM_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include "treedefs.h"
#include "ClassAtom.h"
#include "Exception.h"

/*
 * Class: T2ClassAtom
 *
 * Classifies from any INPUTTYPE to an uli_t OUTPUTTYPE.
 */


#ifndef Tvector
#define Tvector std::vector<T >
#endif
#ifndef CCvector
#define CCvector std::vector<uli_t >
#endif

template <class T>
class T2ClassAtom: public ClassAtom<T, uli_t > {
public:
	T2ClassAtom(): varID(0), threshold(T()) { this->setErrorMissing(2);}
	T2ClassAtom(const T &v ): varID(0), threshold(T(v)) { this->setErrorMissing(2);}
	T2ClassAtom(T val): varID(0), threshold(val) { this->setErrorMissing(2);}

	virtual ~T2ClassAtom(void) { };

  /*
	 * Classifies an input value
	 */
  virtual uli_t classify(Tvector &sample) const {
    if (sample[varID] == this->getMissingCode())
      return this->getErrorMissing();
    return (sample[varID] <= threshold)?0:1;
  }

  virtual uli_t classify(T *sample) const {
    if (rjungle::magicAt(sample, varID) == this->getMissingCode())
      return this->getErrorMissing();
    return (rjungle::magicAt(sample, varID) <= threshold)?0:1;
  }

	/*
	 * Checks the consistency of this vector content which has to be
	 * strictly increasing.
	 */
	bool isConsistent() const {
		return true;
	}

	virtual std::ostream& print(std::ostream& os) const {
		os << "t2class" << "(<=):" << (double)threshold;
		os << "varID:" << varID;
		return os;
	}

	virtual void printXml(std::ostream& os) const {
    os
    << "<classifier id=\"0\" size=\"1\" type=\"t2class\">"
    << "<value varID=\"" << varID << "\">"
    << "</value>" << std::endl
    << "<value id=\"0\">"
    << (double)threshold
    << "</value>" << std::endl
    << "</classifier>" << std::endl;
	}

	inline const T getThreshold() const {return threshold;}
	inline void setThreshold(T val) {this->threshold = val;}

	virtual inline uli_t resultingNodeSize() const {
		return 2;
	}

	inline bool operator==(const T2ClassAtom<T> &t2ClassAtom) {
		return 	t2ClassAtom.threshold == this->threshold;
	}

	virtual inline uli_t getType() const { return ca_T2; }

	inline void setVarID(uli_t varID) {this->varID = varID; }
	inline uli_t getVarID() const {return this->varID;}

	virtual bool isThereVarID(uli_t varID) const {return varID == getVarID();}

	uli_t varID;
	T threshold;
};

/*
 * Source
 */

// ostream
template <class T>
std::ostream& operator<<(std::ostream& os, const T2ClassAtom<T> &t2ClassAtom) {
	return t2ClassAtom.print(os);
}

#endif /*T2CLASSATOM_H_*/
