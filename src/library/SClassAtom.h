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

#ifndef SCLASSATOM_H_
#define SCLASSATOM_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include "treedefs.h"
#include "ClassAtom.h"
#include "Exception.h"

/*
 * Class: SClassAtom
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
class SClassAtom: public ClassAtom<T, uli_t > {
public:
	SClassAtom(): varID(0) { this->setErrorMissing(2);}
  SClassAtom(std::vector< T> &v): varID(0), classSet(v) { this->setErrorMissing(2);}

	virtual ~SClassAtom(void) { };

  /*
	 * Classifies an input
	 */
  virtual uli_t classify(Tvector &sample) const {
    if (sample[varID] == this->getMissingCode())
      return this->getErrorMissing();

    uli_t branch = 1;
    for (uli_t counter = 0; counter < classSet.size(); ++counter)
      if (sample[varID] == classSet[counter]) {
        branch = 0;
        break;
      }

    return branch;
  }

  virtual uli_t classify(T *sample) const {
    if (rjungle::magicAt(sample, varID) == this->getMissingCode())
      return this->getErrorMissing();

    uli_t branch = 1;
    for (uli_t counter = 0; counter < classSet.size(); ++counter)
      if (rjungle::magicAt(sample, varID) == classSet[counter]) {
        branch = 0;
        break;
      }

    return branch;
  }

	/*
	 * Checks the consistency of this vector content which has to be
	 * strictly increasing.
	 */
	bool isConsistent() const {
		return true;
	}

	virtual std::ostream& print(std::ostream& os) const {
    os << "sclass" << "(element of):";
    for (uli_t counter = 0; counter < classSet.size(); ++counter)
      os << classSet[counter] << " ";
    os << std::endl;
		return os;
	}

	virtual void printXml(std::ostream& os) const {
    os
    << "<classifier id=\"" << ca_S
    << "\" size=\"" << 0
    << "\" varID=\"" << varID
    << "\" type=\"sclass\">Not impl. yet</classifier>" << std::endl;
    /*
    os
    << "<classifier id=\"0\" size=\"1\" type=\"t2class\">"
    << "<value id=\"0\">"
    << (double)threshold
    << "</value>" << std::endl
    << "</classifier>" << std::endl;
    */
	}

  inline const std::vector< T> &getClassSet() const {return classSet;}
  inline void setClassSet(std::vector< T> &set) {this->classSet = set;}

	virtual inline uli_t resultingNodeSize() const {
		return 2;
	}

	inline bool operator==(const SClassAtom<T> &sClassAtom) {
    return 	sClassAtom.classSet == this->classSet;
	}

	virtual inline uli_t getType() const { return ca_S; }

	inline void setVarID(uli_t varID) {this->varID = varID; }
	inline uli_t getVarID() const {return this->varID;}

	virtual bool isThereVarID(uli_t varID) const {return varID == getVarID();}

	uli_t varID;
  std::vector<T > classSet;
};

/*
 * Source
 */

// ostream
template <class T>
std::ostream& operator<<(std::ostream& os, const SClassAtom<T> &sClassAtom) {
	return sClassAtom.print(os);
}

#endif /*SCLASSATOM_H_*/
