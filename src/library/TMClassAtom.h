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

#ifndef TMCLASSATOM_H_
#define TMCLASSATOM_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include "treedefs.h"
#include "ClassAtom.h"
#include "Exception.h"

/*
 * Class: TMClassAtom
 *
 * Classifies from any INPUTTYPE to an uli_t OUTPUTTYPE.
 */


#ifndef Tvector
#define Tvector std::vector<T >
#endif
#ifndef CCvector
#define CCvector std::vector<uli_t >
#endif
#ifndef MYSIZETINF
#define MYSIZETINF std::numeric_limits<uli_t>::max()
#endif

template <class T>
class TMClassAtom: public ClassAtom<T, uli_t > {
public:
	TMClassAtom() { this->setErrorMissing(MYSIZETINF);}
	TMClassAtom(std::vector<T > &vec): thresholds(vec) { this->setErrorMissing(MYSIZETINF);}

	virtual ~TMClassAtom(void) { };

	/*
	 * Classifies an input value
	 */
  virtual uli_t classify(Tvector &sample) const {
    if (sample[varID] == this->getMissingCode()) return this->getErrorMissing();

    typename Tvector::const_iterator it(this->thresholds.begin());
    typename Tvector::const_iterator it2(this->thresholds.begin());
    uli_t cla = 0;
    while(it2 != this->thresholds.end()) {
      if(it2 == this->thresholds.begin()) {
        if (sample[varID] <= (*it)) {
          break;
        } else {
          ++cla;
          ++it2;
          continue;
        }
      }
      if ( sample[varID] > (*it) && sample[varID] <= (*it2)) break;
      ++cla;
      ++it;
      ++it2;
      //if (it2 == this->thresholds.end()) --cla;
    }

    return cla;
  }

  virtual uli_t classify(T *sample) const {
    if (rjungle::magicAt(sample, varID) == this->getMissingCode())
      return this->getErrorMissing();

    typename Tvector::const_iterator it(this->thresholds.begin());
    typename Tvector::const_iterator it2(this->thresholds.begin());
    uli_t cla = 0;
    while(it2 != this->thresholds.end()) {
      if(it2 == this->thresholds.begin()) {
        if (rjungle::magicAt(sample, varID) <= (*it)) {
          break;
        } else {
          ++cla;
          ++it2;
          continue;
        }
      }
      if (rjungle::magicAt(sample, varID) > (*it) &&
        rjungle::magicAt(sample, varID) <= (*it2)) break;
      ++cla;
      ++it;
      ++it2;
      //if (it2 == this->thresholds.end()) --cla;
    }

    return cla;
  }



	/*
	 * Checks the consistency of this vector content which has to be
	 * strictly increasing.
	 */
	bool isConsistent() const {
		return true;
	}

	virtual std::ostream& print(std::ostream& os) const {
    typename Tvector::const_iterator it(this->thresholds.begin());

    os << "tmclass(<=):";
    if (it == this->thresholds.end()) return os;
    os << *it;
    ++it;
		while(it != this->thresholds.end()) {
		  os << "," << (double)(*it);
			++it;
		}

		return os;
	}

	virtual void printXml(std::ostream& os) const {
		typename Tvector::const_iterator it(this->thresholds.begin());
		uli_t i = 0;

    os
    << "<classifier size=\"" << this->thresholds.size()
    << "\" type=\"tmclass\">";

    //++it;
		while(it != this->thresholds.end()) {
		  os
      << "<value id=\"" << i << "\">"
      << (double)(*it) << "</value>";
      ++i;
			++it;
		}
		os  << std::endl << "</classifier>" << std::endl;
	}

	inline Tvector &getThresholds() {return thresholds;}
	void setThresholds(Tvector vec) {
    thresholds.assign(vec.begin(), vec.end());
    sort(thresholds.begin(), thresholds.end());
  }

	void setThresholdsNoSort(Tvector vec) {
    thresholds.assign(vec.begin(), vec.end());
  }

	virtual inline uli_t resultingNodeSize() const {
		return thresholds.size() + 1;
	}

	virtual uli_t size() const {
		return resultingNodeSize();
	}

	inline bool operator==(const TMClassAtom<T > &tMClassAtom) {
		return 	tMClassAtom.getThresholds() == this->getThresholds();
	}

	virtual inline uli_t getType() const { return ca_TM; }

  inline void setVarID(uli_t varID) {this->varID = varID; }
  inline uli_t getVarID() const {return this->varID;}

  virtual bool isThereVarID(uli_t varID) const {return varID == getVarID();}

  uli_t varID;
private:
	Tvector thresholds;
};

/*
 * Source
 */

// ostream
template <class T>
std::ostream& operator<<(std::ostream& os, const TMClassAtom<T> &tMClassAtom) {
	return tMClassAtom.print(os);
}

#endif /*TMCLASSATOM_H_*/
