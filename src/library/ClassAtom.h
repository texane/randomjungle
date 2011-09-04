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


#ifndef CLASSATOM_H_
#define CLASSATOM_H_


/*
 * Includes
 */

#include <iostream>
#include <vector>
#include "treedefs.h"

/** 
 * \brief Smallest classifier unit. It classify a "value"
 * with a low level function.
 * I.e. classification via threshold or logic regression.
 */
template <class T, class C>
class ClassAtom{
public:

	/** 
	 * \brief Constructor that sets the internal missing code.
	 */
	ClassAtom(): missingCode(MISSINGCODE) { };

	/** 
	 * \brief Constructor that sets the internal missing code.
	 *
	 * @param Object of classAtom (tree node classifier).
	 */
	ClassAtom(const ClassAtom<T, C > &classAtom): missingCode(MISSINGCODE) { };
	virtual ~ClassAtom() { };

	/** 
	 * \brief Classifies an input value
	 * 
	 * @param vecIn Sample to be classified
	 * 
	 * @return Classification value
	 */
	virtual C classify(const std::vector<T > &vecIn) const {
	  return C();
	}

	/** 
	 * \brief Classifies an input value
	 * 
	 * @param vecIn Sample to be classified
	 * 
	 * @return Classification value
	 */
	virtual C classify(T *vecIn) const {
	  return C();
	}

	/** 
	 * \brief Returns the value of this classifier.
	 * 
	 * @return Classifier value
	 */
	virtual inline T getResult() const {
		return 0;
	}

	/** 
	 * \brief Returns the size of classifier.
	 * 
	 * @return Size of classifier.
	 */
	virtual uli_t size() const {
		return 0;
	}

	/** 
	 * \brief Returns the size of resulting node.
	 * 
	 * @return Size of resulting node.
	 */
	virtual inline uli_t resultingNodeSize() const {return 0; }

	/** 
	 * \brief Returns if the classifier is consistent.
	 * 
	 * @return If the classifier is consistent.
	 */
	virtual bool isConsistent() const {
		return false;
	}

	/** 
	 * \brief Prints the classifier to output stream.
	 * 
	 * @param os Output stream.
	 * 
	 * @return Output stream.
	 */
	virtual std::ostream& print(std::ostream& os) const {
		return os;
	}

	/** 
	 * \brief Prints the classifier to output stream in XML format.
	 * 
	 * @param os Output stream.
	 * 
	 * @return Output stream.
	 */
	virtual void printXml(std::ostream& os) const {
	}

	/** 
	 * \brief Returns the type of current classifier.
	 * 
	 * @return Classifier type.
	 */
	virtual inline uli_t getType() const { return ca_BASE; }

	/** 
	 * \brief Sets the value that means "missing value".
	 * 
	 * @param val Missing value coding.
	 */
	virtual inline void setMissingCode(T val) { missingCode = val; }

	/** 
	 * \brief Returns the missing value coding.
	 * 
	 * @return Missing value coding.
	 */
	virtual inline T getMissingCode() const { return missingCode; }

	/** 
	 * \brief Sets the value that means "error there is a missing value".
	 * 
	 * @param val Error missing value coding.
	 */
	virtual inline void setErrorMissing(C val) { errorMissing = val; }

	/** 
	 * \brief Returns the error missing value coding.
	 * 
	 * @return Error missing value coding.
	 */
	virtual inline C getErrorMissing() const { return errorMissing; }

	/** 
	 * \brief Return if classifier has got a specific variable ID.
	 * 
	 * @param varID ID of variable.
	 * 
	 * @return 
	 */
	virtual bool isThereVarID(uli_t varID) const { return false; }

	/// Missing code
	T missingCode;

	/// Error missing code
  C errorMissing;
};

/*
 * Source
 */

/** 
 * \brief Prints the classifier to output stream.
 * 
 * @param os Output stream.
 * @param classAtom Classifier object.
 * 
 * @return Output stream.
 */
template <class T, class C>
std::ostream& operator<<(std::ostream& os, const ClassAtom<T, C> &classAtom) {
	return classAtom.print(os);
}

#endif /*CLASSATOM_H_*/
