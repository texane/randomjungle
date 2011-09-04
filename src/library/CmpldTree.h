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

#ifndef CMPLDTREE_H_
#define CMPLDTREE_H_

/*
 * Includes
 */

#include <vector>
#include <limits>

#include "init.h"
#include "treedefs.h"

/*
 * CmpldTree
 */
template<class T>
class CmpldTree {
public:

	/** 
	 * \brief Construct that sets several tree parameters.
	 * 
	 * @param varSize Size of variable structure.
	 * @param valueWidth Width of value structure.
	 * @param valueSize Size of value structure.
	 * @param branchSize Size of branch structure.
	 */
	CmpldTree(uli_t varSize, uli_t valueWidth, uli_t valueSize, uli_t branchSize) :
		varID(std::vector<std::vector<uli_t> >(varSize, std::vector<uli_t>())),
				varSize(varSize), valueSize(valueSize), valueWidth(valueWidth),
				branchSize(branchSize), termID(std::numeric_limits<uli_t>::max()) {
	}
	;

	/** 
	 * \brief Destructor that does nothing.
	 */
	virtual ~CmpldTree() {
	}
	;

	/** 
	 * \brief Prints the tree to output stream.
	 * 
	 * @param output Output stream.
	 */
	void printXml(std::ostream &output) {
		variableToXML<std::vector<std::vector<uli_t> > > (varID, "varID", output);
		variableToXML<std::vector<std::vector<std::vector<T> > > > (values,
				"values", output);
		variableToXML<std::vector<std::vector<uli_t> > > (branches, "branches",
				output);
		variableToXML<std::vector<T> > (classes, "classes", output);
		variableToXML<std::vector<uli_t> > (indexes, "indexes", output);
		variableToXML<std::vector<std::vector<double> > > (extra, "extra", output);

		variableToXML<uli_t> (varSize, "varSize", output);
		variableToXML<uli_t> (valueWidth, "valueWidth", output);
		variableToXML<uli_t> (valueSize, "valueSize", output);
		variableToXML<uli_t> (branchSize, "branchSize", output);
		variableToXMLWithMax<uli_t> (termID, "termID", output);
	}

	/** 
	 * \brief Prints the tree to output stream.
	 * 
	 * @param os Output stream.
	 * @param isXML Is it in XML format.
	 * 
	 * @return Output stream.
	 */
	std::ostream &print(std::ostream &os, bool isXML = false) const {
		uli_t i, j;

		// print nodes
		for (i = 0; i < values.size(); ++i) {
			// print var ID
			if (varID.size() != 0)
				os << varID[0][j];
			for (j = 1; j < varID.size(); ++j) {
				os << "\t" << varID[j][i];
			}

			// seperator
			os << "\t|\t";

			// print class values
			os << this->values[i];

			// seperator
			os << "\t|\t";

			// print branches
			if (branches[i].size() != 0)
				os << branches[i][0];
			for (j = 1; j < branches[i].size(); ++j) {
				os << "\t" << branches[i][j];
			}

			isXML ? (os << ";") : (os << std::endl);
		}

		// classes
		os << "classes:";
		for (j = 0; j < classes.size(); ++j) {
			os << " " << classes[j];
		}
		isXML ? (os << ";") : (os << std::endl);

		// indexes
		os << "indexes:";
		for (j = 0; j < indexes.size(); ++j) {
			os << " " << indexes[j];
		}
		isXML ? (os << ";") : (os << std::endl);

		// extra
		os << "extra: (elements not shown)";
		isXML ? (os << ";") : (os << std::endl);

		return os;
	}

	/** 
	 * \brief Appends a clean node to flattend structure.
	 */
	void pushBackNode() {
		for (uli_t i = 0; i < varSize; ++i)
			varID[i].push_back(0);
		values.push_back(std::vector<std::vector<T> >(valueWidth, std::vector<T>(
				valueSize, 0)));
		branches.push_back(std::vector<uli_t>(branchSize, 0));
	}

	/** 
	 * \brief Sets variable ID of last node.
	 * 
	 * @param id The ID.
	 */
	inline void toLastVarID(uli_t id) {
		*(varID[0].rbegin()) = id;
	}

	/** 
	 * \brief Sets multidimensional variable of last node.
	 * 
	 * @param i Dimension position.
	 * @param id The ID.
	 */
	inline void toLastVarID(uli_t i, uli_t id) {
		*(varID[i].rbegin()) = id;
	}

	/** 
	 * \brief Sets value of last node.
	 * 
	 * @param i Index of value.
	 * @param val Value.
	 */
	inline void toLastValues(uli_t i, T val) {
		values.rbegin()->at(0).at(i) = val;
	}

	/** 
	 * \brief Appends value to last value vector.
	 * 
	 * @param val
	 */
	inline void pushBackToLastValues(T val) {
		values.rbegin()->at(0).push_back(val);
	}

	/** 
	 * \brief Appends value to last value vector i.
	 * 
	 * @param i Index of value in vector.
	 * @param val
	 */
	inline void pushBackToLastValues(uli_t i, T val) {
		values.rbegin()->at(i).push_back(val);
	}

	/** 
	 * \brief Sets a branch value to a specific value.
	 * 
	 * @param i Index in last branch.
	 * @param val Value.
	 */
	inline void toLastBranches(uli_t i, uli_t val) {
		branches.rbegin()->at(i) = val;
	}

	/** 
	 * \brief Sets ID of a variable
	 * 
	 * @param at Position.
	 * @param id ID of variable.
	 */
	inline void toVarID(uli_t at, uli_t id) {
		varID[0][at] = id;
	}

	/** 
	 * \brief Sets value of a node.
	 * 
	 * @param at Position in list of nodes.
	 * @param i Index in value vector.
	 * @param val Value.
	 */
	inline void toValues(uli_t at, uli_t i, T val) {
		values[at][0][i] = val;
	}

	/** 
	 * \brief Sets branch of a node.
	 * 
	 * @param at Position in list of nodes.
	 * @param i Index in branch vector.
	 * @param val Value.
	 */
	inline void toBranches(uli_t at, uli_t i, uli_t val) {
		branches[at][i] = val;
	}

	/** 
	 * \brief Gets size of variable ID vector.
	 * 
	 * @return Size of variable ID vector.
	 */
	inline uli_t size() {
		return varID[0].size();
	}

	/** 
	 * \brief Returns if a specific variable ID is within the tree.
	 * 
	 * @param id The ID of variable.
	 * 
	 * @return If variable is within the tree.
	 */
	inline bool hasVarID(unsigned int id) {
		return (std::find(varID[0].begin(), varID[0].end(), id) != varID[0].end());
	}

	/** 
	 * \brief Prints a variable to output stream.
	 * 
	 * @param val Value.
	 * @param name Name of variable.
	 * @param out Output stream.
	 */
	template<class TT>
	void variableToXML(TT &val, std::string name, std::ostream &out) {
		out << "<variable name=\"" << name << "\">" << val << "</variable>";
	}

	/** 
	 * \brief Prints a variable to output stream.
	 * The maximum will be printed as "MAX" string.
	 * 
	 * @param val Value.
	 * @param name Name of variable.
	 * @param out Output stream.
	 */
	template<class TT>
	void variableToXMLWithMax(TT &val, std::string name, std::ostream &out) {
		if (val == std::numeric_limits<TT>::max())
			out << "<variable name=\"" << name << "\">MAX</variable>";
		else
			out << "<variable name=\"" << name << "\">" << val << "</variable>";
	}

	/// Vector of multi dimensional variables IDs (prefix tree).
	std::vector<std::vector<uli_t> > varID;

	/// Vector of multi dimensional values (prefix tree).
	std::vector<std::vector<std::vector<T> > > values;

	/// Vector of multi dimensional branches (prefix tree).
	std::vector<std::vector<uli_t> > branches;

	/// Vector of all classes.
	std::vector<T> classes;

	/// Vector of all indexes.
	std::vector<uli_t> indexes;

	/// Vector of multi doubles (e.g. for term node infos).
	std::vector<std::vector<double> > extra;

	/// Size of multi dimensional variables
	uli_t varSize;

	/// Size of multi dimensional values
	uli_t valueSize;

	/// Width of a value
	uli_t valueWidth;

	/// Number of branches per node
	uli_t branchSize;

	/// ID of a terminal node
	uli_t termID;
};

/** 
 * \brief Prints compiled tree to output stream.
 * 
 * @param os Output stream.
 * @param cmpldTree A compiled tree.
 * 
 * @return Output stream.
 */
template<class T>
inline std::ostream & operator<<(std::ostream &os, CmpldTree<T> &cmpldTree) {
	return cmpldTree.print(os);
}

#endif /*CMPLDTREE_H_*/
