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

#ifndef INODE_H_
#define INODE_H_

/*
 * Includes
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>
#include "ClassAtom.h"
#include "Node.h"
#include "treedefs.h"

/*
 * Declaration and Definition
 */

/*
 * Class INode:
 * - name comes from "output is _int_"-Node
 */

template<class T>
class INode: public Node<T, uli_t> {
public:
	INode() :
		Node<T, uli_t> (), classifier(new ClassAtom<T, uli_t> ()) {
	}
	;

	INode(ClassAtom<T, uli_t> *classifier, uli_t outcomeEmergency =
			std::numeric_limits<uli_t>::max()) :
		Node<T, uli_t> (), classifier(classifier), outcomeEmergency(
				outcomeEmergency) {
	}
	;

	// copy constructor
	INode(const INode<T> & iNode) :
		Node<T, uli_t> (iNode), classifier(iNode.classifier) {
	}
	;

	virtual ~INode() {
		if (this->classifier != NULL)
			delete this->classifier;
	}

	ClassAtom<T, uli_t> *getClassifier() const {
		return this->classifier;
	}

	void setClassifier(ClassAtom<T, uli_t> *classifier) {
		return this->classifier = classifier;
	}

	virtual T classify(const std::vector<T> &sample) const {
		if (this->classifier->getType() == TERMINALNODE) {
			return this->classifier->getResult();
		}

		uli_t cla = this->classifier->classify(sample);

		// if sample can not be classified
		if (cla == this->classifier->getErrorMissing())
			return this->classifier->getMissingCode();

		// if classification is erroneous
		if (cla + 1 > this->kids.size()) {
			throw Exception(ERRORCODE_14);
		}

		if (this->kids[cla] == NULL) {
			throw Exception(ERRORCODE_15);
			// what should be returned here?!
			// return cla;
			// return std::numeric_limits<uli_t >::infinity();
			// I recomment to choose the most frequent outcome value as
			// the outcomeEmergency value.
			// return outcomeEmergency;
		}

		return this->kids[cla]->classify(sample);
	}

	virtual bool isThereVarID(uli_t varID) const {
		uli_t i;

		if (this->classifier->isThereVarID(varID))
			return true;

		for (i = 0; i < this->kids.size(); ++i) {
			if (i < this->kids.size() && this->kids.at(i) != NULL) {
				if (this->kids.at(i)->isThereVarID(varID))
					return true;
			}
		}
		return false;
	}

	// if kids smaller than
	virtual std::ostream &print(std::ostream &os) const {
		uli_t i, claSize;

		os << "Node";
		os << "<cla:" << *(this->classifier) << ">";
		os << "<size:" << this->kids.size() << ">";

		claSize = this->kids.size();

		if (claSize > 0)
			os << "(";
		for (i = 0; i < claSize; ++i) {
			if (i >= this->kids.size()) {
				os << "NULL";
			} else {
				if (this->kids.at(i) == NULL) {
					os << "NULL";
				} else {
					this->kids.at(i)->print(os);
				}
			}
			if (i == claSize - 1)
				os << ")";
			else
				os << ",";
		}
		return os;
	}

	// if kids smaller than
	virtual void printXml(std::ostream &os) const {
		uli_t i, claSize;

		this->classifier->printXml(os);

		claSize = this->kids.size();

		for (i = 0; i < claSize; ++i) {
			if (i >= this->kids.size()) {
				os << "<node></node>" << std::endl;
				;
			} else {
				if (this->kids.at(i) == NULL) {
					os << "<node></node>" << std::endl;
					;
				} else {
					os << "<node id=\"" << i << "\" " << "size=\""
							<< this->kids.at(i)->size() << "\">" << std::endl;
					this->kids.at(i)->printXml(os);
					os << "</node>" << std::endl;
				}
			}
		}
	}

	virtual uli_t getNumOfLeafs() {
		typename std::vector<Node<T, uli_t> *>::const_iterator it(
				this->kids.begin());
		uli_t num = 0;

		if (this->classifier->getType() == TERMINALNODE)
			return 1;

		while (it != this->kids.end()) {
			if (*it != NULL) {
				num += (*it)->getNumOfLeafs();
			} else {
				num++;
			}

			++it;
		}
		return num;
	}

	virtual double getLeafSum() {
		typename std::vector<Node<T, uli_t> *>::const_iterator it(
				this->kids.begin());
		double num = 0.0;

		if (this->classifier->getType() == TERMINALNODE)
			return (double) this->classifier->getResult();

		while (it != this->kids.end()) {
			if (*it != NULL) {
				num += (*it)->getLeafSum();
			}
			++it;
		}
		return num;
	}

	virtual uli_t getMaxDepth() {
		return 0;
	}

	virtual uli_t getNumOfNodes() {
		typename std::vector<Node<T, uli_t> *>::const_iterator it(
				this->kids.begin());
		uli_t num = 1;

		if (this->classifier->getType() == TERMINALNODE)
			return 1;

		while (it != this->kids.end()) {
			if (*it != NULL) {
				num += (*it)->getNumOfNodes();
			} else {
				num++;
			}
			++it;
		}
		return num;
	}

	ClassAtom<T, uli_t> *classifier;
	uli_t outcomeEmergency;
	//double infoGain;
};

#endif /*INODE_H_*/
