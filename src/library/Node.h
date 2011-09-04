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

#ifndef NODE_H_
#define NODE_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include "ClassAtom.h"

#ifndef TCNodeVector
#define TCNodeVector std::vector<Node<T, C > *>
#endif


template <class T, class C>
class Node {
public:
	Node() : kids(TCNodeVector()) {  };

	// copy constructor
	Node(const Node<T, C > &node ) : TCNodeVector(node) { };

	virtual ~Node() {
		typename TCNodeVector::const_iterator it(this->kids.begin());

		while(it != this->kids.end()) {
			if (*it != NULL) delete *it;
			++it;
		}

	}

	inline const Node<T, C > *getParent() const { return this->parent; }
	inline void setParent(Node<T, C > * parent) { this->parent = parent; }

	inline Node<T, C > *operator[](uli_t i) const { return this->kids[i]; }
	inline Node<T, C > *operator[](uli_t i) { return this->kids[i]; } //pray to god, that current content is NULL

	inline Node<T, C > *at(uli_t i) const { return this->kids[i]; }
	inline Node<T, C > *at(uli_t i) { return this->kids[i]; } //pray to god, that current content is NULL


	inline void resize(uli_t size) { this->kids.resize(size); }

	inline void push_back(Node<T, C >* node) { this->kids.push_back(node); }

	inline uli_t size() { return this->kids.size(); }

	virtual uli_t getNumOfLeafs() { return 0; }
	virtual double getLeafSum() { return 0.0; }
	virtual uli_t getMaxDepth() { return 0; }
	virtual uli_t getNumOfNodes() { return 0; }

	/*
	 * classify recursivly
	 */
	virtual C getPath(const T& val) const {
		return C();
	}

	virtual T classify(const std::vector<T > &sample) const {
		return T();
	}

	virtual std::ostream& print(std::ostream& os) const {
		return os;
	}

	virtual bool isThereVarID(uli_t varID) const { return false; }

  virtual void printXml(std::ostream &os) const {
  }
	// ID of the corresponding Variable
	Node<T, C > *parent;
	TCNodeVector kids;
};

template <class T, class C>
std::ostream &operator<<(std::ostream &os, Node<T, C > &node) {
	return node.print(os);
}

#endif /*NODE_H_*/
