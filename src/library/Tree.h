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

#ifndef TREE_H_
#define TREE_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include "Helper.h"
#include "Node.h"


#ifndef NULL
#define NULL 0
#endif


/*
 * Source (Def. & Decl.)
 */

template <class T, class C >
class Tree {
public:
	Tree() : root(NULL) {};
	Tree(Node<T, C> *root) : root(root) {};
	virtual ~Tree() {
		if (this->root != NULL) delete this->root;
	};

	inline void setRoot(Node<T, C> *root) {this->root = root;};
	inline Node<T, C> *getRoot() const {return this->root;};


	virtual T classify(const std::vector<T > &sample) const {
		return this->root->classify(sample);
	}

	virtual std::ostream& print(std::ostream& os) const {
		this->root->print(os);
		return os;
	}

	virtual bool isThereVarID(uli_t varID) const {
    return this->root->isThereVarID(varID);
	}

	virtual void printXml(std::ostream& os) const {
    os
    << "<node id=\"0\" "
    << "size=\"" << this->root->size() << "\">" << std::endl;
		this->root->printXml(os);
		os << "</node>" << std::endl;
	}

	virtual void summary() const {
		std::ostream& os = std::cout;

		os	<< "Number of leafs: "
			<< this->root->getNumOfLeafs() << std::endl;

		os	<< "Number of nodes: "
			<< this->root->getNumOfNodes() << std::endl;
        }

	Node<T, C> *root;
};

template <class T, class C>
std::ostream &operator<<(std::ostream &os, Tree<T, C > &tree) {
	return tree.print(os);
}

#endif /*TREE_H_*/
