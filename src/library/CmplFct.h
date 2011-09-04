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

#ifndef CMPLFCT_H_
#define CMPLFCT_H_

/*
 * Includes
 */
#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <boost/dynamic_bitset.hpp>

#include "ErrorCodes.h"
#include "CmpldTree.h"
#include "INode.h"
#include "T2ClassAtom.h"
#include "SClassAtom.h"
#include "LotusTermClassAtom.h"
#include "Tree.h"
#include "Helper.h"
#include "Exception.h"
#include "DataFrame.h"
#include "Proximities.h"
#include "lr.h"

/** 
 * \brief Representation of a set of compiler functions.
 * This binder serves serveral functions to different tree types
 * as follows:
 *  - compiling tree that consits of pointers
 *  - classify samples
 *  - get cutoffs of the tree
 *  - add proximities
 *
 */
class CmplFct {
public:
	/** 
	 * \brief Constructor of CmplFct
	 */
	CmplFct();

	/** 
	 * \brief Destructor of CmplFct
	 */
	virtual ~CmplFct();

	/** 
	 * \brief Compile a T2 tree (see T2ClassAtom class)
	 * 
	 * @param tree Pointer based tree.
	 * 
	 * @return flat data structure which represents the tree
	 */
	template<class T>
	static CmpldTree<T> *T2(Tree<T, uli_t> &tree) {

		CmpldTree<T> *cmpldTree = new CmpldTree<T> (1, 1, 1, 2); //CART

		// compile tree
		makeT2((INode<T>*) tree.root, cmpldTree);

		return cmpldTree;
	}

	/** 
	 * \brief Recursive wrapper function for compiling
	 * flat trees.
	 * 
	 * @param node Node which shoulb be compiled.
	 * @param cmpldTree Container which contains resulting tree.
	 */
	template<class T>
	static void makeT2(INode<T> *node, CmpldTree<T> *cmpldTree) {

		// add node
		cmpldTree->pushBackNode();

		uli_t at = cmpldTree->size() - 1;

		// is it a term node?
		if (node->classifier->getType() == TERMINALNODE) {
			cmpldTree->toLastValues(0, node->classifier->getResult());
			return;
		}

		// add variable ID
		cmpldTree->toLastVarID(((T2ClassAtom<T>*) node->classifier)->getVarID());

		// add threshold
		cmpldTree->toLastValues(0, ((T2ClassAtom<T>*) node->classifier)->threshold);

		// add kids
		if (node->kids.size() == 0)
			throw Exception("ERROR_CODE_27");
		for (uli_t i = 0; i < node->kids.size(); ++i) {
			cmpldTree->toBranches(at, i, cmpldTree->size());
			if (node->kids.at(i) == NULL)
				throw Exception(ERRORCODE_26);
			makeT2((INode<T>*) node->kids.at(i), cmpldTree);
		}

	}

	/** 
	 * \brief Classify a sample utilizing a compiled tree.
	 * 
	 * @param sample Sample that should be classified.
	 * @param cmpldTree Classifier that should be used.
	 * 
	 * @return Terminal node value.
	 */
	template<class T>
	static T classifyT2(T *sample, CmpldTree<T> *cmpldTree) {

		uli_t at;

		// search down the nodes
		at = 0;
		while (true) {
			// if terminal node break
			if (cmpldTree->branches[at][0] == 0)
				break;

			// jump to correct branch
			if (rjungle::magicAt(sample, cmpldTree->varID[0][at])
					<= cmpldTree->values[at][0][0]) {
				at = cmpldTree->branches[at][0];
			} else {
				at = cmpldTree->branches[at][1];
			}
		}
		return cmpldTree->values[at][0][0]; // return terminal node value
	}

	/** 
	 * \brief Classify a sample utilizing a compiled tree.
	 * 
	 * @param sample Sample that should be classified.
	 * @param cmpldTree Classifier that should be used.
	 * 
	 * @return Node ID of terminal node.
	 */
	template<class T>
	static unsigned int classifyT2nodeID(T *sample, CmpldTree<T> *cmpldTree) {

		uli_t at;

		// search down the nodes
		at = 0;
		while (true) {
			// if terminal node break
			if (cmpldTree->branches[at][0] == 0)
				break;

			// jump to correct branch
			if (rjungle::magicAt(sample, cmpldTree->varID[0][at])
					<= cmpldTree->values[at][0][0]) {
				at = cmpldTree->branches[at][0];
			} else {
				at = cmpldTree->branches[at][1];
			}
		}

		return at; // node ID
	}

	/** 
	 * \brief Get cutoffs of a specific variable in tree
	 *
	 * @param cmpldTree A compiled tree.
	 * @param var Index of variable.
	 * @param cutoffs Cutoffs of a variable in tree. (Return value).
	 */
	template<class T>
	static void getCutoffsT2(CmpldTree<T> &cmpldTree, uli_t var,
			std::vector<T> &cutoffs) {
		// map of cutoffs
		std::map<T, int> cutoffMap;
		typename std::map<T, int>::iterator iter;

		// free vector
		cutoffs.clear();

		uli_t at;

		// search down the nodes

		for (at = 0; at < cmpldTree.branches.size(); ++at) {
			// if terminal node break
			if (cmpldTree.branches[at][0] == 0)
				continue;

			// jump to correct branch
			if (cmpldTree.varID[0][at] == var) {
				cutoffMap[cmpldTree.values[at][0][0]] = 1;
			}
		}

		// get keys (= cutoffs)
		for (iter = cutoffMap.begin(); iter != cutoffMap.end(); ++iter) {
			cutoffs.push_back(iter->first);
		}
	}

	// return 2D proximity probability matrix

	/** 
	 * \brief Add sample proximities of current tree to proximity matrix.
	 * 
	 * @param data Data set which contains sample data.
	 * @param proxi Sample proximity matrix.
	 * @param cmpldTree Compiled tree.
	 * @param colMask Mask of pre selected variables.
	 */
	template<class T>
	static void addVarProxT2(DataFrame<T> &data, Proximities<T> &proxi,
			CmpldTree<T> *cmpldTree, std::vector<uli_t> *colMask) {

		uli_t i, at, parent, kid, numOfNodePairs;
		std::vector<uli_t> atStack;
		std::map<uli_t, uli_t> varIDMap;
		Proximities<T> proxiLocal(proxi.nrow, proxi.ncol);
		proxiLocal.reset();

		// create index map
		if (colMask == NULL) {
			for (i = 0; i < data.par.ncol; ++i) {
				varIDMap[i] = i;
			}
		} else {
			for (i = 0; i < colMask->size(); ++i) {
				varIDMap[colMask->at(i)] = i;
			}
		}

		// search down the nodes
		at = parent = numOfNodePairs = 0;
		while (at < cmpldTree->branches.size()) {
			/*
			 std::cout << at << ">" << cmpldTree->branches[at][0] << ","
			 << cmpldTree->branches[at][1] << " " << "(" << atStack.size() << ")";
			 if (atStack.size() != 0)
			 std::cout << atStack.back();
			 std::cout << " max:" << cmpldTree->branches.size() << std::endl;
			 */
			// if terminal node break
			if (cmpldTree->branches[at][0] == 0) {
				if (atStack.size() == 0)
					break;

				// go down the last right branch
				at = atStack.back();
				atStack.pop_back();

				parent = varIDMap[cmpldTree->varID[0][at]];
				at = cmpldTree->branches[at][1];
			} else {

				// save varID of parent
				kid = varIDMap[cmpldTree->varID[0][at]];

				// do not count root
				if (at > 0) {
					++numOfNodePairs;
					proxiLocal.add(kid, parent, 1);
					proxiLocal.add(parent, kid, 1);
				}

				parent = kid;

				// save node
				atStack.push_back(at);

				// walk down the node (left)
				at = cmpldTree->branches[at][0];
			}
		}

		// normalize
		proxiLocal.div(numOfNodePairs);

		// save proximity data
		proxi.add(proxiLocal);
	}

	////////////////////////////////////////////////////////////////////////////

	/** 
	 * \brief Compile a S tree (see SClassAtom class)
	 * 
	 * @param tree Pointer based tree.
	 * 
	 * @return flat data structure which represents the tree
	 */
	template<class T>
	static CmpldTree<T> *S(Tree<T, uli_t> &tree) {

		//Nominal (unknown value size)
		CmpldTree<T> *cmpldTree = new CmpldTree<T> (1, 1, 0, 2);

		// compile tree
		makeS((INode<T>*) tree.root, cmpldTree);

		return cmpldTree;
	}

	/** 
	 * \brief Recursive wrapper function for compiling
	 * flat trees.
	 * 
	 * @param node Node which shoulb be compiled.
	 * @param cmpldTree Container which contains resulting tree.
	 */
	template<class T>
	static void makeS(INode<T> *node, CmpldTree<T> *cmpldTree) {

		uli_t i;

		// add node
		cmpldTree->pushBackNode();

		uli_t at = cmpldTree->size() - 1;

		// is it a term node?
		if (node->classifier->getType() == TERMINALNODE) {
			// delete next line?
			cmpldTree->pushBackToLastValues(node->classifier->getResult());
			return;
		}

		// add variable ID
		cmpldTree->toLastVarID(((SClassAtom<T>*) node->classifier)->getVarID());

		// add elements in set
		std::vector<T> &classSet = ((SClassAtom<T>*) node->classifier)->classSet;
		for (i = 0; i < classSet.size(); ++i) {
			cmpldTree->pushBackToLastValues(classSet[i]);
		}

		// add kids
		if (node->kids.size() == 0)
			throw Exception(ERRORCODE_27);
		for (i = 0; i < node->kids.size(); ++i) {
			cmpldTree->toBranches(at, i, cmpldTree->size());
			if (node->kids.at(i) == NULL)
				throw Exception(ERRORCODE_26);
			makeS((INode<T>*) node->kids.at(i), cmpldTree);
		}

	}

	/** 
	 * \brief Classify a sample utilizing a compiled tree.
	 * 
	 * @param sample Sample that should be classified.
	 * @param cmpldTree Classifier that should be used.
	 * 
	 * @return Terminal node value.
	 */
	template<class T>
	static T classifyS(T *sample, CmpldTree<T> *cmpldTree) {

		uli_t at, atPreFetched, i;

		// search down the nodes
		at = 0;
		while (true) {
			// if terminal node break
			if (cmpldTree->branches[at][0] == 0)
				break;

			// jump to correct branch
			atPreFetched = cmpldTree->branches[at][1];

			for (i = 0; i < cmpldTree->values[at][0].size(); ++i) {
				if (rjungle::magicAt(sample, cmpldTree->varID[0][at])
						== cmpldTree->values[at][0][i]) {
					atPreFetched = cmpldTree->branches[at][0];
					break;
				}
			}

			at = atPreFetched;
		}

		return cmpldTree->values[at][0][0]; // return terminal node value
	}

	/** 
	 * \brief Classify a sample utilizing a compiled tree.
	 * 
	 * @param sample Sample that should be classified.
	 * @param cmpldTree Classifier that should be used.
	 * 
	 * @return Node ID of terminal node.
	 */
	template<class T>
	static unsigned int classifySnodeID(T *sample, CmpldTree<T> *cmpldTree) {

		uli_t at, atPreFetched, i;

		// search down the nodes
		at = 0;
		while (true) {
			// if terminal node break
			if (cmpldTree->branches[at][0] == 0)
				break;

			// jump to correct branch
			atPreFetched = cmpldTree->branches[at][1];

			for (i = 0; i < cmpldTree->values[at][0].size(); ++i) {
				if (rjungle::magicAt(sample, cmpldTree->varID[0][at])
						== cmpldTree->values[at][0][i]) {
					atPreFetched = cmpldTree->branches[at][0];
					break;
				}
			}

			at = atPreFetched;
		}

		return at; // node ID
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	/** 
	 * \brief Compile a Lotus tree (see LotusTermClassAtom class)
	 * 
	 * @param tree Pointer based tree.
	 * 
	 * @return flat data structure which represents the tree
	 */
	template<class T>
	static CmpldTree<T> *Lotus(Tree<T, uli_t> &tree) {

		CmpldTree<T> *cmpldTree = new CmpldTree<T> (1, 1, 1, 2); //CART

		// compile tree
		makeLotus((INode<T>*) tree.root, cmpldTree);

		return cmpldTree;
	}

	/** 
	 * \brief Recursive wrapper function for compiling
	 * flat trees.
	 * 
	 * @param node Node which shoulb be compiled.
	 * @param cmpldTree Container which contains resulting tree.
	 */
	template<class T>
	static void makeLotus(INode<T> *node, CmpldTree<T> *cmpldTree) {

		// add node
		cmpldTree->pushBackNode();

		uli_t at = cmpldTree->size() - 1;
		uli_t atExtra = 0;

		// is it a term node?
		if (node->classifier->getType() == TERMINALNODE) {
			// add extra element for logistic regression parameters
			cmpldTree->extra.push_back(std::vector<double>());
			atExtra = cmpldTree->extra.size() - 1;

			// add parameters
			LotusTermClassAtom<T, uli_t>* ltca =
					(LotusTermClassAtom<T, uli_t>*) node->classifier;

			cmpldTree->toLastVarID(ltca->varID);
			cmpldTree->extra[atExtra].push_back(ltca->betas[0]);
			cmpldTree->extra[atExtra].push_back(ltca->betas[1]);

			// add index of corresponding extra vector to branches 
			cmpldTree->toLastBranches(1, atExtra);

			return;
		}

		// add variable ID
		cmpldTree->toLastVarID(((T2ClassAtom<T>*) node->classifier)->getVarID());

		// add threshold
		cmpldTree->toLastValues(0, ((T2ClassAtom<T>*) node->classifier)->threshold);

		// add kids
		if (node->kids.size() == 0)
			throw Exception("ERROR_CODE_27");
		for (uli_t i = 0; i < node->kids.size(); ++i) {
			cmpldTree->toBranches(at, i, cmpldTree->size());
			if (node->kids.at(i) == NULL)
				throw Exception(ERRORCODE_26);
			makeLotus((INode<T>*) node->kids.at(i), cmpldTree);
		}

	}

	/** 
	 * \brief Classify a sample utilizing a compiled tree.
	 * 
	 * @param sample Sample that should be classified.
	 * @param cmpldTree Classifier that should be used.
	 * 
	 * @return Terminal node value.
	 */
	template<class T>
	static T predictLotus(T *sample, CmpldTree<T> *cmpldTree) {

		uli_t at;

		// search down the nodes
		at = 0;
		while (true) {
			// if terminal node break
			if (cmpldTree->branches[at][0] == 0)
				break;

			// jump to correct branch
			if (rjungle::magicAt(sample, cmpldTree->varID[0][at])
					<= cmpldTree->values[at][0][0]) {
				at = cmpldTree->branches[at][0];
			} else {
				at = cmpldTree->branches[at][1];
			}
		}

		// get lr parameters

		uli_t varID = cmpldTree->varID[0][at]; // get index of variable
		uli_t atExtra = cmpldTree->branches[at][1]; // get index of extra vector
		double xval = (double) rjungle::magicAt(sample, varID); // get x value
		double b0 = (cmpldTree->extra[atExtra]).at(0); // get beta0
		double b1 = (cmpldTree->extra[atExtra]).at(1); // get beta1
		double en = exp(b0 + b1 * xval); // calc. logit prediction
		double result = en / (1.0 + en); // Compute model value, u = e^n / (1.0 + e^n) 

		// avoid 1 or 0
		if (result == 1)
			result = 1 - std::numeric_limits<double>::min();

		if (result == 0)
			result = std::numeric_limits<double>::min();

		return result;
	}

};

#endif /*CMPLFCT_H_*/
