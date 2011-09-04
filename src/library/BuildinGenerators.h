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

#ifndef BUILDINGENERATORS_H_
#define BUILDINGENERATORS_H_

/*
 * Includes
 */

#include "treedefs.h"
#include "RJungleHelper.h"
#include "TImportance.h"
#include "FittingFct.h"

/** 
 * \brief A binder for all tree specific functions.
 * This is kind of configuration that will be used to
 * handle several functionalities.
 */
template<class T>
class BuildinGenFct {
public:

	/** 
	 * \brief Constructor the resets all pointer to NULL.
	 */
	BuildinGenFct() :
		newImportanceObject(NULL), getBestVariable(NULL), getNodePerformance(NULL),
				getTermResult(NULL), cmplTree(NULL), classifyCmpldTree(NULL),
				classifyIDcmpldTree(NULL), getAccuracyOfCmpldTrees(NULL),
				checkUsability(NULL), getCutoffs(NULL) {
	}
	;

	/** 
	 * \brief Destructor that does nothing.
	 */
	virtual ~BuildinGenFct() {
	}
	;

	/// Importance score structure and functions.
	Importance<T>* (*newImportanceObject)(void);

	/// Fitting function that determines the "best" variable in data set.
	ClassAtom<T, uli_t>* (*getBestVariable)(FittingFctPar<T> &fitFctPar);

	/// Function that calculates the performance of the current node.
	double (*getNodePerformance)(DataFrame<T>&data, uli_t depVar, std::vector<
			uli_t>&rowMaskVec, std::vector<double> &rowWeight);

	/// Fitting function that determines the "best" variable in data set.
	//TermClassAtom<T, uli_t>* (*getBestVariable)(FittingFctPar<T> &fitFctPar);

	/// Function that determines the classifier of current terminal node.
	TermClassAtom<T, uli_t> *(*getTermResult)(DataFrame<T> &data,
			RJunglePar &par, uli_t col, std::vector<uli_t> &rowMaskVec, std::vector<
					uli_t> &colMaskVec, std::vector<double> &rowWeight);

	/// Function that compiles a pointer based tree to a flat data structure.
	CmpldTree<T>* (*cmplTree)(Tree<T, uli_t> &tree);

	/// Classifier function that classifies a sample using a compiled tree.
	T (*classifyCmpldTree)(T *sample, CmpldTree<T> *cmpldTree);

	/// Classifier function that classifies a sample using a compiled tree and returns the id of the result.
	unsigned int (*classifyIDcmpldTree)(T *sample, CmpldTree<T> *cmpldTree);

	/// Function that returns the accuracy of current tree(s).
	double (*getAccuracyOfCmpldTrees)(DataFrame<T> *data, T* pred, uli_t depVar,
			DataTreeSet &dataSet, long int treeId, bool one, bool showInfos,
			RJungleIO &io);

	/// Pre function that checks before tree growing if data is suitable for specific tree type.
	void (*checkUsability)(DataFrame<T> &data);

	/// Function that returns the center of a vector such as mean or median.
	double (*getCenter)(std::vector<T> &vec);

	/// Function that updates the proximity matrix.
	void (*addVarProx)(DataFrame<T>&data, Proximities<T> &proxi,
			CmpldTree<T> *cmpldTree, std::vector<uli_t> *colMask);

	/// Function that calculates the default mtry value.
	uli_t (*setMtry)(uli_t ncol);

	/// Function that returns all cutoffs within a tree. (E.g. for conditional importance).
	void (*getCutoffs)(CmpldTree<T> &cmpldTree, uli_t var,
			std::vector<T> &cutoffs);
};

/*
 * Includes
 */

#include "RJungleAcc.h"
#include "TermResult.h"

/*
 * * Buildin generator functions
 * for all trees
 *
 */

namespace rjungle {
namespace Generator {

/*
 * CART
 */

namespace CART {
template<class T>

/**
 * \brief Assign all CART (y=nominal/x=numeric) specific functions.
 *
 * @param genFct Generator object to be set.
 */
void assignFunctions(BuildinGenFct<T> *genFct) {
	genFct->getBestVariable = FittingFct::CARTcateCont<T>;
	genFct->getNodePerformance = Helper::purityGini<T>;
	genFct->getTermResult = TermResult::getTermMostFreq<T>;
	genFct->getAccuracyOfCmpldTrees = RJungleAcc<T>::getAccuracyCmpld;
	genFct->cmplTree = CmplFct::T2<T>;
	genFct->classifyCmpldTree = CmplFct::classifyT2<T>;
	genFct->classifyIDcmpldTree = CmplFct::classifyT2nodeID<T>;
	genFct->newImportanceObject = TImportance<T>::newImportanceObject;
	genFct->checkUsability = NULL;
	genFct->getCenter = Helper::getMedian<T>;
	genFct->addVarProx = CmplFct::addVarProxT2<T>;
	genFct->setMtry = RJungleHelper<T>::setMtryClassification;
	genFct->getCutoffs = CmplFct::getCutoffsT2<T>;
}
;
}

namespace CARTcateCate {
template<class T>

/**
 * \brief Assign all CART (y=nominal/x=nominal) specific functions.
 *
 * @param genFct Generator object to be set.
 */
void assignFunctions(BuildinGenFct<T> *genFct) {
	genFct->getBestVariable = FittingFct::CARTcateCate<T>;
	genFct->getNodePerformance = Helper::purityNominal<T>;
	genFct->getTermResult = TermResult::getTermMostFreq<T>;
	genFct->getAccuracyOfCmpldTrees = RJungleAcc<T>::getAccuracyCmpld;
	genFct->cmplTree = CmplFct::S<T>;
	genFct->classifyCmpldTree = CmplFct::classifyS<T>;
	genFct->classifyIDcmpldTree = CmplFct::classifySnodeID<T>;
	genFct->newImportanceObject = TImportance<T>::newImportanceObject;
	genFct->checkUsability = NULL;
	genFct->getCenter = Helper::getMode<T>;
	genFct->addVarProx = NULL;
	genFct->setMtry = RJungleHelper<T>::setMtryClassification;
	genFct->getCutoffs = NULL;
}
;
}

namespace CARTcontCont {
template<class T>

/**
 * \brief Assign all CART (y=numeric/x=numeric) specific functions.
 *
 * @param genFct Generator object to be set.
 */
void assignFunctions(BuildinGenFct<T> *genFct) {
	genFct->getBestVariable = FittingFct::CARTcontCont<T>;
	genFct->getNodePerformance = Helper::puritySos<T>;
	genFct->getTermResult = TermResult::getTermMean<T>;
	genFct->getAccuracyOfCmpldTrees = RJungleAcc<T>::getAccuracyCmpldSos;
	genFct->cmplTree = CmplFct::T2<T>;
	genFct->classifyCmpldTree = CmplFct::classifyT2<T>;
	genFct->classifyIDcmpldTree = CmplFct::classifyT2nodeID<T>;
	genFct->newImportanceObject = TImportance<T>::newImportanceObject;
	genFct->checkUsability = NULL;
	genFct->getCenter = Helper::getMedian<T>;
	genFct->addVarProx = CmplFct::addVarProxT2<T>;
	genFct->setMtry = RJungleHelper<T>::setMtryRegression;
	genFct->getCutoffs = CmplFct::getCutoffsT2<T>;
}
;
}

namespace CARTcontCate {
template<class T>

/**
 * \brief Assign all CART (y=numeric/x=nominal) specific functions.
 *
 * @param genFct Generator object to be set.
 */
void assignFunctions(BuildinGenFct<T> *genFct) {
	genFct->getBestVariable = FittingFct::CARTcontCate<T>;
	genFct->getNodePerformance = Helper::puritySos<T>;
	genFct->getTermResult = TermResult::getTermMean<T>;
	genFct->getAccuracyOfCmpldTrees = RJungleAcc<T>::getAccuracyCmpldSos;
	genFct->cmplTree = CmplFct::S<T>;
	genFct->classifyCmpldTree = CmplFct::classifyS<T>;
	genFct->classifyIDcmpldTree = CmplFct::classifySnodeID<T>;
	genFct->newImportanceObject = TImportance<T>::newImportanceObject;
	genFct->checkUsability = NULL;
	genFct->getCenter = Helper::getMode<T>;
	genFct->addVarProx = NULL;
	genFct->setMtry = RJungleHelper<T>::setMtryRegression;
	genFct->getCutoffs = NULL;
}
;
}

/*
 * CART extensions
 */

/**
 * \brief Assign all CART (y=nominal/x=numeric) specific functions.
 *
 * @param genFct Generator object to be set.
 */
namespace CARTsrt {
template<class T>
void assignFunctions(BuildinGenFct<T> *genFct) {
	genFct->getBestVariable = FittingFct::CARTcateContSrt<T>;
	genFct->getNodePerformance = Helper::purityGini<T>;
	genFct->getTermResult = TermResult::getTermMostFreq<T>;
	genFct->getAccuracyOfCmpldTrees = RJungleAcc<T>::getAccuracyCmpld;
	genFct->cmplTree = CmplFct::T2<T>;
	genFct->classifyCmpldTree = CmplFct::classifyT2<T>;
	genFct->classifyIDcmpldTree = CmplFct::classifyT2nodeID<T>;
	genFct->newImportanceObject = TImportance<T>::newImportanceObject;
	genFct->checkUsability = NULL;
	genFct->getCenter = Helper::getMedian<T>;
	genFct->addVarProx = NULL;
	genFct->setMtry = RJungleHelper<T>::setMtryClassification;
	genFct->getCutoffs = CmplFct::getCutoffsT2<T>;
}
;
}

/*
 * Experimental
 */

/**
 * \brief Experimental fitting function. Implements LOTUS trees.
 *
 * @param genFct Generator object to be set.
 */
namespace LOTUS {
template<class T>
void assignFunctions(BuildinGenFct<T> *genFct) {
	genFct->getBestVariable = FittingFct::LOTUS<T>;
	genFct->getNodePerformance = Helper::purityLotus<T>;
	genFct->getTermResult = TermResult::getTermLotus<T>;
	genFct->getAccuracyOfCmpldTrees = RJungleAcc<T>::getDevianceCmpld;
	genFct->cmplTree = CmplFct::Lotus<T>;
	genFct->classifyCmpldTree = CmplFct::predictLotus<T>;
	genFct->classifyIDcmpldTree = CmplFct::classifyT2nodeID<T>;
	genFct->newImportanceObject = TImportance<T>::newImportanceObject;
	genFct->checkUsability = FittingFct::checkUsabilityLOTUS<T>;
	genFct->getCenter = Helper::getMedian<T>;
	genFct->addVarProx = NULL;
	genFct->setMtry = RJungleHelper<T>::setMtryClassification;
	genFct->getCutoffs = CmplFct::getCutoffsT2<T>;
}
;
}
}
}
#endif /* BUILDINGENERATORS_H */

