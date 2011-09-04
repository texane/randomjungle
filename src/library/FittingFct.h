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

#ifndef FITTINGFCT_H_
#define FITTINGFCT_H_

/*
 * Includes
 */
#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <boost/dynamic_bitset.hpp>

#include "DataFrame.h"
#include "ClassAtom.h"
#include "SClassAtom.h"
#include "T2ClassAtom.h"
#include "TMClassAtom.h"
#include "TermClassAtom.h"
#include "TImportance.h"
#include "Tree.h"
#include "Helper.h"
#include "INode.h"

#include "lr.h"

/*
 * Source (Def. & Decl.)
 */

template<class T>
class FittingFctPar {
public:
	FittingFctPar() :
		parent(NULL), data(NULL), depVar(0), rowMaskVec(NULL), colMaskVec(NULL),
				rowWeight(NULL), nodePerformance(0), iteration(0), depth(0), rng(NULL),
				treeImprovement(NULL), bestVar(NULL), stopGrowing(false),
				twoWayInteraction(NULL), io(NULL) {
	}

	FittingFctPar(
			INode<T> *parent,
			DataFrame<T> *data,
			uli_t depVar,
			std::vector<uli_t> *rowMaskVec,
			std::vector<uli_t> *colMaskVec,
			std::vector<double> *rowWeight,
			double nodePerformance,
			uli_t iteration,
			uli_t depth,
			gsl_rng *rng,
			//OUTPUT
			Importance<T> *importance,
			double *treeImprovement,
			uli_t *bestVar,
			bool stopGrowing,
			std::vector<std::pair<double, std::pair<uli_t, uli_t> > > *twoWayInteraction,
			RJungleIO *io) :
		parent(parent), data(data), depVar(depVar), rowMaskVec(rowMaskVec),
				colMaskVec(colMaskVec), rowWeight(rowWeight), nodePerformance(
						nodePerformance), iteration(iteration), depth(depth), rng(rng),
				importance(importance), treeImprovement(treeImprovement), bestVar(
						bestVar), stopGrowing(stopGrowing), twoWayInteraction(
						twoWayInteraction), io(io) {
	}

	virtual ~FittingFctPar() {
	}

	INode<T> *parent;
	DataFrame<T> *data;
	uli_t depVar;
	std::vector<uli_t> *rowMaskVec;
	std::vector<uli_t> *colMaskVec;
	std::vector<double> *rowWeight;
	double nodePerformance;
	uli_t iteration;
	uli_t depth;
	gsl_rng *rng;
	//OUTPUT
	Importance<T> *importance;
	double *treeImprovement;
	uli_t *bestVar;
	bool stopGrowing;
	std::vector<std::pair<double, std::pair<uli_t, uli_t> > > *twoWayInteraction;
	RJungleIO *io;

};

class FittingFct {
public:
	FittingFct();
	virtual ~FittingFct();

	template<class T>
	static ClassAtom<T, uli_t> *CARTcontCont(FittingFctPar<T> &fFP) {

		// declaration and definition
		uli_t k, mt; // counter(s)
		uli_t nsp; // index of sample
		uli_t nbestt, nbest; // index of best split
		uli_t non = 0; // number of variables which did have only one kind of values
		std::vector<uli_t> &rowMask = *fFP.rowMaskVec; // selection of samples
		double &decsplit = *fFP.treeImprovement; // final decrease of sum of square
		DataFrame<T> &data = *fFP.data; // input data
		uli_t resp = fFP.depVar; // index of response variable
		T d; // response variable of a sample
		uli_t smpl; // index of sample
		uli_t kv; // index of variable
		uli_t &msplit = *fFP.bestVar; // index of best variable
		double ss; // sum of square and average value of node
		double critmax = 0; // maximal sum of square
		double crit = 0; // sum of square of split
		double ubestt = 0, ubest = 0; // result values of best split
		std::vector<uli_t> &colMask = *fFP.colMaskVec; // selection of columns
		std::vector<T> yl; // values of response variable
		std::vector<T> xt, ut; // values of indep. variable
		std::vector<std::pair<T, uli_t> > validx; // pair of variable value and rank
		T2ClassAtom<T> *classifier = new T2ClassAtom<T> (); // classifier
		double suml, sumr, sumnode; // sums of new left, right and current node
		double npopl, npopr, nodecnt; // number of samples of new left, ...
		double av, critvar;
		av = ss = critvar = 0;

		if (rowMask.size() == 1) {// stop it, node size == 1
			fFP.stopGrowing = true;
			return classifier;
		}

		// get values of indep. variable
		data.getRowMaskedCol(resp, rowMask, yl);

		// get initial standard deviation and average of current node
		for (k = 0; k < rowMask.size(); ++k) {

			// get index of current sample
			smpl = rowMask[k];

			// get response value of current sample
			d = yl[k];

			// add amount to sum of square
			ss = ss + (double) k * pow(av - d, 2) / (double) (k + 1);

			// add amount to average
			av = ((double) k * av + d) / (double) (k + 1);

			// get values of response variable
			yl[k] = data.at(smpl, resp);
		}

		// number of samples in current node
		nodecnt = rowMask.size();

		// sum of all indep. variable values in current node
		sumnode = av * nodecnt;

		// find best split/variable
		for (mt = 0; mt < colMask.size(); ++mt) {

			// get index of current sample
			kv = colMask[mt];

			// reset decrease of sum of squares of current variable
			// critvar = 0;

			// get values of indep. variable
			data.getRowMaskedCol(kv, rowMask, xt);

			// paste indexes to indep. variable
			validx.clear();
			for (k = 0; k < rowMask.size(); ++k) {
				validx.push_back(std::make_pair(xt[k], k));
			}

			// sort indep. variable
			sort(validx.begin(), validx.end());

			// if there is only one kind of values in this variables then process next
			if ((validx.begin())->first == (validx.rbegin())->first) {
				++non;

				if (non > 3 * data.par.ncol) { // too many variables were pure, so stop it
					fFP.stopGrowing = true;
					return classifier;
				}

				continue;
			}

			// init sums and number of samples
			suml = npopl = 0;
			sumr = sumnode;
			npopr = nodecnt;

			// get best split of current indep. variable
			for (nsp = 0; nsp < rowMask.size() - 1; ++nsp) {

				// get response value of current sample (nsp)
				d = yl[validx[nsp].second];

				// adjust sums and sample counts
				suml += d;
				sumr -= d;
				++npopl;
				--npopr;

				if (validx[nsp].first < validx[nsp + 1].first) // value did change
					crit = (pow(suml, 2) / npopl) + (pow(sumr, 2) / npopr);

				if (crit > critvar) {// crit is the best split result til now

					// get split value
					ubestt = (validx[nsp].first + validx[nsp + 1].first) / 2.0;

					// store best value
					critvar = crit;

					// store index of sample
					nbestt = nsp;
				}
			}

			// current indep. variable had got the best split
			if (critvar > critmax) {

				// store best split value
				ubest = ubestt;

				// store best index of sample
				nbest = nbestt;

				// store index of variable
				msplit = kv;

				// store best crit value
				critmax = critvar;
			}
		}

		// store final decrease of sum of squares
		decsplit = critmax - (pow(sumnode, 2) / nodecnt);

		// build classifier
		classifier->setVarID(msplit);
		classifier->setThreshold((T) ubest);
		classifier->setMissingCode(data.getMissingCode());

		//add importance
		((TImportance<T> *) fFP.importance)->add(msplit, decsplit / data.par.ntree);

		return classifier;
	}

	template<class T>
	static ClassAtom<T, uli_t> *CARTcontCate(FittingFctPar<T> &fFP) {
		DataFrame<T> &data = *fFP.data;
		uli_t depVar = fFP.depVar;
		std::vector<uli_t> &rowMaskVec = *fFP.rowMaskVec;
		std::vector<uli_t> &colMaskVec = *fFP.colMaskVec;
		//OUTPUT
		double &decOfSos = *fFP.treeImprovement;
		uli_t &sosBestIdx = *fFP.bestVar;

		uli_t j, k;
		uli_t idxCell;
		T valCell;
		std::vector<T> colVec;
		std::vector<uli_t> setCur;
		std::vector<uli_t> bestSubSet;
		std::vector<uli_t> classVarVec;
		std::vector<double> nl, nr, nrinit;
		std::vector<uli_t> tmpclass;
		std::vector<uli_t> tclasspop;
		std::vector<uli_t> idxVec;
		std::vector<uli_t> idxRowVec;
		std::vector<uli_t> maskVec;
		std::vector<uli_t> classMaskVec;
		std::vector<uli_t> rowMaskForDepVec;
		std::vector<uli_t> classOutcomeVec;
		typename std::vector<uli_t>::iterator idxRow;
		typename std::vector<uli_t>::iterator idxCol;
		SClassAtom<T> *classifier = new SClassAtom<T> ();
		boost::dynamic_bitset<> powerSet;

		// declaration and definition
		uli_t non = 0; // number of variables which did have only one kind of values
		std::vector<uli_t> &rowMask = *fFP.rowMaskVec; // selection of samples
		uli_t resp = fFP.depVar; // index of response variable
		T d; // response variable of a sample
		uli_t smpl; // index of sample
		double ss; // sum of square and average value of node
		double critmax = 0; // maximal sum of square
		std::vector<T> yl; // values of response variable
		std::vector<T> xt, ut; // values of indep. variable
		std::vector<std::pair<T, uli_t> > validx; // pair of variable value and rank
		double suml, sumr, sumnode; // sums of new left, right and current node
		double npopl, npopr, nodecnt; // number of samples of new left, ...
		double av, critvar;
		av = ss = critvar = 0;

		T val;

		if (rowMaskVec.size() == 1) {// stop it, node size == 1
			fFP.stopGrowing = true;
			return classifier;
		}

		// get values of indep. variable
		data.getRowMaskedCol(resp, rowMask, yl);

		// get initial standard deviation and average of current node
		for (k = 0; k < rowMaskVec.size(); ++k) {

			// get index of current sample
			smpl = rowMaskVec[k];

			// get response value of current sample
			d = yl[k];

			// add amount to sum of square
			ss = ss + (double) k * pow(av - d, 2) / (double) (k + 1);

			// add amount to average
			av = ((double) k * av + d) / (double) (k + 1);

			// get values of response variable
			yl[k] = data.at(smpl, resp);
		}

		// number of samples in current node
		nodecnt = rowMaskVec.size();

		// sum of all indep. variable values in current node
		sumnode = av * nodecnt;

		idxCol = colMaskVec.begin();
		while (idxCol != colMaskVec.end()) {
			//skip dep var
			if (*idxCol == depVar) {
				++idxCol;
				continue;
			}

			// create class indexes of current row
			idxRowVec.clear();
			for (j = 0; j < rowMaskVec.size(); ++j) {
				valCell = data.at(rowMaskVec[j], *idxCol);
				idxCell = data.find(valCell, data.varClassIndexer[*idxCol]);
				idxRowVec.push_back(idxCell);
			}

			// get variable's classes
			Helper::getClasses<uli_t>(idxRowVec, classVarVec);

			// break if there is just one class
			if (classVarVec.size() == 1) {
				++non;

				if (non > 3 * data.par.ncol) { // too many variables were pure, so stop it
					fFP.stopGrowing = true;
					return classifier;
				}

				++idxCol;
				continue;
			}

			// check all combinations
			// => check all sets in power set 2^(numberofclasses-1) - 1
			// (avoiding duplicates)

			// get powerset (binary representation)
			powerSet.reset();
			powerSet.resize(classVarVec.size(), false);

			// check all sets
			do {

				// init sums and number of samples
				suml = npopl = 0;
				sumr = sumnode;
				npopr = nodecnt;

				// get the differences
				idxRow = rowMaskVec.begin();
				j = 0;
				while (idxRow != rowMaskVec.end()) {

					if (powerSet[idxRowVec[j]]) {

						val = data.at(*idxRow, depVar);

						// adjust sums and sample counts
						suml += val;
						sumr -= val;
						++npopl;
						--npopr;

					}
					++j;
					++idxRow;
				}

				critvar = (pow(suml, 2) / npopl) + (pow(sumr, 2) / npopr);

				// current indep. variable had got the best split
				if (critvar > critmax) {
					Helper::getPowerSet(powerSet, classVarVec, setCur);
					critmax = critvar;
					sosBestIdx = *idxCol;
					bestSubSet = setCur;
				}

				// get next set
				Helper::inc(powerSet);

			} while (!powerSet[powerSet.size() - 1]);

			++idxCol;
		}

		// output dec. of sos
		if (critmax == 0) {
			fFP.stopGrowing = true;
			decOfSos = 0;
		} else {
			// store final decrease of sum of squares
			decOfSos = critmax - (pow(sumnode, 2) / nodecnt);

			// get output vector
			Helper::getIndexedVector(data.varClassIndexer[sosBestIdx], bestSubSet,
					colVec);

			// build classifier
			classifier->setVarID(sosBestIdx);
			classifier->setClassSet(colVec);
			classifier->setMissingCode(data.getMissingCode());

			//add importance
			((TImportance<T> *) fFP.importance)->add(classifier->getVarID(), decOfSos
					/ data.par.ntree);

		}
		return classifier;
	}

	template<class T>
	static ClassAtom<T, uli_t> *CARTcateCate(FittingFctPar<T> &fFP) {
		DataFrame<T> &data = *fFP.data;
		uli_t depVar = fFP.depVar;
		std::vector<uli_t> &rowMaskVec = *fFP.rowMaskVec;
		std::vector<uli_t> &colMaskVec = *fFP.colMaskVec;
		std::vector<double> &rowWeight = *fFP.rowWeight;
		//OUTPUT
		double &decOfNominal = *fFP.treeImprovement;
		uli_t &nominalBestIdx = *fFP.bestVar;

		uli_t j, k;
		uli_t idxCell;
		T valCell;
		std::vector<T> colVec;
		std::vector<uli_t> setCur;
		std::vector<uli_t> bestSubSet;
		std::vector<uli_t> classVarVec;
		std::vector<double> nl, nr, nrinit;
		std::vector<uli_t> tmpclass;
		std::vector<uli_t> tclasspop;
		std::vector<uli_t> idxVec;
		std::vector<uli_t> idxRowVec;
		std::vector<uli_t> maskVec;
		std::vector<uli_t> classMaskVec;
		std::vector<uli_t> rowMaskForDepVec;
		std::vector<uli_t> classOutcomeVec;
		typename std::vector<uli_t>::iterator idxRow;
		typename std::vector<uli_t>::iterator idxCol;
		SClassAtom<T> *classifier = new SClassAtom<T> ();
		boost::dynamic_bitset<> powerSet;
		double a, b, ad, bd, Nt; // Nt = number of samples in current node
		double nominal;
		uli_t depVarClassSize = data.varCategories[depVar].size();
		//double stdDecOfNominal = 0.0;
		double ntwmofbest;
		double pln, pld, prn, prd, pno;

		a = b = ad = bd = 0.0;
		decOfNominal = -std::numeric_limits<double>::infinity();

		// initial count

		k = j = 0;
		nl.clear();
		nl.resize(depVarClassSize, 0.0);
		nrinit = nl;
		idxVec.clear();
		idxRow = rowMaskVec.begin();
		j = 0;
		Nt = pno = 0;
		while (idxRow != rowMaskVec.end()) {
			// get index for next class in y
			idxVec.push_back(data.depIdxVec[*idxRow]);

			// number of samples (weighted) of node
			nrinit[idxVec[j]] += rowWeight[*idxRow];
			Nt += rowWeight[*idxRow];

			++j;
			++idxRow;
		}
		ntwmofbest = Nt;

		pno = std::inner_product(nrinit.begin(), nrinit.end(), nrinit.begin(), 0);

		idxCol = colMaskVec.begin();
		while (idxCol != colMaskVec.end()) {
			//skip dep var
			if (*idxCol == depVar) {
				++idxCol;
				continue;
			}

			nominal = 0.0;

			// create class indexes of current row
			idxRowVec.clear();
			for (j = 0; j < rowMaskVec.size(); ++j) {
				valCell = data.at(rowMaskVec[j], *idxCol);
				idxCell = data.find(valCell, data.varClassIndexer[*idxCol]);
				idxRowVec.push_back(idxCell);
			}

			// get variable's classes
			Helper::getClasses<uli_t>(idxRowVec, classVarVec);

			// break if there is just one class
			if (classVarVec.size() == 1) {
				//nominal = stdDecOfNominal;
				nominal = 0;
				if (nominal > decOfNominal) {
					setCur.clear();
					decOfNominal = nominal;
					nominalBestIdx = *idxCol;
					bestSubSet = setCur;
				}
				++idxCol;
				continue;
			}

			// check all combinations
			// => check all sets in power set 2^(numberofclasses-1) - 1
			// (avoiding duplicates)

			// get powerset (binary representation)
			powerSet.reset();
			powerSet.resize(classVarVec.size(), false);

			// check all sets
			do {
				// count number of samples per class in y in sub set
				tmpclass.clear();
				tmpclass.resize(data.varClassIndexer[data.par.depVar].size(), 0);
				for (j = 0; j < idxVec.size(); ++j)
					if (powerSet[idxRowVec[j]])
						tmpclass[idxVec[j]] += (uli_t) rowWeight[rowMaskVec[j]];

				prn = pln = pld = prd = 0.0;

				for (j = 0; j < tmpclass.size(); ++j) {
					pln += tmpclass[j] * tmpclass[j];
					pld += tmpclass[j];
				}

				for (j = 0; j < tmpclass.size(); ++j) {
					tmpclass[j] = (uli_t) nrinit[j] - tmpclass[j];
					prn += tmpclass[j] * tmpclass[j];
				}

				prd = (Nt - pld);

				a = (pld == 0) ? (0) : (pln / pld);
				b = (prd == 0) ? (0) : (prn / prd);

				nominal = a + b;

				// is it a good nominal?
				if (nominal > decOfNominal) {
					Helper::getPowerSet(powerSet, classVarVec, setCur);
					decOfNominal = nominal;
					nominalBestIdx = *idxCol;
					bestSubSet = setCur;
				}

				// get next set
				Helper::inc(powerSet);

			} while (!powerSet[powerSet.size() - 1]);

			++idxCol;
		}

		decOfNominal = decOfNominal - pno / Nt;

		// get output vector
		Helper::getIndexedVector(data.varClassIndexer[nominalBestIdx], bestSubSet,
				colVec);

		// build classifier
		classifier->setVarID(nominalBestIdx);
		classifier->setClassSet(colVec);
		classifier->setMissingCode(data.getMissingCode());

		//add importance
		((TImportance<T>*) fFP.importance)->add(nominalBestIdx, decOfNominal
				/ data.par.ntree);

		return classifier;
	}

	template<class T>
	static ClassAtom<T, uli_t> *CARTcateCont(FittingFctPar<T> &fFP) {
		DataFrame<T> &data = *fFP.data;
		uli_t depVar = fFP.depVar;
		std::vector<uli_t> &rowMaskVec = *fFP.rowMaskVec;
		std::vector<uli_t> &colMaskVec = *fFP.colMaskVec;
		std::vector<double> &rowWeight = *fFP.rowWeight;
		double nodeGini = fFP.nodePerformance;
		//OUTPUT
		double &decOfGini = *fFP.treeImprovement;
		uli_t &giniBestIdx = *fFP.bestVar;

		uli_t j, k;
		std::vector<T> colVec;
		std::vector<T> classVarVec;
		std::vector<double> nl;
		std::vector<double> nr;
		std::vector<double> nrinit;
		std::vector<uli_t> idxVec;
		typename std::vector<uli_t>::iterator idxRow;
		typename std::vector<uli_t>::iterator idxCol;
		T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
		double bestCutoff = 0;
		double a, b, ad, bd, Nt; // Nt = number of samples in current node
		T cutoff;
		double krit;
		uli_t depVarClassSize = data.varCategories[depVar].size();
		double stdDecOfGini = 0.0;
		double ntwm, ntwmofbest;

		a = b = ad = bd = 0.0;
		decOfGini = -std::numeric_limits<double>::infinity();

		// initial count

		k = j = 0;
		nl.clear();
		nl.resize(depVarClassSize, 0.0);
		nrinit.clear();
		nrinit.resize(depVarClassSize, 0.0);
		idxVec.clear();
		idxRow = rowMaskVec.begin();
		j = 0;
		Nt = 0;
		while (idxRow != rowMaskVec.end()) {
			idxVec.push_back(data.depIdxVec[*idxRow]);
			nrinit[idxVec[j]] += rowWeight[*idxRow];
			++j;
			++idxRow;
		}

		// 1 + gini (without normalisation) of the whole node
		b  = std::inner_product(nrinit.begin(), nrinit.end(), nrinit.begin(), 0);
		bd = std::accumulate(nrinit.begin(), nrinit.end(), 0);
		(bd != 0.0) ? (b /= bd) : (b = 0.0);

		stdDecOfGini = b;
		Nt = bd;
		ntwmofbest = Nt;
		ntwm = Nt;

		idxCol = colMaskVec.begin();
		while (idxCol != colMaskVec.end()) {
			//skip dep var
			if (*idxCol == depVar) {
				++idxCol;
				continue;
			}

			krit = 0.0;

			// get variable's classes
			data.getRowMaskedCol(*idxCol, rowMaskVec, colVec);
			Helper::getClasses<T>(colVec, classVarVec);
			sort(classVarVec.begin(), classVarVec.end());
			//data.getRowMaskedClassesInCol(*idxCol, rowMaskVec, classVarVec);

			// break if there is just one class
			if (classVarVec.size() == 1) {
				krit = stdDecOfGini;
				if (krit > decOfGini) {
					decOfGini = krit;
					giniBestIdx = *idxCol;
					bestCutoff = (double) classVarVec[0];
				}
				++idxCol;
				continue;
			}

			// initial reset of the left node
			std::fill(nl.begin(), nl.end(), 0.0);
			nr = nrinit;

			// iterate over all cutoffs
			for (j = 0; j < classVarVec.size() - 1; ++j) {
				// intrinsic cutoff between two value which is
				// (v_i + v_{i+1})/2
				cutoff = (classVarVec[j]);

				k = 0;
				idxRow = rowMaskVec.begin();
				// get the differences
				while (idxRow != rowMaskVec.end()) {
					if (data.at(*idxRow, *idxCol) == cutoff) {
						nl[idxVec[k]] += rowWeight[*idxRow];
						nr[idxVec[k]] -= rowWeight[*idxRow];
					}
					++k;
					++idxRow;
				}

				// calculate updated gini
				// left child.
				a  = std::inner_product(nl.begin(), nl.end(), nl.begin(), 0);
				ad = std::accumulate(nl.begin(), nl.end(), 0);
				(ad != 0.0) ? (a /= ad) : (a = 0.0);

				// right child.
				b  = std::inner_product(nr.begin(), nr.end(), nr.begin(), 0);
				bd = std::accumulate(nr.begin(), nr.end(), 0);
				(bd != 0.0) ? (b /= bd) : (b = 0.0);

				krit = a + b;

				// is it a good gini?
				if ((krit > decOfGini) || ((gsl_rng_uniform_int(data.par.rng, 2) == 0)
						&& (krit == decOfGini))) {
					decOfGini = krit;
					giniBestIdx = *idxCol;
					bestCutoff = ((double) cutoff + (double) classVarVec[j + 1]) / 2.0;
					ntwmofbest = ntwm;
				}

			}
			++idxCol;
		}

		// weighting with sample size
		// what rf5 of breiman/cutler do
		decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

		// build classifier
		classifier->setVarID(giniBestIdx);
		classifier->setThreshold((T) bestCutoff);
		classifier->setMissingCode(data.getMissingCode());

		//add importance
		((TImportance<T>*) fFP.importance)->add(classifier->getVarID(), decOfGini
				/ data.par.ntree);

		return classifier;
	}


	template<class T>
	static void checkUsabilityCARTcateContSnp(DataFrame<T> &data) {
		uli_t depVar = data.par.depVar;

		if (data.par.memMode > 1) {
			throw Exception(ERRORCODE_64);
		}

		if (data.varClassIndexer.size() <= depVar) {
			throw Exception(ERRORCODE_19);
		}

		if ((data.varClassIndexer[depVar][0] == 0
				&& data.varClassIndexer[depVar][1] == 1)
				|| (data.varClassIndexer[depVar][0] == 1
						&& data.varClassIndexer[depVar][1] == 0)) {
		} else {
			throw Exception(ERRORCODE_20);
		}
	}


	template<class T>
	static ClassAtom<T, uli_t> *CARTcateContSnp(FittingFctPar<T> &fFP) {
		DataFrame<T> &data = *fFP.data;
		uli_t depVar = fFP.depVar;
		std::vector<uli_t> &rowMaskVec = *fFP.rowMaskVec;
		std::vector<uli_t> &colMaskVec = *fFP.colMaskVec;
		std::vector<double> &rowWeight = *fFP.rowWeight;
		double nodeGini = fFP.nodePerformance;
		//OUTPUT
		double &decOfGini = *fFP.treeImprovement;
		uli_t &giniBestIdx = *fFP.bestVar;

		uli_t j, k;
		std::vector<T> colVec;
		std::vector<T> classVarVec;
		std::vector<double> nl;
		std::vector<double> nr;
		std::vector<double> nrinit;
		std::vector<uli_t> idxVec;
		typename std::vector<uli_t>::iterator idxRow;
		typename std::vector<uli_t>::iterator idxCol;
		T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
		double bestCutoff = 0;
		double a, b, ad, bd, Nt; // Nt = number of samples in current node
		T cutoff;
		double krit;
		uli_t depVarClassSize = data.varCategories[depVar].size();
		double stdDecOfGini = 0.0;
		double ntwm, ntwmofbest;

		a = b = ad = bd = 0.0;
		decOfGini = -std::numeric_limits<double>::infinity();

		// initial count

		k = j = 0;
		nl.clear();
		nl.resize(depVarClassSize, 0.0);
		nrinit.clear();
		nrinit.resize(depVarClassSize, 0.0);
		idxVec.clear();
		idxRow = rowMaskVec.begin();
		j = 0;
		Nt = 0;
		while (idxRow != rowMaskVec.end()) {
			idxVec.push_back(data.depIdxVec[*idxRow]);
			nrinit[idxVec[j]] += rowWeight[*idxRow];
			++j;
			++idxRow;
		}

		// 1 + gini (without normalisation) of the whole node
		b  = std::inner_product(nrinit.begin(), nrinit.end(), nrinit.begin(), 0);
		bd = std::accumulate(nrinit.begin(), nrinit.end(), 0);
		(bd != 0.0) ? (b /= bd) : (b = 0.0);

		stdDecOfGini = b;
		Nt = bd;
		ntwmofbest = Nt;
		ntwm = Nt;

		// fill classVarVec
		classVarVec.push_back(0);
		classVarVec.push_back(1);

		idxCol = colMaskVec.begin();
		while (idxCol != colMaskVec.end()) {
			//skip dep var
			if (*idxCol == depVar) {
				++idxCol;
				continue;
			}

			krit = 0.0;

			// get variable's classes
			data.getRowMaskedCol(*idxCol, rowMaskVec, colVec);

			// initial reset of the left node
			std::fill(nl.begin(), nl.end(), 0.0);
			nr = nrinit;

			// iterate over all cutoffs
			for (j = 0; j < classVarVec.size() - 1; ++j) {
				// intrinsic cutoff between two value which is
				// (v_i + v_{i+1})/2
				cutoff = (classVarVec[j]);

				k = 0;
				idxRow = rowMaskVec.begin();
				// get the differences
				while (idxRow != rowMaskVec.end()) {
					if (data.at(*idxRow, *idxCol) == cutoff) {
						nl[idxVec[k]] += rowWeight[*idxRow];
						nr[idxVec[k]] -= rowWeight[*idxRow];
					}
					++k;
					++idxRow;
				}

				// calculate updated krit 
				// left child.
				a  = std::inner_product(nl.begin(), nl.end(), nl.begin(), 0);
				ad = std::accumulate(nl.begin(), nl.end(), 0);
				(ad != 0.0) ? (a /= ad) : (a = 0.0);

				// right child.
				b  = std::inner_product(nr.begin(), nr.end(), nr.begin(), 0);
				bd = std::accumulate(nr.begin(), nr.end(), 0);
				(bd != 0.0) ? (b /= bd) : (b = 0.0);

				krit = a + b;

				// is it a good gini?
				if ((krit > decOfGini) || ((gsl_rng_uniform_int(data.par.rng, 2) == 0)
						&& (krit == decOfGini))) {
					decOfGini = krit;
					giniBestIdx = *idxCol;
					bestCutoff = ((double) cutoff + (double) classVarVec[j + 1]) / 2.0;
					ntwmofbest = ntwm;
				}

			}
			++idxCol;
		}

		// weighting with sample size
		// what rf5 of breiman/cutler do
		decOfGini = ntwmofbest * (nodeGini - 1) + decOfGini;

		// build classifier
		classifier->setVarID(giniBestIdx);
		classifier->setThreshold((T) bestCutoff);
		classifier->setMissingCode(data.getMissingCode());

		//add importance
		((TImportance<T>*) fFP.importance)->add(classifier->getVarID(), decOfGini
				/ data.par.ntree);

		return classifier;
	}


	template<class T>
	static ClassAtom<T, uli_t> *CARTcateContSrt(FittingFctPar<T> &fFP) {
		DataFrame<T> &data = *fFP.data;
		uli_t depVar = fFP.depVar;
		std::vector<uli_t> &rowMaskVec = *fFP.rowMaskVec;
		std::vector<uli_t> &colMaskVec = *fFP.colMaskVec;
		std::vector<double> &rowWeight = *fFP.rowWeight;
		//OUTPUT
		double &decOfGini = *fFP.treeImprovement;
		uli_t &giniBestIdx = *fFP.bestVar;

		uli_t j, k;
		std::vector<T> colVec;
		std::vector<T> classVarVec;
		std::vector<T> giniVec;
		std::vector<T> outcomeVec;
		std::vector<std::pair<T, uli_t> > rowPairVec;
		typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
		std::vector<double> nl, nr, nrinit;
		std::vector<uli_t> idxVec;
		std::vector<uli_t> idxMissVec;
		std::vector<uli_t> maskVec;
		std::vector<uli_t> classMaskVec;
		std::vector<uli_t> rowMaskForDepVec;
		std::vector<uli_t> classOutcomeVec;
		typename std::vector<uli_t>::iterator idxRow;
		typename std::vector<uli_t>::iterator idxCol;
		T2ClassAtom<T> *classifier = new T2ClassAtom<T> ();
		double bestCutoff;
		double w, a, b, ad, bd, Nt, Nt2, noml, denoml, nomr, denomr; // Nt = number of samples in current node
		T cutoffT, cutoffTlast;
		double krit;
		uli_t depVarClassSize = data.varCategories[depVar].size();
		double stdDecOfGini;
		double ntwm, ntwmofbest;

		w = a = b = ad = bd = stdDecOfGini = bestCutoff = Nt = ntwmofbest = 0.0;

		k = j = 0;
		cutoffT = cutoffTlast = 0;
		decOfGini = -std::numeric_limits<double>::infinity();
		T missingcode = data.getMissingCode();
		nl.clear();
		nl.resize(depVarClassSize, 0.0);
		nrinit.clear();
		nrinit.resize(depVarClassSize, 0.0);

		// init count
		idxRow = rowMaskVec.begin();
		while (idxRow != rowMaskVec.end()) {
			nrinit[data.depIdxVec[*idxRow]] += rowWeight[*idxRow];
			++idxRow;
		}

		// 1 + gini (without normalisation) of the whole node
		b  = std::inner_product(nrinit.begin(), nrinit.end(), nrinit.begin(), 0);
		bd = std::accumulate(nrinit.begin(), nrinit.end(), 0);
		(bd != 0.0) ? (b /= bd) : (b = 0.0);

		Nt = bd;
		Nt2 = b;
		stdDecOfGini = b;

		idxCol = colMaskVec.begin();
		while (idxCol != colMaskVec.end()) {
			//skip dep var
			if (*idxCol == depVar) {
				++idxCol;
				continue;
			}

			krit = 0.0;

			// get column data information
			data.getRowMaskedColNoMissing2AndPair(*idxCol, rowMaskVec, rowPairVec,
					idxMissVec);

			// if data contains only missings then skip the current column
			if (rowPairVec.size() == 0) {
				++idxCol;
				continue;
			}

			// decount missings
			nr = nrinit;
			ntwmofbest = ntwm = Nt;

			sort(rowPairVec.begin(), rowPairVec.end());

			// initial reset of the left node
			std::fill(nl.begin(), nl.end(), 0.0);
			nr = nrinit;

			noml = 0.0;
			denoml = 0.0;
			nomr = Nt2;
			denomr = Nt;

			// iterate over all cutoffs in rowPairVec
			cutoffT = missingcode;
			itRowPairVec = rowPairVec.begin();
			while (itRowPairVec != rowPairVec.end()) {
				// change cutoff if necessary
				// and calc gini index
				if (cutoffT != itRowPairVec->first) {
					// first iteration check
					(missingcode == cutoffT) ? cutoffTlast = itRowPairVec->first
							: cutoffTlast = cutoffT;
					cutoffT = itRowPairVec->first;

					krit = (noml / denoml) + (nomr / denomr);

					// is it a good gini check
					if ((krit > decOfGini)
							|| ((gsl_rng_uniform_int(data.par.rng, 2) == 0) && (krit
									== decOfGini))) {
						decOfGini = krit;
						giniBestIdx = *idxCol;
						bestCutoff = ((double) cutoffT + (double) cutoffTlast) / 2;
						ntwmofbest = ntwm;
					}
				}

				// get the differences
				k = data.depIdxVec[itRowPairVec->second];
				w = rowWeight[itRowPairVec->second];

				// update contibution using binomial theorem
				noml += w * (2 * nl[k] + w);
				nomr += w * (-2 * nr[k] + w);

				// update denominators
				denoml += w;
				denomr -= w;

				// update class weights
				nl[k] += w;
				nr[k] -= w;

				++itRowPairVec;
			}
			++idxCol;
		}

		// weighting with sample size
		// what rf5 of breiman/cutler do
		decOfGini = decOfGini - Nt2 / Nt;

		// build classifier
		classifier->setVarID(giniBestIdx);
		classifier->setThreshold((T) bestCutoff);
		classifier->setMissingCode(data.getMissingCode());

		//add importance
		((TImportance<T> *) fFP.importance)->add(classifier->getVarID(), decOfGini
				/ data.par.ntree);

		return classifier;
	}

	template<class T>
	static ClassAtom<T, uli_t> *LOTUS(FittingFctPar<T> &fFP) {
		DataFrame<T> &data = *fFP.data;
		uli_t depVar = fFP.depVar;
		std::vector<uli_t> &rowMaskVec = *fFP.rowMaskVec;
		std::vector<uli_t> &colMaskVec = *fFP.colMaskVec;
		double nodeDev = fFP.nodePerformance;
		//OUTPUT
		double &decOfPurity = *fFP.treeImprovement;
		uli_t &bestIdx = *fFP.bestVar;

		uli_t j;
		std::vector<std::pair<T, uli_t> > rowPairVec;
		typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
		typename std::vector<std::pair<T, uli_t> >::reverse_iterator
				itRevRowPairVec;
		std::vector<T> colVec;
		std::vector<T> classVarVec;
		std::vector<T> giniVec;
		std::vector<T> outcomeVec;
		std::vector<double> nl, nr, nrinit;
		std::vector<T> leftNodeVec, rightNodeVec, leftOutcomeVec, rightOutcomeVec;
		std::vector<uli_t> idxVec;
		std::vector<uli_t> maskVec;
		std::vector<uli_t> classMaskVec;
		std::vector<uli_t> rowMaskForDepVec;
		std::vector<uli_t> classOutcomeVec;
		typename std::vector<uli_t>::iterator idxRow;
		typename std::vector<uli_t>::iterator idxCol;
		T2ClassAtom<T> *classifier;
		classifier = new T2ClassAtom<T> ();
		double bestCutoff = 0;
		double a, b, bno, bde, ad, bd;
		uli_t parentsBest = bestIdx;
		typename std::vector<T>::iterator itT;
		typename std::vector<T>::reverse_iterator itRevT;

		a = b = ad = bd = 0.0;
		//		decOfGini = -std::numeric_limits<double>::infinity();

		// 2 x 5 conting. table (see LOTUS paper by Loh)
		// (aff. and unaff.) x (counts of X-values grouped in 5 sets)
		uli_t J = 5;
		std::vector<T> quantVec;

		// quantile indexer
		uli_t i;
		typename std::vector<T>::iterator itq;
		std::vector<double> nReset(J, 0);
		std::vector<double> nAff;
		std::vector<double> nUnaff;
		std::vector<T> xj;
		double n = 0;
		double XL, bestXL;
		double xweighted;

		bestXL = -1.0;

		// the name "aff." or "unaff." must not accord to the real aff. status
		T afftype = data.at(rowMaskVec[0], depVar);

		// 0) search best variable
		idxCol = colMaskVec.begin();
		while (idxCol != colMaskVec.end()) {
			//skip dep var
			if (*idxCol == depVar) {
				++idxCol;
				continue;
			}

			// get variable's quantiles
			data.getRowMaskedColNoMissing(*idxCol, rowMaskVec, colVec, idxVec);

			// if data contains only missings then skip the current column
			if (colVec.size() == 0) {
				++idxCol;
				continue;
			}

			Helper::makePairVec(colVec, idxVec, rowPairVec);
			sort(rowPairVec.begin(), rowPairVec.end());

			// count samples
			n = 0;
			itRowPairVec = rowPairVec.begin();
			while (itRowPairVec != rowPairVec.end()) {
				//n += rowWeight[itRowPairVec->second];
				++n;
				++itRowPairVec;
			}

			// get quantiles
			quantVec.clear();
			quantVec.push_back(0);
			i = 0;
			j = 0;
			itRowPairVec = rowPairVec.begin();
			while (itRowPairVec != rowPairVec.end()) {
				//if (i > (uli_t) floor(n * (j + 1) / J)) {
				if (i > (uli_t) floor(rowPairVec.size() * (j + 1) / J)) {
					//quantVec[j] = colVec[i];
					quantVec[j] = itRowPairVec->first;
					quantVec.push_back(0);
					++j;
				}
				//i += (uli_t) rowWeight[itRowPairVec->second];
				++i;
				++itRowPairVec;
			}

			quantVec[j] = *(colVec.rbegin());

			// calc trend test score
			xj = quantVec;
			nAff = nUnaff = nReset;
			itRowPairVec = rowPairVec.begin();
			while (itRowPairVec != rowPairVec.end()) {
				for (i = 0; i < quantVec.size(); ++i) {
					if (i == 0) {
						if (data.at(itRowPairVec->second, *idxCol) <= quantVec[0]) {
							if (data.at(itRowPairVec->second, depVar) == afftype)
								//nAff[0] += rowWeight[itRowPairVec->second];
								++nAff[0];
							else
								//nUnaff[0] += rowWeight[itRowPairVec->second];
								++nUnaff[0];
							break;
						}
					} else {
						if (data.at(itRowPairVec->second, *idxCol) > quantVec[i - 1]
								&& data.at(itRowPairVec->second, *idxCol) <= quantVec[i]) {
							if (data.at(itRowPairVec->second, depVar) == afftype)
								//nAff[i] += rowWeight[itRowPairVec->second];
								++nAff[i];
							else
								//nUnaff[i] += rowWeight[itRowPairVec->second];
								++nUnaff[i];
							break;
						}
					}
				}
				++itRowPairVec;
			}

			xweighted = 0.0;
			for (j = 0; j < J; ++j) {
				xweighted += xj[j] * (nAff[j] + nUnaff[j]);
			}

			bno = 0.0;
			for (j = 0; j < J; ++j) {
				bno += (nAff[j] + nUnaff[j]) * (nUnaff[j] - Helper::getSum<double>(
						nUnaff)) / n * (xj[j] - xweighted / n);
			}

			bde = 0.0;
			for (j = 0; j < J; ++j) {
				bde += (nAff[j] + nUnaff[j]) * (xj[j] - xweighted / n) * (xj[j]
						- xweighted / n);
			}

			b = bno / bde;

			XL = 0.0;
			for (j = 0; j < J; ++j) {
				bno = nUnaff[j] / n - (Helper::getSum<double>(nUnaff) / n + b * (xj[j]
						- xweighted / n));
				XL += (nAff[j] + nUnaff[j]) * bno * bno * n * n
						/ Helper::getSum<double>(nAff) / Helper::getSum<double>(nUnaff);
			}

			if ((XL > bestXL) || ((gsl_rng_uniform_int(data.par.rng, 2) == 0) && (XL
					== bestXL))) {
				bestXL = XL;
				bestIdx = *idxCol;
				bestCutoff = xj[(J - 1) / 2];
			}
			++idxCol;

		}

		if (parentsBest == bestIdx) {
			decOfPurity = 0.0;
			return NULL;
		}

		// 1) search split point

		// get quantiles
		data.getRowMaskedColNoMissing2AndPair(bestIdx, rowMaskVec, rowPairVec);
		sort(rowPairVec.begin(), rowPairVec.end());
		Helper::getQuantilesLotusSplit(rowPairVec, quantVec);
		quantVec.resize(unique(quantVec.begin(), quantVec.end()) - quantVec.begin());

		// prepare for split search
		Helper::putPairFirstToVec(rowPairVec, leftNodeVec);
		data.getRowMaskedColNoMissing(depVar, rowMaskVec, leftOutcomeVec);

		// 1.1) foreach 0.3, 0.4, 0.5, 0.6, and 0.7 quantiles of X
		//FILE *outputFile = fopen("stuff", "w");
		//fclose(outputFile);
		double dev, devBest;
		T cutoffBest;
		dev = i = cutoffBest = 0;
		devBest = std::numeric_limits<double>::infinity();
		itRevRowPairVec = rowPairVec.rbegin();
		itRevT = quantVec.rbegin();

		while (itRevRowPairVec != rowPairVec.rend()) {

			// next split/quantile? then ->
			// LOTUS(S)                                
			// a best simple linear regression model to every node
			// no multiple regression, because of huge data
			if (itRevRowPairVec->first <= *itRevT) {

				//   1.1.2) perform lr on left data and save deviance
				//          if it doesnt converged then jump to step 2)
				dev = Helper::performLrAndGetDeviance(leftNodeVec, leftOutcomeVec);

				//   1.1.3) perform lr on right data and save deviance
				//          if it doesnt converged then jump to step 2)
				dev += Helper::performLrAndGetDeviance(rightNodeVec, rightOutcomeVec);

				//   1.1.4) add both left and right deviance
				//          if it is the smallest til now
				//          then store it as the best
				if ((dev > devBest) || ((gsl_rng_uniform_int(data.par.rng, 2) == 0)
						&& (dev == devBest))) {
					devBest = dev;
					cutoffBest = *itRevT;
				}

				++itRevT;
				if (itRevT == quantVec.rend())
					break;
			}

			// put one sample from left to right node
			rightNodeVec.push_back(leftNodeVec.back());
			leftNodeVec.pop_back();
			rightOutcomeVec.push_back(leftOutcomeVec.back());
			leftOutcomeVec.pop_back();

			++itRevRowPairVec;
		}

		// 2) if dev==inf then perform lr to current node and convert it to terminal node
		if (devBest == std::numeric_limits<double>::infinity()) {
			fFP.stopGrowing = true;
			return classifier;
		}

		// 3) return split point etc and leave the fitting function

		// build classifier
		classifier->setVarID(bestIdx);
		classifier->setThreshold(cutoffBest);
		classifier->setMissingCode(data.getMissingCode());

		//add importance
		decOfPurity = nodeDev - devBest;
		((TImportance<T> *) fFP.importance)->add(classifier->getVarID(),
				decOfPurity / data.par.ntree);

		return classifier;
	}

	/*
	 * check data, if LOTUSginiHybrid is usable
	 */
	template<class T>
	static void checkUsabilityLOTUS(DataFrame<T> &data) {
		uli_t depVar = data.par.depVar;

		if (data.par.memMode > 1) {
			throw Exception(ERRORCODE_64);
		}

		if (data.varClassIndexer.size() <= depVar) {
			throw Exception(ERRORCODE_19);
		}

		if ((data.varClassIndexer[depVar][0] == 0
				&& data.varClassIndexer[depVar][1] == 1)
				|| (data.varClassIndexer[depVar][0] == 1
						&& data.varClassIndexer[depVar][1] == 0)) {
		} else {
			throw Exception(ERRORCODE_20);
		}
	}

};

#endif /*FITTINGFCT_H_*/
