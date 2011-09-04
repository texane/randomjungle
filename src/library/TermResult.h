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

#ifndef TERMRESULT_H_
#define TERMRESULT_H_

/*
 * Includes
 */
#include <vector>
#include <map>

#include "DataFrame.h"


/** 
 * \brief A set of functions which yield TermClassAtoms for terminal nodes.
 */
namespace TermResult {

	/** 
	 * \brief Get most frequent value in vector. Choose value at random 
	 * if multiple values are most frequent.
	 * 
	 * @param data global data
	 * @param par random jungle paramter
	 * @param col input vector
	 * @param rowMaskVec selection of samples
	 * @param rowWeight sample weights
	 * 
	 * @return 
	 */
	template<class T>
		TermClassAtom<T, uli_t> *getTermMostFreq(
				DataFrame<T> &data, RJunglePar &par, uli_t col,
				std::vector<uli_t> &rowMaskVec, std::vector<uli_t> &colMaskVec,
				std::vector<double> &rowWeight) {

			T mostCommonVal = getMostFreq(data, par, col, rowMaskVec, rowWeight);

			return new TermClassAtom<T, uli_t> (mostCommonVal);
		}

	/** 
	 * \brief Get mean of vector
	 * 
	 * @param dataVec input vector
	 * @param missingCode missing code
	 * 
	 * @return mean value
	 */
	template<class T>
		TermClassAtom<T, uli_t> *getTermMean(
				DataFrame<T> &data, RJunglePar &par, uli_t col,
				std::vector<uli_t> &rowMaskVec,  std::vector<uli_t> &colMaskVec,
				std::vector<double> &rowWeight) {

			T mostCommonVal = getMean(data, par, col, rowMaskVec, rowWeight);

			return new TermClassAtom<T, uli_t> (mostCommonVal);
		}


	template<class T>
		TermClassAtom<T, uli_t> *getTermLotus(
				DataFrame<T> &data, RJunglePar &par, uli_t col,
				std::vector<uli_t> &rowMaskVec, std::vector<uli_t> &colMaskVec,
				std::vector<double> &rowWeight) {

			double devBest = std::numeric_limits<double>::infinity();
			std::vector<T > var, resp; // response and predictor varaible vectors
			std::vector<uli_t > idxVec; // dump vector
			typename std::vector<uli_t>::iterator idxCol;
			LotusTermClassAtom<T, uli_t> *term = new LotusTermClassAtom<T, uli_t> ();
			double dev = std::numeric_limits<double>::infinity();
			lr_options *opts = mk_lr_options(); // init lr options
			lr_train *lrt = NULL;
			dym *factors;
			dyv *outputs;
			uli_t j;
			dev = devBest = std::numeric_limits<double>::infinity();

			// init variables & copy variables
			factors = mk_dym(rowMaskVec.size(), 1);
			outputs = mk_dyv(rowMaskVec.size());
			term->betas.clear();
			term->betas.resize(2);

			// get response variable
			data.getRowMaskedColNoMissing(data.par.depVar, rowMaskVec, resp, idxVec);

			// search best LOTUS(S) simple logistic regression model 
			idxCol = colMaskVec.begin();
			while (idxCol != colMaskVec.end()) {

				// do not take dep var
				if (*idxCol == data.par.depVar)
					continue;

				// get predictor variable
				data.getRowMaskedColNoMissing(*idxCol, rowMaskVec, var, idxVec);

				// copy vectors
				for (j = 0; j < var.size(); ++j) {
					factors->tdarr[j][0] = var[j];
					outputs->farr[j] = resp[j];
				}

				// perform logistic regression
				lrt = mk_lr_train(NULL, factors, outputs, NULL, opts);

				if (lrt->lrs->converged) {

					// get deviance of lr
					dev = lr_train_deviance(lrt);

					if (devBest > dev) {
						devBest = dev;

						// save variable ID
						term->varID = *idxCol;

						// save betas
						term->betas[0] = lrt->lrs->b0;
						term->betas[1] = lrt->lrs->b->farr[0];
					}


				}

				// if deviance of all variable's LR is Inf then take
				// one at random.
				if (devBest == std::numeric_limits<double>::infinity()) {
					// save variable ID
					term->varID = *idxCol;

					// save betas
					term->betas[0] = lrt->lrs->b0;
					term->betas[1] = lrt->lrs->b->farr[0];
				}
				free_lr_train(lrt);

				++idxCol;
			}


			// free memory
			free_lr_options(opts);
			free_dyv(outputs);
			free_dym(factors);

			return term;
		}




	/** 
	 * \brief Get most frequent value in vector. Choose value at random 
	 * if multiple values are most frequent.
	 * 
	 * @param data global data
	 * @param par random jungle paramter
	 * @param col input vector
	 * @param rowMaskVec selection of samples
	 * @param rowWeight sample weights
	 * 
	 * @return 
	 */
	template<class T>
		T getMostFreq(DataFrame<T> &data, RJunglePar &par, uli_t col,
				std::vector<uli_t> &rowMaskVec, std::vector<double> &rowWeight) {

			std::vector<T> winner;
			std::map<T, uli_t> freqMap;
			std::vector<std::pair<uli_t, T> > freqPairs;
			typename std::map<T, uli_t>::iterator itMap;
			typename std::vector<std::pair<uli_t, T> >::reverse_iterator itFreq;
			size_t i;

			if (rowMaskVec.size() == 0)
				return data.par.missingcode;

			// get class counts
			for (i = 0; i < rowMaskVec.size(); ++i) {
				if (data.at(rowMaskVec[i], col) != data.par.missingcode) {
					freqMap[data.at(rowMaskVec[i], col)]
						+= (uli_t) rowWeight[rowMaskVec[i]];
				}
			}
			if (freqMap.size() == 0) return data.par.missingcode; 
			itMap = freqMap.begin();
			while (itMap != freqMap.end()) {
				freqPairs.push_back(std::pair<uli_t, T>(itMap->second, itMap->first));
				++itMap;
			}

			sort(freqPairs.begin(), freqPairs.end());

			itFreq = freqPairs.rbegin();
			while (itFreq != freqPairs.rend() && itFreq->first
					== freqPairs.rbegin()->first) {
				winner.push_back(itFreq->second);
				++itFreq;
			}

			// Liaw's randomForest in R has proposed this:
			assert(par.rng != NULL);

			return winner[gsl_rng_uniform_int(par.rng, winner.size())];
		}

	/** 
	 * \brief Get mean of vector
	 * 
	 * @param dataVec input vector
	 * @param missingCode missing code
	 * 
	 * @return mean value
	 */
	template<class T>
		static T getMean(DataFrame<T> &data, RJunglePar &par, uli_t col,
				std::vector<uli_t> &rowMaskVec, std::vector<double> &rowWeight) {
			double mean = 0;
			double siz = 0;

			if (rowMaskVec.size() == 0)
				return data.par.missingcode;

			// get class counts
			for (size_t i = 0; i < rowMaskVec.size(); ++i) {
				if (data.at(rowMaskVec[i], col) != data.par.missingcode) {
					mean += data.at(rowMaskVec[i], col);
					++siz;
				}
			}

			return (T) (mean / siz);
		}


}

#endif /*TERMRESULT_H_*/
