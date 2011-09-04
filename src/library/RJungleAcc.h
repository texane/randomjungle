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

#ifndef RJUNGLEACC_H_
#define RJUNGLEACC_H_
/*
 * Includes
 */

#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <boost/dynamic_bitset.hpp>

#include "RJungleIO.h"
#include "Proximities.h"
#include "DataFrame.h"
#include "DataTreeSet.h"
#include "ClassAtom.h"
#include "TermClassAtom.h"
#include "Tree.h"
#include "CmpldTree.h"
#include "treedefs.h"
#include "gzstream.h"

/*
 * Source (Def. & Decl.)
 */

template<class T>
class RJungleAcc {
public:
  RJungleAcc();
  virtual ~RJungleAcc();

  static double getAccuracyCmpld(DataFrame<T> *data, T* pred, uli_t depVar,
      DataTreeSet &dataSet, long int treeId, bool one, bool showInfos,
      RJungleIO &io) {

    T *classMeVec;
    std::vector<T> classRes;
    std::vector<T> classVec;
    std::vector<T> colVec;
    uli_t i, k, jj, j, m;
    T majorityVote;
    std::vector<uli_t>::iterator it;
    typename std::vector<T>::iterator itT;
    typename std::vector<double>::iterator itDbl;
    typename std::vector<std::vector<double> > tbl;
    typename std::vector<std::vector<double> >::iterator itTbl;
    typename std::map<T, uli_t> indexer;
    typename std::vector<T> indexerInv;
    typename std::map<T, uli_t>::iterator itIdx;
    typename std::vector<double> depClassCount;

    data->getCol(depVar, colVec);
    Helper::getClasses<T>(colVec, classVec);
    sort(classVec.begin(), classVec.end());

    // build index table
    jj = 0;
    depClassCount = std::vector<double>(classVec.size(), 0.0);
    itT = classVec.begin();
    while (itT != classVec.end()) {
      if (*itT == data->getMissingCode()) {
        ++itT;
        continue;
      }

      indexer.insert(std::pair<T, uli_t>(*itT, jj));
			indexerInv.push_back(*itT);
      ++itT;
      ++jj;
    }

    // fill cont. table
    for (i = 0; i < jj; i++) {
      tbl.push_back(std::vector<double>(jj, 0.0));
    }

    // classify samples
    for (j = 0; j < data->getnrow(); ++j) {

      if (j >= dataSet.getNsmpl())
        continue;

      // get sample
      data->getRow(j, classMeVec);

      if (data->at(j, depVar) == data->getMissingCode())
        continue;

      // classify
      if (treeId != -1) {
        // one tree
        majorityVote = (dataSet.at(j, treeId)) ? (one ? pred[j] : pred[j
            + data->getnrow() * treeId]) : data->getMissingCode();

      } else {
        // forest
        classRes.clear();
        for (m = 0; m < dataSet.getNtree(); ++m) {
          if (dataSet.at(j, m))
            classRes.push_back(pred[j + data->getnrow() * m]);
				}

        majorityVote = Helper::getMostFreqProp<T>(data->par, classRes,
            data->getMissingCode(), data, showInfos, io.outVotes, j);
      }

      // store results
      if (majorityVote != data->getMissingCode()) {
        depClassCount[indexer.find(data->at(j, depVar))->second]++;

				// get index of majority vote
				itIdx = indexer.find(majorityVote);

				// what if class was not found? (important for prediction)
				if (itIdx == indexer.end()) {

					// add class to indexer
					indexer.insert(std::pair<T, uli_t>(majorityVote, jj));
					indexerInv.push_back(majorityVote);
					k = jj;
					++jj;

					// extend table
					for (i = 0; i < jj; i++) {
						if (i < tbl.size())
							tbl[i].push_back(0);
						else
							tbl.push_back(std::vector<double>(jj, 0));
					}
				} else {
					k = itIdx->second;
				}

				// get index of reponse of sample j
        i = indexer.find(data->at(j, depVar))->second;

				// add sample to table
				tbl[i][k]++;
      }

    }

    double allsamples, truesamples;
    allsamples = truesamples = 0.0;

    // get infos
    for (j = 0; j < tbl.size(); ++j) {
      truesamples += tbl[j][j];
      allsamples += depClassCount[j];
    }

    if (allsamples == 0)
      return 0;

    // print first row
    if (showInfos) {
      // header of first file
      *io.outConfusion
          << "(real outcome == rows / predicted outcome == columns )\t\n";
      *io.outConfusion << "\t";

			// print classes of y
			itIdx = indexer.begin();
			while(itIdx != indexer.end()) {
        *io.outConfusion << (double) itIdx->first << "\t";
				++itIdx;
			}
			*io.outConfusion << "error" << std::endl;

			// print remaining rows
      j = 0;
      itTbl = tbl.begin();
      while (itTbl != tbl.end()) {
        *io.outConfusion << (double) indexerInv[j] << "\t";

        itDbl = itTbl->begin();
        while (itDbl != itTbl->end()) {
          *io.outConfusion << *itDbl << "\t";
          ++itDbl;
        }

				if (depClassCount[j] == 0) {
					*io.outConfusion << "-";
				} else {
					*io.outConfusion << 1.0 - (tbl[j][j] / depClassCount[j]);
				}

				if (io.isClassTree(data->par)) {
					if (depClassCount[j] == 0) {
						*io.outConfusion2 << "-" << data->par.delimiter;
					} else {
						*io.outConfusion2 << 1.0 - (tbl[j][j] / depClassCount[j]) << data->par.delimiter;
					}
				}

        *io.outConfusion << std::endl;
        ++j;
        ++itTbl;
      }
      itTbl = tbl.begin();
      while (itTbl != tbl.end()) {
        *io.outConfusion << "\t";
        ++itTbl;
      }
      if (io.isClassTree(data->par))
        *io.outConfusion2 << 1.0 - truesamples / allsamples << std::endl;
      *io.outConfusion << "\t" << 1.0 - truesamples / allsamples << std::endl;
    }

    return truesamples / allsamples;
  }

  static double getAccuracyCmpldSos(DataFrame<T> *data, T* pred, uli_t depVar,
      DataTreeSet &dataSet, long int treeId, bool one, bool showInfos,
      RJungleIO &io) {

    std::vector<uli_t>::iterator it;
    //T *classMeVec;
    double yptr, ytr, ssTotal, avy, ssRes;
    uli_t nout, j, i, m;
    T y;

    ssTotal = ssRes = avy = i = 0;

    //data->getRow((uli_t)3, classMeVec);

    // all samples
    for (j = 0; j < data->getnrow(); ++j) {
      // sample in a tree?
      if (j >= dataSet.getNsmpl())
        throw Exception(ERRORCODE_31);

      // missing at y?
      if (data->at(j, depVar) == data->getMissingCode())
        continue;

      //data->getRow(j, classMeVec);

      y = data->at(j, depVar);

      // add your part to overall variance
      ssTotal += i * pow(y - avy, 2) / (i + 1);
      avy = (i * avy + y) / (i + 1);
      ++i;

      // init. variables
      nout = 0;
      yptr = 0;

      if (treeId != -1) {
        if (dataSet.at(j, treeId)) {
          // one tree

          // predicted result
          ytr = one ? pred[j] : pred[j + data->getnrow() * treeId];

          if (ytr != data->getMissingCode()) {
            // get mean of prediction
            yptr = ((nout * yptr) + ytr) / (nout + 1);
            ++nout;
          }
        }
      } else {
        // all trees
        for (m = 0; m < dataSet.getNtree(); ++m) {
          if (!dataSet.at(j, m))
            continue;

          // predicted result
          ytr = pred[j + data->getnrow() * m];

          if (ytr != data->getMissingCode()) {
            // get mean of prediction
            yptr = ((nout * yptr) + ytr) / (nout + 1); //Y^
            ++nout;
          }
        }
      }

      // get error estimate for samples
      if (nout != 0) {
        ssRes += pow(y - yptr, 2);
        /*
         if (showInfos)
         *io.outVerbose << y << " " << yptr << " " << nout << std::endl;
         */
      }
    }

    if (showInfos) {
      // get overall error and variance
      // variance of y^
      *io.outConfusion << "Variance of y^ (SS_Residual) = " << ssRes / i
          << " (" << ssRes << ")" << std::endl;

      // variance of y
      *io.outConfusion << "Variance of y (SS_Total) = " << ssTotal / i << " ("
          << ssTotal << ")" << std::endl;

      // explained variance of y
      *io.outConfusion << "Explained Variance (SS_Predictor) = " << (ssTotal
          - ssRes) / i << " (" << (ssTotal - ssRes) << ")" << std::endl;

      // R^2
      *io.outConfusion
          << "Accuracy estimate = R^2 = SS_Predictor / SS_Total = " << (ssTotal
          - ssRes) / ssTotal << std::endl;

      // The bagged error estimates are known to be
      // accurate estimates of the true error values.
      *io.outConfusion << "Error estimate = 1 - Accuracy estimate = " << ssRes
          / ssTotal << std::endl;

      *io.outConfusion << std::endl;
    }

    // get overall accuracy
    if (ssTotal <= 0)
      return -1;
    else
      return 1 - ssRes / ssTotal;

  }

  static double getDevianceCmpld(DataFrame<T> *data, T* pred, uli_t depVar,
      DataTreeSet &dataSet, long int treeId, bool one, bool showInfos,
      RJungleIO &io) {

    std::vector<uli_t>::iterator it;
		
		uli_t j, m;
    double nout, dev, devTotal, y, phat, phatx, ter1, ter2, mse;

		nout = dev = devTotal = y = phat = ter1 = ter2 = mse = 0;

    // all samples
    for (j = 0; j < data->getnrow(); ++j) {
      // sample in a tree?
      if (j >= dataSet.getNsmpl())
        throw Exception(ERRORCODE_31);

      // missing at y?
      if (data->at(j, depVar) == data->getMissingCode())
        continue;

      y = data->at(j, depVar);

      // init. variables
      nout = 0;
			phat = 0;

      if (treeId != -1) {
        if (dataSet.at(j, treeId)) {
          // one tree

          // predicted result
          phat = one ? pred[j] : pred[j + data->getnrow() * treeId];

          if (phat != data->getMissingCode()) {
						// mse
						mse += pow(y - phat, 2);

            // get mean of prediction
						ter1 = (y!=0)?(y*log(phat/y)):0;
						ter2 = (y!=1)?(log((1-phat)/(1-y))):0;
            dev += ter1 + ter2;
						//std::cout << phat << "(phat) " << y << "(y) " << ter1 << "(ter1) " << ter2 << "(ter2) " << dev << std::endl;
						++nout;
          }
        }
      } else {
        // all trees
				for (m = 0; m < dataSet.getNtree(); ++m) {
					if (!dataSet.at(j, m))
						continue;

					// predicted result
					phatx = pred[j + data->getnrow() * m];

					if (phat != data->getMissingCode())
						phat += phatx;
				}

				if (phat != data->getMissingCode()) {
					// mse
					mse += pow(y - phat, 2) / data->getnrow();

					// get mean of prediction
					phat /= dataSet.getNtree();
					ter1 = (y!=0)?(y*log(phat/y)):0;
					ter2 = (y!=1)?(log((1-phat)/(1-y))):0;
					dev += (ter1 + ter2) / data->getnrow();
					++nout;
				}

      }

    }
		dev *= -2;

    if (showInfos) {
      // get overall deviance
      *io.outConfusion << "Mean Deviance (Deviance / number of samples) = " << dev << std::endl;
      *io.outConfusion << "Mean MSE (MSE / number of samples) = " << mse << std::endl;
      *io.outConfusion << std::endl;
    }

    // get overall accuracy
		return -dev;
  }

};

#endif /* RJUNGLEACC_H_ */
