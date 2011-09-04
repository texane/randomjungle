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

#ifndef HELPER_H_
#define HELPER_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/variant.hpp>

#include "config.h"
#include "RJungleIO.h"
#include "Proximities.h"
#include "DataFrame.h"
#include "DataTreeSet.h"
#include "ClassAtom.h"
#include "TermClassAtom.h"
#include "Tree.h"
#include "CmpldTree.h"
#include "gzstream.h"

#include "lr.h"
#include "treedefs.h"

/*
 * Source (Def. & Decl.)
 */

class Helper {
public:
	Helper() {
	}
	;
	virtual ~Helper() {
	}
	;

	/**
	 * \brief Converts seconds to human readable format
	 * and prints it to output stream.
	 *
	 * @param timeEst Seconds
	 * @param outVerbose Output stream.
	 */
	static void printTime(double timeEst, std::ostream &outVerbose) {
		if ((uli_t) floor(timeEst) < 60)
			outVerbose << "~" << (uli_t) floor(timeEst) << " sec.";

		if ((uli_t) floor(timeEst) >= 60 && (uli_t) floor(timeEst) < 3600)
			outVerbose << "~" << (uli_t) floor(timeEst / 60) << " min.";

		if ((uli_t) floor(timeEst) >= 3600 && (uli_t) floor(timeEst) < 86400)
			outVerbose << "~" << (uli_t) floor(timeEst / 3600) << " hours";

		if ((uli_t) floor(timeEst) >= 86400)
			outVerbose << "~" << (uli_t) floor(timeEst / 86400) << " days";
	}

	/**
	 * \brief Reads column names from a file.
	 *
	 * @param colSelection File that contains the column names.
	 * @param delimiter Delimeter (";" is default or use new line).
	 *
	 * @return Column names.
	 */
	static std::vector<std::string> getColSelection(std::string colSelection,
			char delimiter = ';') {

		try {
			std::vector<std::string> varNames;
			std::ifstream infile(colSelection.c_str());
			std::istringstream linestream;
			std::string line;
			std::string token;

			if (!infile)
				throw Exception(ERRORCODE_7);

			// read colnames
			while (!infile.eof()) {

				// fetch line
				std::getline(infile, line);
				if (line.empty())
					continue;

				// convert line content
				linestream.clear();
				linestream.str(line);

				while (!linestream.eof()) {

					// fetch items
					std::getline(linestream, token, delimiter);

					if (token.size() == 0)
						continue;

					// *nix
					if (token.at(token.size() - 1) == '\r')
						token = token.substr(0, token.size() - 1);

					// ms win., *nix
					if (token.at(token.size() - 1) == '\n')
						token = token.substr(0, token.size() - 1);

					if (token.size() == 0)
						continue;

					varNames.push_back(token);
				}
			}
			infile.close();
			return varNames;
		} catch (std::exception &e) {
			std::string *out = new std::string("Helper::getColSelection:");
			out->append(e.what());
			throw Exception(out->c_str());
		}
		return std::vector<std::string>();
	}

	/**
	 * \brief Returns the dimensions of data in input file.
	 *
	 * @param fileName File name of input data file.
	 * @param delimiter Delimeter that seperates data cells. (";" is default).
	 *
	 * @return Dimension.
	 */
	static std::pair<uli_t, uli_t> getCSVdim(char *fileName, char delimiter = ';') {
		try {
			igzstream infile(fileName);
			std::istringstream linestream;
			std::string line;
			std::string token;
			uli_t irows;
			uli_t jcols, j;

			if (!infile)
				throw Exception(ERRORCODE_7);

			std::getline(infile, line);
			linestream.clear();
			linestream.str(line);
			jcols = 0;
			while (!linestream.eof()) {
				std::getline(linestream, token, delimiter);
				++jcols;
			}

			//read
			irows = 1;
			while (!infile.eof()) {
				std::getline(infile, line);
				linestream.clear();
				linestream.str(line);
				j = 0;
				while (!linestream.eof()) {
					std::getline(linestream, token, delimiter);
					++j;
				}
				if (j == jcols)
					++irows;
			}

			infile.close();

			// irows - 1: without variable names
			return std::pair<uli_t, uli_t>(irows - 1, jcols);

		} catch (std::exception &e) {
			std::string *out = new std::string("Helper::getCSVdim:");
			out->append(e.what());
			throw Exception(out->c_str());
		}
	}

	/**
	 * \brief Gets number of non missing values in vector.
	 *
	 * @param vec Vector.
	 * @param missingCode Missing code.
	 *
	 * @return Number of non missing values in vector.
	 */
	template<class T>
	static uli_t getSize(std::vector<T> &vec, T missingCode) {
		typename std::vector<T>::iterator iter;
		uli_t i = 0;

		for (iter = vec.begin(); iter != vec.end(); ++iter)
			if (*iter != missingCode)
				i++;

		return i;
	}

	/**
	 * \brief Purity of lotus node.
	 *
	 * @param data global data
	 * @param col variable index
	 * @param rowMaskVec selection of samples
	 * @param rowWeight weights of samples
	 *
	 * @return deviance
	 */
	template<class T>
	static double purityLotus(
			//INPUT
			DataFrame<T> &data, uli_t col, std::vector<uli_t> &rowMaskVec,
			std::vector<double> &rowWeight) {

		std::vector<T> var, resp;
		std::vector<uli_t> idxVec;

		// get predictor variable
		data.getRowMaskedColNoMissing(col, rowMaskVec, var, idxVec);

		// get response variable
		data.getRowMaskedColNoMissing(data.par.depVar, rowMaskVec, resp, idxVec);

		return performLrAndGetDeviance(var, resp);
	}

	template<class T>
	static double performLrAndGetDeviance(std::vector<T> nodeVec,
			std::vector<T> outcomeVec) {

		double dev = std::numeric_limits<double>::infinity();
		lr_options *opts = mk_lr_options(); // init lr options
		lr_train *lrt;
		dym *factors;
		dyv *outputs;
		uli_t j;

		factors = mk_dym(nodeVec.size(), 1);
		outputs = mk_dyv(nodeVec.size());

		for (j = 0; j < nodeVec.size(); ++j) {
			factors->tdarr[j][0] = nodeVec[j];
			outputs->farr[j] = outcomeVec[j];
		}

		lrt = mk_lr_train(NULL, factors, outputs, NULL, opts);

		if (lrt->lrs->converged)
			dev = lr_train_deviance(lrt);

		free_lr_options(opts);
		free_lr_train(lrt);
		free_dyv(outputs);
		free_dym(factors);

		return dev;
	}

	/**
	 * \brief Returns the Gini index of a specific data column.
	 *
	 * @param data Data frame.
	 * @param col Index of column in data frame.
	 * @param rowMaskVec Selection of samples.
	 * @param rowWeight Weights of samples.
	 *
	 * @return Gini index.
	 */
	template<class T>
	static double purityGini(
			//INPUT
			DataFrame<T> &data, uli_t col, std::vector<uli_t> &rowMaskVec,
			std::vector<double> &rowWeight) {

		std::map<T, double> classVec;
		typename std::map<T, double>::iterator iter;
		std::vector<uli_t>::iterator it = rowMaskVec.begin();
		double gini = 1;
		double n = 0.0;
		double dummy;
		T val = 0;

		while (it != rowMaskVec.end()) {
			val = data.at(*it, col);

			if (val != data.getMissingCode()) {
				classVec[val] += rowWeight[*it];
				n += rowWeight[*it];

			}

			++it;
		}

		gini = 1;

		for (iter = classVec.begin(); iter != classVec.end(); ++iter) {
			dummy = (double) iter->second / n;
			gini -= dummy * dummy;
		}

		return gini;
	}

	/**
	 * \brief Returns the purity of a specific data column.
	 * If perfect clean then return 0 otherwise 0.5.
	 *
	 * @param data Data frame.
	 * @param col Index of column in data frame.
	 * @param rowMaskVec Selection of samples.
	 * @param rowWeight Weights of samples.
	 *
	 * @return Purity.
	 */
	template<class T>
	static double purityNominal(
			// There does not exists a consistent definition of
			// purity in a nominal tree (one should use gini instead)
			// Will be used in this tool because it serve the purpose here.

			//INPUT
			DataFrame<T> &data, uli_t col, std::vector<uli_t> &rowMaskVec,
			std::vector<double> &rowWeight) {

		std::vector<uli_t>::iterator it = rowMaskVec.begin();
		T val = 0;
		bool pure = true;

		if (it != rowMaskVec.end()) {
			val = data.at(*it, col);
			while (it != rowMaskVec.end()) {
				if (val != data.at(*it, col)) {
					pure = false;
					break;
				}
				++it;
			}
		}

		return pure ? 0 : 0.5;
	}

	/**
	 * \brief Returns the purity of a specific data column.
	 * Purity equals sum of squares divided by sample size.
	 *
	 * @param data Data frame.
	 * @param col Index of column in data frame.
	 * @param rowMaskVec Selection of samples.
	 * @param rowWeight Weights of samples.
	 *
	 * @return Purity.
	 */
	template<class T>
	static double puritySos(
			//INPUT
			DataFrame<T> &data, uli_t depVar, std::vector<uli_t> &rowMaskVec,
			std::vector<double> &rowWeight) {
		double ynrinit, ynr2init, n, dummy;
		std::vector<uli_t>::iterator it;
		T val;

		ynrinit = ynr2init = n = 0.0;
		it = rowMaskVec.begin();
		while (it != rowMaskVec.end()) {
			val = data.at(*it, depVar);
			if (val != data.getMissingCode()) {
				n++;
				dummy = val;
				ynrinit += dummy;
				ynr2init += dummy * dummy;
			}
			++it;
		}

		return ynr2init - ynrinit * ynrinit / n;
	}

	/**
	 * \brief Returns classification of samples sorted by
	 * prediction class.
	 *
	 * @param data Data frame.
	 * @param rowMaskVec Selection of samples.
	 * @param classifier Classifier which classifies the samples.
	 * @param classVec Classes in vectors.
	 * @param classMaskMultiVec Sample index of different prediction classes.
	 */
	template<class T>
	static void getClassifiedVectors(
			//INPUT
			DataFrame<T> &data, std::vector<uli_t> &rowMaskVec,
			ClassAtom<T, uli_t> *classifier,
			//OUTPUT
			std::vector<uli_t> &classVec,
			std::vector<std::vector<uli_t> *> &classMaskMultiVec) {

		std::vector<uli_t> maskVec;
		T *classMeVec;
		std::vector<uli_t> classMaskVec;
		std::vector<uli_t> *newClassMaskVec;
		classMaskMultiVec.clear();

		// classify
		for (uli_t i = 0; i < rowMaskVec.size(); ++i) {
			data.getRow(rowMaskVec[i], classMeVec);
			classMaskVec.push_back(classifier->classify(classMeVec));
		}
		getClasses<uli_t> (classMaskVec, classVec);

		for (uli_t j = 0; j < classVec.size(); ++j) {
			newClassMaskVec = new std::vector<uli_t>();
			classMaskMultiVec.push_back(newClassMaskVec);

			getMask<uli_t> (classVec[j], classMaskVec, maskVec);
			getIndexedVector<uli_t> (rowMaskVec, maskVec, *newClassMaskVec);
		}
	}

	template<class T>
	static double getSum(std::vector<T> &dataVec) {
		double sum = 0;
		typename std::vector<T>::iterator it = dataVec.begin();

		while (it != dataVec.end()) {
			sum += (double) *it;
			++it;
		}
		return sum;
	}

	template<class T>
	static std::vector<T> getAbs(std::vector<T> &dataVec) {
		typename std::vector<T>::iterator it = dataVec.begin();
		std::vector<T> out;

		while (it != dataVec.end()) {
			out.push_back(fabs(*it));
			++it;
		}

		return out;
	}

	template<class T>
	static std::vector<T> getDiff(std::vector<T> &dataVec, T val) {
		typename std::vector<T>::iterator it = dataVec.begin();
		std::vector<T> out;

		while (it != dataVec.end()) {
			out.push_back(*it - val);
			++it;
		}

		return out;
	}

	/*
	 * Median Absolute Deviation
	 *
	 * More robust than IQR.
	 * adjust by a factor for asymptotically normal consistency (1.4826)
	 * The default 'constant = 1.4826' (approximately 1/ Phi^(-1)(3/4) =
	 *   '1/qnorm(3/4)') ensures consistency, i.e.,
	 *
	 *                   E[mad(X_1,...,X_n)] = sigma
	 *
	 * for X_i distributed as N(mu,sigma^2) and large n.
	 *
	 */
	template<class T>
	static T getMad(std::vector<T> &dataVec, double constant = 1.4826) {
		typename std::vector<T>::iterator it = dataVec.begin();
		std::vector<T> vec;
		T res;

		res = getMedian<T> (dataVec);
		vec = getDiff<T> (dataVec, res);
		vec = getAbs<T> (vec);
		res = getMedian<T> (vec);

		return (T) (res * constant);
	}

	/*
	 * Clac. Mean
	 */
	template<class T>
	static double getMean(std::vector<T> &dataVec) {
		double mean = 0;
		typename std::vector<T>::iterator it = dataVec.begin();

		while (it != dataVec.end()) {
			mean += (double) *it;
			++it;
		}
		return (mean / (double) dataVec.size());
	}

	template<class T>
	static T getMean2(std::vector<T> &dataVec, T missingCode) {
		double mean = 0;
		T siz = dataVec.size();
		typename std::vector<T>::iterator it = dataVec.begin();

		while (it != dataVec.end()) {
			if (*it != missingCode)
				mean += (double) *it;
			else
				--siz;
			++it;
		}
		return (T) (mean / (double) siz);
	}

	template<class T>
	static T getMean4(DataFrame<T> &data, RJunglePar &par, uli_t col,
			std::vector<uli_t> &rowMaskVec, std::vector<double> &rowWeight) {
		double mean = 0;
		double siz = 0;

		if (rowMaskVec.size() == 0)
			return data.par.missingcode;

		// get class counts
		for (size_t i = 0; i < rowMaskVec.size(); ++i) {
			if (data.at(rowMaskVec[i], col) != data.par.missingcode) {
				mean += data.at(rowMaskVec[i], col);
				siz += rowWeight[rowMaskVec[i]];
			}
		}

		return (T) (mean / siz);
	}

	template<class T>
	static T getMean3(DataFrame<T> &data, RJunglePar &par, uli_t col,
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

	template<class T>
	static T getMeanProx(unsigned int row1, std::vector<T> &dataVec,
			T missingCode, Proximities<T> &proxi) {
		unsigned int row2 = 0;
		double mean = 0;
		typename std::vector<T>::iterator it = dataVec.begin();

		while (it != dataVec.end()) {
			if (*it != missingCode)
				mean += (double) *it * proxi.at(row1, row2);
			++it;
			++row2;
		}
		return (T) (mean / (double) dataVec.size());
	}

	template<class T>
	static T getMeanProxX(unsigned int row1, std::vector<T> &dataVec,
			T missingCode, Proximities<T> &proxi) {
		unsigned int row2 = 0;
		double mean = 0;
		typename std::vector<T>::iterator it = dataVec.begin();

		while (it != dataVec.end()) {
			if (*it != missingCode)
				mean += proxi.at(row1, row2);
			++it;
			++row2;
		}
		return (T) (mean / (double) dataVec.size());
	}

	/*
	 * get order
	 */
	template<class T>
	static std::vector<size_t> getRanks(std::vector<T> &in, bool ascending = true) {
		std::vector<size_t> out;

		getRanks(in, out, ascending);

		return out;
	}

	/*
	 * get order
	 */
	template<class T>
	static void getRanks(std::vector<T> &in, std::vector<size_t> &out,
			bool ascending = true) {
		std::vector<std::pair<T, size_t> > tmp;
		size_t i;

		out.clear();

		// paste indexes to input vector
		for (i = 0; i < in.size(); ++i) {
			tmp.push_back(std::make_pair(in[i], i));
		}

		// sort it
		sort(tmp.begin(), tmp.end());

		// get ranks
		if (ascending) {
			for (i = 0; i < tmp.size(); ++i) {
				out.push_back(tmp[i].second);
			}
		} else {
			for (i = 0; i < tmp.size(); ++i) {
				out.push_back(tmp[tmp.size() - i - 1].second);
			}
		}
	}

	/*
	 * get order
	 */
	template<class T>
	static void getRanks(DataFrame<T> &data, std::vector<uli_t> &idx, uli_t j,
			std::vector<size_t> &out, bool ascending = true) {
		std::vector<std::pair<T, size_t> > tmp;
		size_t i;

		out.clear();

		// paste indexes to input vector
		for (i = 0; i < idx.size(); ++i) {
			tmp.push_back(std::make_pair(data.at(idx[i], j), i));
		}

		// sort it
		sort(tmp.begin(), tmp.end());

		// get ranks
		if (ascending) {
			for (i = 0; i < tmp.size(); ++i) {
				out.push_back(tmp[i].second);
			}
		} else {
			for (i = 0; i < tmp.size(); ++i) {
				out.push_back(tmp[tmp.size() - i - 1].second);
			}
		}
	}

	/*
	 * Clac. Median
	 */
	//  template<class T>
	//  static T getMedian(std::vector<T> dataVec) {
	//    sort(dataVec.begin(), dataVec.end());
	//    if (Helper::isOdd(dataVec.size()))
	//      return dataVec[(dataVec.size() - 1) / 2];
	//    else
	//      return static_cast<T> (static_cast<double> (dataVec[dataVec.size() / 2]
	//        + dataVec[dataVec.size() / 2 - 1]) / 2.0);
	//
	//  }

	template<class T>
	static double getMedian(std::vector<T> &dataVec) {
		sort(dataVec.begin(), dataVec.end());
		if (Helper::isOdd(dataVec.size()))
			return dataVec[(dataVec.size() - 1) / 2];
		else
			return (dataVec[dataVec.size() / 2] + dataVec[dataVec.size() / 2 - 1])
					/ 2.0;

	}

	template<class T>
	static T getMedianX2(std::vector<T> &dataVec) {
		sort(dataVec.begin(), dataVec.end());
		if (Helper::isOdd(dataVec.size()))
			return dataVec[(dataVec.size() - 1) / 2];
		else
			return (T) ((double) (dataVec[dataVec.size() / 2]
					+ dataVec[dataVec.size() / 2 - 1]) / 2);

	}

	template<class T>
	static T getMedianX(std::vector<T> &dataVec) {
		sort(dataVec.begin(), dataVec.end());
		return dataVec[(dataVec.size() - 1) / 2];
	}

	template<class T>
	static T getOne(std::vector<T> &dataVec, T missingCode) {
		return 0;
	}

	static uli_t getMostFreqIndex(std::vector<uli_t> &idxVec, uli_t termCode) {
		uli_t max = *(std::max_element(idxVec.begin(), idxVec.end()));
		std::vector<uli_t> freqVec(max, 0);
		std::vector<uli_t>::iterator it = idxVec.begin();

		while (it != idxVec.end()) {
			if (*it != termCode)
				++freqVec[*it];
			++it;
		}
		return std::max_element(freqVec.begin(), freqVec.end()) - freqVec.begin();

	}

	template<class T>
	static T getMostFreq(std::vector<T> &dataVec, T missingCode) {
		std::vector<T> winner;
		std::map<T, uli_t> freqMap;
		std::vector<std::pair<uli_t, T> > freqPairs;
		typename std::vector<T>::iterator it;
		typename std::map<T, uli_t>::iterator itMap;
		typename std::vector<std::pair<uli_t, T> >::reverse_iterator itFreq;

		if (dataVec.size() == 0)
			return missingCode;

		it = dataVec.begin();
		while (it != dataVec.end()) {
			if (*it != missingCode)
				++freqMap[*it];
			++it;
		}
		if (freqMap.size() == 0)
			return missingCode;

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

		return (winner.size() > 1) ? missingCode : winner[0];
		//return winner[0];
	}

	template<class T>
	static T getMostFreqProp(DataFrame<T> &data, RJunglePar &par, uli_t col,
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
		if (freqMap.size() == 0)
			return data.par.missingcode;

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

	template<class T>
	static T getMostFreqProp(RJunglePar &par, std::vector<T> &dataVec,
			T missingCode, DataFrame<T> *data = NULL, bool showInfos = false,
			std::ostream *outStream = NULL, uli_t sample = 0) {

		std::vector<T> winner;
		std::map<T, uli_t> freqMap;
		std::vector<std::pair<uli_t, T> > freqPairs;
		typename std::vector<T>::iterator it;
		typename std::map<T, uli_t>::iterator itMap;
		typename std::vector<std::pair<uli_t, T> >::reverse_iterator itFreq;

		if (dataVec.size() == 0)
			return missingCode;

		// get frequencies
		it = dataVec.begin();
		while (it != dataVec.end()) {
			if (*it != missingCode)
				++freqMap[*it];
			++it;
		}
		if (freqMap.size() == 0)
			return missingCode;

		// rearrange result
		itMap = freqMap.begin();
		while (itMap != freqMap.end()) {
			freqPairs.push_back(std::pair<uli_t, T>(itMap->second, itMap->first));
			++itMap;
		}

		sort(freqPairs.begin(), freqPairs.end());

		// show detailed infos
		if (showInfos && data->par.votes_flag) {
			assert(data != NULL);
			assert(outStream != NULL);

			std::vector<double>
					rates(data->varCategories[data->par.depVar].size(), 0);
			double sum, margin;
			size_t idx = 0;
			sum = 0;
			margin = -std::numeric_limits<T>::infinity();

			// get sum
			itFreq = freqPairs.rbegin();
			while (itFreq != freqPairs.rend()) {
				sum += itFreq->first;
				++itFreq;
			}

			// print
			itFreq = freqPairs.rbegin();
			while (itFreq != freqPairs.rend()) {
				idx = data->find(itFreq->second,
						data->varClassIndexer[data->par.depVar]);
				rates[idx] = itFreq->first;
				++itFreq;
			}

			// get margin
			idx = data->find(data->at(sample, data->par.depVar),
					data->varClassIndexer[data->par.depVar]);

			for (size_t i = 0; i < rates.size(); ++i)
				if ((i != idx) && (margin < rates[i]))
					margin = rates[i];

			margin = rates[idx] - margin;
			*outStream << margin / sum;

			// print rates and margin
			for (idx = 0; idx < rates.size(); ++idx) {
				*outStream << data->par.delimiter << rates[idx] / sum;
			}
			*outStream << std::endl;

		}

		// get winner(s)
		itFreq = freqPairs.rbegin();
		while (itFreq != freqPairs.rend() && itFreq->first
				== freqPairs.rbegin()->first) {
			winner.push_back(itFreq->second);
			++itFreq;
		}

		// Liaw's randomForest in R has proposed this:
		return winner[gsl_rng_uniform_int(par.rng, winner.size())];
	}

	template<class T>
	static T getMostFreqX(std::vector<T> &dataVec, T missingCode) {
		std::vector<T> winner;
		std::map<T, uli_t> freqMap;
		std::vector<std::pair<uli_t, T> > freqPairs;
		typename std::vector<T>::iterator it;
		typename std::map<T, uli_t>::iterator itMap;
		typename std::vector<std::pair<uli_t, T> >::reverse_iterator itFreq;

		if (dataVec.size() == 0)
			return missingCode;

		it = dataVec.begin();
		while (it != dataVec.end()) {
			if (*it != missingCode)
				++freqMap[*it];
			++it;
		}
		if (freqMap.size() == 0)
			return missingCode;

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

		return winner[0];
	}

	template<class T>
	static double getMode(std::vector<T> &dataVec) {
		std::vector<T> winner;
		std::map<T, uli_t> freqMap;
		std::vector<std::pair<uli_t, T> > freqPairs;
		typename std::vector<T>::iterator it;
		typename std::map<T, uli_t>::iterator itMap;
		typename std::vector<std::pair<uli_t, T> >::reverse_iterator itFreq;

		if (dataVec.size() == 0)
			return 0;

		it = dataVec.begin();
		while (it != dataVec.end()) {
			++freqMap[*it];
			++it;
		}

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

		return (double) winner[0];
	}

	template<class T>
	static T getMostFreqProx(unsigned int row1, std::vector<T> &dataVec,
			T missingCode, Proximities<T> &proxi) {

		unsigned int row2;
		std::vector<T> winner;
		std::map<T, uli_t> freqMap;
		std::vector<std::pair<uli_t, T> > freqPairs;
		typename std::vector<T>::iterator it;
		typename std::map<T, uli_t>::iterator itMap;
		typename std::vector<std::pair<uli_t, T> >::reverse_iterator itFreq;

		if (dataVec.size() == 0)
			return missingCode;

		row2 = 0;
		it = dataVec.begin();
		while (it != dataVec.end()) {
			if (*it != missingCode)
				freqMap[*it] += (uli_t) proxi.at(row1, row2);
			++it;
			++row2;
		}
		if (freqMap.size() == 0)
			return missingCode;

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

		return winner[0];
	}

	template<class T>
	static std::pair<T, double> getMostFreqProx2(unsigned int row1,
			std::vector<T> &dataVec, T missingCode, Proximities<T> &proxi) {

		uli_t row2;
		double count;
		std::vector<T> winner;
		std::map<T, uli_t> freqMap;
		std::vector<std::pair<uli_t, T> > freqPairs;
		typename std::vector<T>::iterator it;
		typename std::map<T, uli_t>::iterator itMap;
		typename std::vector<std::pair<uli_t, T> >::reverse_iterator itFreq;

		if (dataVec.size() == 0)
			return std::pair<T, double>(missingCode, 0);

		row2 = 0;
		count = 0;
		it = dataVec.begin();
		while (it != dataVec.end()) {
			if (*it != missingCode) {
				freqMap[*it] += (uli_t) proxi.at(row1, row2);
				count += proxi.at(row1, row2);
			}
			++it;
			++row2;
		}
		if (freqMap.size() == 0)
			return std::pair<T, double>(missingCode, 0);

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

		return std::pair<T, double>(winner[0], freqPairs.rbegin()->first / count);
	}

	template<class T>
	static std::pair<T, double> getMostFreqLD(unsigned int row1,
			std::vector<T> &dataVec, uli_t missingCode, std::vector<double> &ldVec) {

		uli_t col;
		double count;
		std::vector<T> winner;
		std::map<T, double> freqMap;
		std::vector<std::pair<double, T> > freqPairs;
		typename std::vector<T>::iterator it;
		typename std::map<T, double>::iterator itMap;
		typename std::vector<std::pair<double, T> >::reverse_iterator itFreq;

		if (dataVec.size() == 0)
			return std::pair<uli_t, double>(missingCode, 0);

		// count votes
		col = 0;
		count = 0;
		it = dataVec.begin();
		while (it != dataVec.end()) {
			if (*it != (T) missingCode) {
				freqMap[*it] += ldVec[col];
				count += ldVec[col];
			}
			++it;
			++col;
		}
		if (freqMap.size() == 0)
			return std::pair<T, double>(missingCode, 0);

		// create the top list
		itMap = freqMap.begin();
		while (itMap != freqMap.end()) {
			freqPairs.push_back(std::pair<double, T>(itMap->second, itMap->first));
			++itMap;
		}

		sort(freqPairs.begin(), freqPairs.end());

		// get all winners (same frequence)
		itFreq = freqPairs.rbegin();
		while (itFreq != freqPairs.rend() && itFreq->first
				== freqPairs.rbegin()->first) {
			winner.push_back(itFreq->second);
			++itFreq;
		}

		return std::pair<T, double>(winner[0], freqPairs.rbegin()->first / count);
	}

	template<class T>
	static T getMedianWithSort(std::vector<T> &dataVec) {
		sort(dataVec.begin(), dataVec.end());
		if (Helper::isOdd(dataVec.size()))
			return dataVec[(dataVec.size() - 1) / 2];
		else
			return static_cast<T> (static_cast<double> (dataVec[dataVec.size() / 2]
					+ dataVec[dataVec.size() / 2 - 1]) / 2.0);

	}

	/*
	 * Quantile
	 *
	 * quantVec_i in [0, 1] real numbers
	 */
	template<class T>
	static void getCountOfInterval(std::vector<T> &dataVec,
			std::vector<double> &quantVec, std::vector<uli_t> &outVec) {

	}

	/*
	 * Clac. Max
	 */
	template<class T>
	static T getMax(std::vector<T> &dataVec) {
		T max(-std::numeric_limits<T>::infinity());
		for (uli_t i = 0; i < dataVec.size(); ++i)
			if (max < dataVec[i])
				max = dataVec[i];
		return max;
	}

	/*
	 * Clac. Min
	 */
	template<class T>
	static T getMin(std::vector<T> &dataVec) {
		T min(std::numeric_limits<T>::infinity());
		for (uli_t i = 0; i < dataVec.size(); ++i)
			if (min > dataVec[i])
				min = dataVec[i];
		return min;
	}

	/*
	 * Sort values of a map and return keys
	 */
	template<class T, class C>
	static void inverseMap(const std::map<T, C>& o, std::map<C, T>& result) {
		result.clear();
		for (typename std::map<T, C>::const_iterator begin(o.begin()); begin
				!= o.end(); ++begin)
			result.insert(std::make_pair(begin->second, begin->first));
	}

	template<class T, class C>
	static void makePairVec(std::vector<T>& colVec, std::vector<C>& idxVec,
			std::vector<std::pair<T, C> >& result) {
		result.clear();
		typename std::vector<T>::iterator itCol(colVec.begin());
		typename std::vector<C>::iterator itIdx(idxVec.begin());

		while (itCol != colVec.end()) {
			result.push_back(std::make_pair(*itCol, *itIdx));
			++itCol;
			++itIdx;
		}
	}

	template<class T, class C>
	static void makeMap(std::vector<T>& colVec, std::vector<C>& idxVec, std::map<
			T, C>& result) {
		result.clear();
		typename std::vector<T>::iterator itCol(colVec.begin());
		typename std::vector<C>::iterator itIdx(idxVec.begin());

		while (itCol != colVec.end()) {
			result.insert(std::make_pair(*itCol, *itIdx));
			++itCol;
			++itIdx;
		}
	}

	template<class T>
	static void printVec(std::vector<T> &vec) {
		for (uli_t j = 0; j < vec.size(); ++j)
			std::cout << (double) vec[j] << " ";
		std::cout << std::endl;
	}

	template<class T>
	static void printVec(T *vec, size_t len) {
		for (uli_t j = 0; j < len; ++j)
			std::cout << (double) vec[j] << " ";
		std::cout << std::endl;
	}

	template<class T>
	static void printVec(const std::vector<T> &vec) {
		for (uli_t j = 0; j < vec.size(); ++j)
			std::cout << (double) vec[j] << " ";
		std::cout << std::endl;
	}

	template<class T>
	static bool isElement(T elem, std::vector<T> &vec) {
		for (uli_t j = 0; j < vec.size(); ++j)
			if (vec[j] == elem)
				return true;
		return false;
	}

	/*
	 * Creates a "outVec" mask vector.
	 *
	 * Algo.:
	 * 1) Clear outVec.
	 * 2) Iter. over inVec and add current index to outVec if classVal occurs.
	 */
	template<class T>
	static void getMask(T classVal, std::vector<T> &inVec,
			std::vector<uli_t> &maskVec) {

		maskVec.clear();

		for (uli_t i = 0; i < inVec.size(); ++i)
			if (classVal == inVec[i])
				maskVec.push_back(i);

	}

	template<class T>
	static void getMaskedVec(std::vector<T> &inVec, std::vector<uli_t> &maskVec,
			std::vector<T> &outVec) {

		if (outVec.size() != maskVec.size())
			outVec.resize(maskVec.size());

		for (uli_t i = 0; i < maskVec.size(); ++i) {
			outVec[i] = inVec[maskVec[i]];
		}
	}

	template<class T>
	static void getIndexedVector(std::vector<T> &inVec,
			std::vector<uli_t> &indexesVec, std::vector<T> &outVec) {

		outVec.clear();

		for (uli_t i = 0; i < indexesVec.size(); ++i)
			outVec.push_back(inVec[indexesVec[i]]);

	}

	/*
	 * Get classes of vector inVec.
	 */
	template<class T>
	static void getClasses(std::vector<T> &inVec, std::vector<T> &classVec) {

		classVec.clear();
		uli_t i;
		typename std::vector<T>::iterator it;

		for (i = 0; i < inVec.size(); ++i) {
			it = std::find(classVec.begin(), classVec.end(), inVec[i]);
			if (it == classVec.end())
				classVec.push_back(inVec[i]);
		}
	}

	template<class T>
	static void getClassesAndRate(std::vector<T> &inVec, typename std::map<T,
			unsigned int> &classMap) {

		uli_t i;
		typename std::vector<T>::iterator it;

		for (i = 0; i < inVec.size(); ++i) {
			classMap[inVec[i]]++;
		}
	}

	/*
	 * prints a tree map
	 */
	template<class C, class T, class D>
	static void printTreeMap(std::map<C, T> &mapInverse, DataFrame<D> &data,
			uli_t ntree) {

		typename std::map<C, T>::iterator itInverse;
		itInverse = mapInverse.begin();
		while (itInverse != mapInverse.end()) {
			if (itInverse->second == data.getncol())
				std::cout << "leafs: " << itInverse->first << std::endl;
			else
				std::cout << "Var: " << data.getVarNames()[itInverse->second]
						<< "\t\t,   val: " << (double) itInverse->first / (double) ntree
						<< std::endl;
			++itInverse;
		}
	}

	template<class D>
	static void printTreePair(std::vector<std::pair<double, uli_t> > &freqPairs,
			DataFrame<D> &data, uli_t ntree, uli_t iteration = 0,
			std::ostream &outstream = std::cout) {

		typename std::vector<std::pair<double, uli_t> >::iterator itFreq;

		if (iteration == 0)
			outstream << "iteration" << data.par.delimiter << "id"
					<< data.par.delimiter << "varname" << data.par.delimiter << "value"
					<< std::endl;

		itFreq = freqPairs.begin();
		while (itFreq != freqPairs.end()) {

			if (itFreq->second != data.getncol() && itFreq->second != data.par.depVar)
				outstream << iteration << data.par.delimiter << itFreq->second
						<< data.par.delimiter << data.getVarNames()[itFreq->second]
						<< data.par.delimiter << (double) itFreq->first // /(double)ntree
						<< std::endl;
			//else outstream  << "leafs: " << itFreq->first << std::endl;
			++itFreq;
		}
	}

	template<class D>
	static void printTreePair2Way(
			std::vector<std::pair<double, uli_t> > &freqPairs, std::vector<std::pair<
					uli_t, uli_t> > &indexer, DataFrame<D> &data, uli_t ntree,
			uli_t iteration = 0, std::ostream &outstream = std::cout) {

		typename std::vector<std::pair<double, uli_t> >::iterator itFreq;

		if (iteration == 0)
			outstream << "iteration" << data.par.delimiter << "id1"
					<< data.par.delimiter << "id2" << data.par.delimiter << "varname1"
					<< data.par.delimiter << "varname2" << data.par.delimiter << "value"
					<< std::endl;

		itFreq = freqPairs.begin();
		while (itFreq != freqPairs.end()) {

			outstream << iteration;
			outstream << data.par.delimiter;
			outstream << indexer[itFreq->second].first;
			outstream << data.par.delimiter;
			outstream << indexer[itFreq->second].second;
			outstream << data.par.delimiter;
			outstream << data.getVarNames()[indexer[itFreq->second].first];
			outstream << data.par.delimiter;
			outstream << data.getVarNames()[indexer[itFreq->second].second];
			outstream << data.par.delimiter << (double) itFreq->first;
			outstream << std::endl;

			++itFreq;
		}
	}

	static void makeMap2SortedInvPair(std::map<uli_t, double> &freqMap,
			std::vector<std::pair<double, uli_t> > &freqPairs) {
		std::map<uli_t, double>::iterator itMap;
		freqPairs.clear();

		itMap = freqMap.begin();
		while (itMap != freqMap.end()) {
			freqPairs.push_back(std::pair<double, uli_t>(itMap->second, itMap->first));
			++itMap;
		}

		sort(freqPairs.begin(), freqPairs.end());
	}

	static inline bool isOdd(uli_t val) {
		return 1 & val;
	}

	static inline bool isEven(uli_t val) {
		return !(1 & val);
	}

	template<class T>
	static double getPredScore(std::vector<T> &outcomeVec,
			std::vector<T> &classMeVec) {
		uli_t i;
		double truePos = 0;
		double posRate = 0;
		double trueNeg = 0;
		double negRate = 0;
		/*
		 for(i = 0; i < outcomeVec.size(); ++i) {
		 if (classMeVec[i] == 0)
		 ++posRate;
		 else
		 ++negRate;

		 if (classMeVec[i] == 0 && classMeVec[i] == classRes)
		 ++trueNeg;
		 else if (classMeVec[i] == 1 && classMeVec[i] == classRes)
		 ++truePos;
		 }
		 */
		return truePos / posRate * trueNeg / negRate;
	}

	// get propability that a variable has the chance to be chosen at least 1 time
	// in a worst case scenario (every tree has got only one node)
	static double getProp(uli_t ntree, uli_t mtry, uli_t numOfVars) {
		double p;

		p = 1 - (double) mtry / (double) numOfVars;
		p = 1 - pow(p, (double) ntree);

		return p;
	}

	// get a reduced tree size (dep. from mtry/numOfVars) that produce
	// the same propability given in double getProp(...)
	static uli_t getTreeSizeX(uli_t ntree, //old
			uli_t mtry1, //old
			uli_t numOfVars1, //old
			uli_t mtry2, //new
			uli_t numOfVars2 //new
	) {

		double nominator = log(1 - (double) mtry1 / (double) numOfVars1);
		double denominator = log(1 - (double) mtry2 / (double) numOfVars2);

#ifndef HAVE_ROUND
		return (uli_t)round((double)ntree * nominator / denominator);
#else
		return (uli_t) floor((double) ntree * nominator / denominator);
#endif
	}

	// get tree size (dep. from mtry/numOfVars) with p=...
	static uli_t getTreeSizeProb(double p, uli_t mtry, uli_t numOfVars) {

		double nominator = log(1 - p);
		double denominator = log(1 - (double) mtry / (double) numOfVars);

#ifndef HAVE_ROUND
		return (uli_t)round(nominator / denominator);
#else
		return (uli_t) floor(nominator / denominator);
#endif
	}

	static void getCSVdim(RJunglePar &par) {
		try {
			igzstream infile(par.filename);
			std::istringstream linestream;
			std::string line;
			std::string token;
			uli_t irows;
			uli_t jcols, j;

			if (!infile)
				throw Exception(ERRORCODE_7);

			std::getline(infile, line);
			linestream.clear();
			linestream.str(line);
			jcols = 0;
			while (!linestream.eof()) {
				std::getline(linestream, token, par.delimiter);
				++jcols;
			}

			//read
			irows = 1;
			while (!infile.eof()) {
				std::getline(infile, line);
				linestream.clear();
				linestream.str(line);
				j = 0;
				while (!linestream.eof()) {
					std::getline(linestream, token, par.delimiter);
					++j;
				}

				if (j == jcols) {
					++irows;
				} else if (j == 1) {
					break;
				} else {
					throw Exception(ERRORCODE_57);
				}
			}

			infile.close();

			// irows - 1: without variable names
			// get dimensions from input file
			if (par.nrow == 0)
				par.nrow = irows - 1 - par.skipRow;
			if (par.ncol == 0)
				par.ncol = jcols - par.skipCol;

		} catch (std::exception &e) {
			std::string *out = new std::string("Helper::getCSVdim:");
			out->append(e.what());
			throw Exception(out->c_str());
		}
	}

	static void add(std::vector<double> &to, std::vector<double> &from) {
		for (unsigned int i = 0; i < to.size(); ++i)
			to[i] += from[i];
	}

	static void addMasked(std::vector<double> &to, std::vector<uli_t> &maskVec,
			std::vector<double> &from) {
		for (unsigned int i = 0; i < maskVec.size(); ++i)
			to[maskVec[i]] += from[i];
	}

	static void vecToPairs(std::vector<std::pair<double, uli_t> > &pairs,
			std::vector<double> &vec) {
		pairs.clear();
		for (unsigned int vrs = 0; vrs < vec.size(); ++vrs)
			pairs.push_back(std::make_pair(vec[vrs], vrs));
	}

	static void makeToFreq(std::vector<double> &vec) {
		double total = 0;
		unsigned int i;

		for (i = 0; i < vec.size(); ++i)
			total += vec[i];
		for (i = 0; i < vec.size(); ++i)
			vec[i] /= total;
	}

	template<class T>
	static void getPowerSet(boost::dynamic_bitset<> &powerSet,
			std::vector<T> &from, std::vector<T> &to) {
		size_t i = powerSet.find_first();
		to.clear();

		// check all bits
		if (i != powerSet.npos)
			do {
				to.push_back(from[i]);
				i = powerSet.find_next(i);
			} while (i != powerSet.npos);
	}

	static void inc(boost::dynamic_bitset<> &bitset) {
		uli_t i = 0;

		do {
			bitset.flip(i);
			i++;
		} while (!bitset[i - 1]);
	}

	template<class T>
	static void getSNPsInHighLD(DataFrame<T> &data,
			std::vector<double> &ldVecFull, std::vector<double> &ldVec, std::vector<
					uli_t> *&colMaskFullVec, std::vector<uli_t> *&colMaskVec) {

		// init
		uli_t i;
		std::vector<uli_t> idxVec;
		std::vector<std::pair<double, uli_t> > ldPair;
		typename std::vector<std::pair<double, uli_t> >::iterator itPair;

		ldVec.clear();
		ldVec.resize(ldVecFull.size(), 0.0);

		// make index vector
		for (i = 0; i < ldVecFull.size(); i++) {
			idxVec.push_back(i);
		}

		// combine index with content
		Helper::makePairVec<double, uli_t>(ldVecFull, idxVec, ldPair);

		// sort
		sort(ldPair.begin(), ldPair.end());

		// copy ld vector
		colMaskVec->clear();
		i = 0;
		itPair = ldPair.begin();
		while (itPair != ldPair.end()) {
			if ((itPair->first > data.par.cutoffHighLD) || (i < 5)) {
				ldVec[itPair->second] = itPair->first;
				colMaskVec->push_back(colMaskFullVec->at(itPair->second));
			}
			++itPair;
			++i;
		}

		std::cout << "size cMFV " << colMaskFullVec->size() << std::endl;
		std::cout << "size lVF " << ldVecFull.size() << std::endl;
		std::cout << "size cMV " << colMaskVec->size() << std::endl;
		std::cout << "size lV " << ldVec.size() << std::endl;

	}

	template<class T>
	static void getSNPsInLD(DataFrame<T> &data, unsigned int col, std::vector<
			uli_t> *colMaskVec, std::vector<double> &ldVec) {
	}

	template<class T>
	static void vectorToXML(std::vector<T> &vec, std::string &name,
			std::string &xmlString) {

		typename std::vector<T>::iterator it;
		xmlString.clear();

		xmlString.append("<vector name=\"");
		xmlString.append(name);
		xmlString.append("\">");

		it = vec.begin();
		while (it != vec.end()) {
			xmlString.append(std::string(*it));
			++it;
			if (it != vec.end())
				xmlString.append(",");
		}

		xmlString.append("</vector>");
	}

	static void strSplit(std::string &str, std::string delim, std::vector<
			std::string> &results) {
		int cutAt;
		results.clear();

		while ((cutAt = str.find_first_of(delim)) != (signed) str.npos) {
			if (cutAt > 0)
				results.push_back(str.substr(0, cutAt));
			str = str.substr(cutAt + 1);
		}
		if (str.length() > 0)
			results.push_back(str);
	}

	template<class T>
	static void copy(std::vector<T> vec1, T* vec2) {
		for (size_t i = 0; i < vec1.size(); ++i)
			vec2[i] = vec1[i];
	}

	template<class T>
	static double one(DataFrame<T> &data, uli_t col,
			std::vector<uli_t> &rowMaskVec, std::vector<double> &rowWeight) {
		return 1;
	}

	static void strSplit(char *in, const char *delim,
			std::vector<std::string> &out) {

		char *result = NULL;

		out.clear();
		result = strtok(in, delim);
		while (result != NULL) {
			out.push_back(result);
			result = strtok(NULL, delim);
		}
	}

	template<class T>
	static bool isIn(T val, std::vector<T> &vec) {
		for (size_t i = 0; i < vec.size(); ++i) {
			if (val == vec[i])
				return true;
		}

		return false;
	}

	/*
	 * Create a vector for iterative search:
	 * Theta is a tuning parameter. The higher the more values.
	 */
	static void getSlicedVec(uli_t ncol, std::vector<uli_t> &mtrys, double theta) {

		// declare
		double diff, step, halfTime;
		size_t i;

		// define
		mtrys.clear();
		i = 0;
		step = ncol * theta;
		halfTime = ncol - step;
		diff = step;
		if (floor(diff) < 2)
			diff = 2;

		while (diff <= ncol) {
			mtrys.push_back((uli_t) floor(diff));
			diff += step;

			if (mtrys[i] == floor(diff)) {
				diff = mtrys[i] + 1;
			}

			++i;
		}
	}

	template<class T>
	static void clearPtrVec(std::vector<T> &vec) {
		for (size_t i = 0; i < vec.size(); ++i) {
			if (vec[i] != NULL) {
				delete vec[i];
				vec[i] = NULL;
			}
		}

		vec.clear();
	}

	template<class T>
	static void seq(std::vector<T> &vec, T from, T to, T step) {
		vec.clear();
		vec.push_back(from);
		while (vec[vec.size() - 1] + step <= to) {
			vec.push_back(vec[vec.size() - 1] + step);
		}
	}

	static bool isClassTree(RJunglePar &par) {
		return (par.treeType == tt_CART) || (par.treeType == tt_CARTcateCate)
				|| (par.treeType == tt_CARTsrt);
	}

	static bool isGuessingLomNumeric(RJunglePar &par) {
		return (par.treeType == tt_CART) || (par.treeType == tt_CARTcontCont)
				|| (par.treeType == tt_CARTsrt);
	}

	/*!
	 * \brief Calculates the sample correlation coefficient, i.e.,
	 * Pearson's correlation coefficient between two variables
	 *
	 * @param idx Index vector of samples
	 * @param j1 Index of 1st variable
	 * @param j2 Index of 2nd variable
	 * @return Sample correlation coefficient
	 */
	template<class T>
	static double getPearsonCor(DataFrame<T> &data, std::vector<uli_t> &idx,
			uli_t j1, uli_t j2) {
		double cov, res, meanj1, meanj2, sd1, sd2, xmx1, xmx2;
		size_t i;

		// defining variables
		cov = res = meanj1 = meanj2 = sd1 = sd2 = xmx1 = xmx2 = 0.0;

		// get means
		for (i = 0; i < idx.size(); ++i) {
			meanj1 += (double) data.at(idx[i], j1);
			meanj2 += (double) data.at(idx[i], j2);
		}
		meanj1 /= (double) idx.size();
		meanj2 /= (double) idx.size();

		// estimate s c coefficient
		for (i = 0; i < idx.size(); ++i) {
			xmx1 = (double) data.at(idx[i], j1) - meanj1;
			xmx2 = (double) data.at(idx[i], j2) - meanj2;

			// estimate cov
			cov += xmx1 * xmx2;

			// estimate s d
			sd1 += pow(xmx1, 2);
			sd2 += pow(xmx2, 2);
		}
		res = cov / sqrt(sd1) / sqrt(sd2);

		return res;
	}

	/*!
	 * /brief Group input vector elements using interval vector. (Assign element to a interval)
	 * @param data Input data set (DataFrame)
	 * @param idx Indexes of input observations
	 * @param j1 Index of variable
	 * @param cutoffs Vector of cutoff values (intervals)
	 * @param groups Resulting interval assignment
	 */
	template<class T>
	static void getIntervalGroup(DataFrame<T> &data, std::vector<uli_t> &idx,
			uli_t j1, std::vector<T> &cutoffs, std::vector<uli_t> &groups) {

		size_t i, j;
		std::vector<size_t> order;

		// define variables
		j = 0;
		assert(cutoffs.size()> 0);

		getRanks(data, idx, j1, order);
		groups.resize(order.size(), 0);

		// assign all vec elements (observations) to a interval group
		for (i = 0; i < order.size(); ++i) {
			if (cutoffs.size() > j) {
				if (data.at(idx[order[i]], j1) > cutoffs[j])
					++j;
			}

			groups[order[i]] = j;
		}

	}

	/*!
	 * /brief Convert oobSet sample information to a traditional rowMaskVec
	 * @param oobSet OOB set
	 * @param treeId ID of current tree
	 * @param rowMaskVec Output vector
	 */
	static void convertOOBtoMask(DataTreeSet &oobSet, long int treeId,
			std::vector<uli_t> &rowMaskVec) {
		uli_t i;

		rowMaskVec.clear();

		for (i = 0; i < oobSet.getNsmpl(); ++i) {
			if (oobSet.at(i, treeId)) {
				rowMaskVec.push_back(i);
			}
		}
	}

	/*!
	 * /brief Expand conditional importance grid by a new variable
	 * @param groups
	 * @param grid
	 */
	static void addToCIGrid(std::vector<uli_t> &groups, std::vector<uli_t>&grid) {

		if (grid.size() == 0) {
			grid = groups;
			return;
		}

		assert(grid.size() == groups.size());

		uli_t i, j;
		std::map<std::pair<uli_t, uli_t>, std::vector<uli_t> > gridMap;
		std::map<std::pair<uli_t, uli_t>, std::vector<uli_t> >::iterator iter;

		// store concordant groups
		for (i = 0; i < grid.size(); ++i) {
			gridMap[std::pair<uli_t, uli_t>(groups[i], grid[i])].push_back(i);
		}

		// get grid
		j = 0;
		iter = gridMap.begin();
		while (iter != gridMap.end()) {

			for (i = 0; i < iter->second.size(); ++i) {
				grid[iter->second[i]] = j;
			}

			++j;
			++iter;
		}
	}

	/*!
	 * /brief Permutes values in vec using conditional grid
	 * @param rng Random number generator (GSL)
	 * @param vec Vector of Observations
	 * @param grid Conditional grid
	 */
	template<class T>
	static void permuteInGrid(gsl_rng *rng, T *vec, std::vector<uli_t> &grid) {

		uli_t i, j;
		std::map<uli_t, std::vector<uli_t> > gridMap;
		std::map<uli_t, std::vector<uli_t> >::iterator iter;
		T *tmpVecPtr;

		// store grid groups in a map
		for (i = 0; i < grid.size(); ++i) {
			gridMap[grid[i]].push_back(i);
		}

		// permute each grid group
		j = 0;
		iter = gridMap.begin();
		while (iter != gridMap.end()) {
			tmpVecPtr = new T[iter->second.size()];

			// collect next group in group
			for (i = 0; i < iter->second.size(); ++i) {
				tmpVecPtr[i] = vec[iter->second[i]];
			}

			// permute group
			if (iter->second.size() > 1)
				gsl_ran_shuffle(rng, tmpVecPtr, iter->second.size(), sizeof(T));

			// save permuted elements
			for (i = 0; i < iter->second.size(); ++i) {
				vec[iter->second[i]] = tmpVecPtr[i];
			}

			delete[] tmpVecPtr;

			++j;
			++iter;
		}
	}

	template<class T>
	static void getQuantiles(std::vector<double> &qvec, std::vector<std::pair<T,
			uli_t> > &rowPairVec, std::vector<T> &quantVec) {

		typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
		uli_t i, j;

		i = j = 0;
		quantVec.clear();

		itRowPairVec = rowPairVec.begin();
		while (itRowPairVec != rowPairVec.end()) {
			if (i > (uli_t) floor(rowPairVec.size() * qvec[j])) {
				quantVec.push_back(itRowPairVec->first);
				if (j == (qvec.size() - 1))
					break;
				++j;
			}
			//i += (uli_t) rowWeight[itRowPairVec->second];
			++i;
			++itRowPairVec;
		}
	}

	template<class T>
	static void getQuantilesLotusSplit(
			std::vector<std::pair<T, uli_t> > &rowPairVec, std::vector<T> &quantVec) {
		double dd[] = { .3, .4, .5, .6, .7 };
		std::vector<double> qvec(dd, dd + 5);

		Helper::getQuantiles(qvec, rowPairVec, quantVec);
	}

	template<class T>
	static void putPairFirstToVec(std::vector<std::pair<T, uli_t> > &rowPairVec,
			std::vector<T> &vec) {

		typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
		vec.clear();

		itRowPairVec = rowPairVec.begin();
		while (itRowPairVec != rowPairVec.end()) {
			vec.push_back(itRowPairVec->first);
			++itRowPairVec;
		}
	}

	template<class T>
	static void putPairSecondToVec(std::vector<std::pair<T, uli_t> > &rowPairVec,
			std::vector<uli_t> &vec) {

		typename std::vector<std::pair<T, uli_t> >::iterator itRowPairVec;
		vec.clear();

		itRowPairVec = rowPairVec.begin();
		while (itRowPairVec != rowPairVec.end()) {
			vec.push_back(itRowPairVec->second);
			++itRowPairVec;
		}
	}

};

#endif /*HELPER_H_*/
