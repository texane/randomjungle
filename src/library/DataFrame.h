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

#ifndef DATAFRAME_H_
#define DATAFRAME_H_

/*
 * Includes
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>
#include <cmath>
#include <map>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#include "ErrorCodes.h"
#include "init.h"
#include "ioLines.h"
#include "Exception.h"
#include "gzstream.h"
#include "BitMatrix.h"
#include "RJunglePar.h"

#ifndef NULL
#define NULL 0
#endif

/*
 * DataFrame
 */
template<class T>
class DataFrame {
public:
	DataFrame(RJunglePar &par) :
		par(par), m_data(NULL), ped_data(NULL), catLimit(5), isUnsupervised(false),
				lom(std::vector<char>(par.ncol, lom_undef)), varClassIndexer() {
		this->setMissingCode(par.missingcode);
	}
	;

	DataFrame(uli_t size) :
		m_data(NULL), ped_data(NULL), catLimit(5), isUnsupervised(false), lom(
				std::vector<char>(size, lom_undef)), varClassIndexer() {
		this->par.ncol = size;
		this->par.nrow = size;
	}
	;

	DataFrame(const DataFrame &dataFrame) {
		// copy content
		varNames = dataFrame.varNames;
		lom = dataFrame.lom;
		varCategories = dataFrame.varCategories;
		varClassIndexer = dataFrame.varClassIndexer;
		depIdxVec = dataFrame.depIdxVec;
		catLimit = dataFrame.catLimit;
		isUnsupervised = dataFrame.isUnsupervised;
		par = dataFrame.par;

		initMatrix();

		uli_t j;
		for (uli_t i = 0; i < this->par.nrow; i++) {
			memcpy(m_data[i], dataFrame.m_data[i], sizeof(T) * par.ncol);

			for (j = 0; j < this->par.nrow; j++) {
				ped_data[i][j] = dataFrame.ped_data[i][j];
			}
		}
	}

	virtual ~DataFrame() {
		if (this->m_data != NULL) {
			for (uli_t i = 0; i < this->par.nrow; i++) {
				if (this->m_data[i] != NULL)
					delete[] this->m_data[i];
			}
			delete[] this->m_data;
			this->m_data = NULL;
		}

		if (this->ped_data != NULL) {
			uli_t nrowReal = this->par.nrow;
			if (isUnsupervised)
				nrowReal /= 2;
			for (uli_t i = 0; i < nrowReal; i++) {
				if (this->ped_data[i] != NULL)
					delete[] this->ped_data[i];
			}
			delete[] this->ped_data;
			this->ped_data = NULL;
		}
	}

	void setDim(uli_t nrow, uli_t ncol) {
		this->par.ncol = ncol;
		this->par.nrow = nrow;
	}

	T getArray(uli_t i, uli_t j) {
		if (i >= 0 && i < par.nrow && j >= 0 && j < par.ncol)
			return at(i, j);
		else
			return 0;
	}

	inline T at(uli_t i, uli_t j) const {
		return rjungle::magicAt(m_data[i], j);
	}
	;

	//ROW FIRST!!!, caused by this data structure
	inline const T* const operator[](uli_t i) {
		return m_data[i];
	}
	;

	void getRowMaskedColNoMissing(uli_t j, std::vector<uli_t> &maskVec,
			std::vector<T> &outVec, std::vector<uli_t> &idxVec) const {
		outVec.clear();
		idxVec.clear();

		for (uli_t i = 0; i < maskVec.size(); ++i)
			if (at(maskVec[i], j) != par.missingcode) {
				outVec.push_back(at(maskVec[i], j));
				idxVec.push_back(maskVec[i]);
			}
	}

	void getRowMaskedColNoMissing(uli_t j, std::vector<uli_t> &maskVec,
			std::vector<T> &outVec) const {
		outVec.clear();

		for (uli_t i = 0; i < maskVec.size(); ++i)
			if (at(maskVec[i], j) != par.missingcode) {
				outVec.push_back(at(maskVec[i], j));
			}
	}

	void getRowMaskedColNoMissing2(uli_t j, std::vector<uli_t> &maskVec,
			std::vector<T> &outVec, std::vector<uli_t> &idxVec,
			std::vector<uli_t> &idxMissVec) const {
		outVec.clear();
		idxVec.clear();
		idxMissVec.clear();

		for (uli_t i = 0; i < maskVec.size(); ++i)
			if (at(maskVec[i], j) != par.missingcode) {
				outVec.push_back(at(maskVec[i], j));
				idxVec.push_back(maskVec[i]);
			} else {
				idxMissVec.push_back(maskVec[i]);
			}
	}

	void getRowMaskedColNoMissing2AndPair(uli_t j, std::vector<uli_t> &maskVec,
			std::vector<std::pair<T, uli_t> >& result, std::vector<uli_t> &idxMissVec) const {
		idxMissVec.clear();
		result.clear();

		for (uli_t i = 0; i < maskVec.size(); ++i)
			if (at(maskVec[i], j) != par.missingcode) {
				result.push_back(std::make_pair(at(maskVec[i], j), maskVec[i]));
			} else {
				idxMissVec.push_back(maskVec[i]);
			}
	}

	void getRowMaskedColNoMissing2AndPair(uli_t j, std::vector<uli_t> &maskVec,
			std::vector<std::pair<T, uli_t> >& result) const {
		result.clear();

		for (uli_t i = 0; i < maskVec.size(); ++i)
			if (at(maskVec[i], j) != par.missingcode) {
				result.push_back(std::make_pair(at(maskVec[i], j), maskVec[i]));
			}
	}

	void getRowMaskedColNoMissing3(size_t j, std::vector<size_t> &maskVec,
			std::vector<double> &outVec) const {
		outVec.clear();

		for (uli_t i = 0; i < maskVec.size(); ++i)
			if (at(maskVec[i], j) != par.missingcode)
				outVec.push_back(at(maskVec[i], j));
	}

	/*
	 * Sorts input data automatically
	 * uses: Binary search and insert
	 */
	void getRowMaskedClassesInCol(uli_t j, std::vector<uli_t> &maskVec,
			std::vector<T> &outVec) const {

		std::deque<T> classes; //for fast insertion
		uli_t len, pos, jmp, nleft, nright;
		T target, val;
		bool isLeft, isThere;

		for (uli_t i = 0; i < maskVec.size(); ++i) {
			len = classes.size();
			target = at(maskVec[i], j);

			if (len == 0) {
				classes.push_back(target);
				continue;
			}

			jmp = floor(len / 2) + 1;
			nleft = jmp - 1;
			nright = len - jmp;
			pos = jmp;
			isLeft = false;

			isThere = false;

			while (true) {
				val = classes[pos - 1];

				if (target == val) {
					isThere = true;
					break;
				} else {
					if (target < val) {
						if (nleft == 0) {
							isLeft = true;
							break;
						} else {
							if (nleft == 2) {
								jmp = 1;
							} else {
								jmp = floor((nleft) / 2) + 1;
							}
						}

						nleft = nleft - jmp;
						nright = jmp - 1;
						pos = pos - jmp;
					} else {
						if (nright == 0) {
							break;
						} else {
							if (nright == 2) {
								jmp = 1;
							} else {
								jmp = floor((nright) / 2) + 1;
							}
						}

						nright = nright - jmp;
						nleft = jmp - 1;
						pos = pos + jmp;
					}
				}
			}

			//insert if not there
			if (!isThere) {
				if (isLeft) {
					classes.insert(classes.begin() + pos - 1, target);
				} else {
					classes.insert(classes.begin() + pos, target);
				}
			}
		}
		outVec.clear();
		outVec.assign(classes.begin(), classes.end());
	}

	void getRowMaskedCol(uli_t j, std::vector<uli_t> &maskVec,
			std::vector<T> &outVec) const {

		if (outVec.size() != maskVec.size())
			outVec.resize(maskVec.size());

		for (uli_t i = 0; i < maskVec.size(); ++i) {
			outVec[i] = at(maskVec[i], j);
		}
	}

	void getRow(uli_t j, std::vector<T> &outVec) const {

		if (outVec.size() != this->par.ncol)
			outVec.resize(this->par.ncol);

		for (uli_t i = 0; i < this->par.ncol; ++i) {
			outVec[i] = at(j, i);
		}
	}

	void getRow(uli_t j, T *&outVec) const {
		outVec = (this->m_data)[j];
	}

	void getRowPieceWise(uli_t j, T *&outVec) const {
		for (uli_t i = 0; i < this->par.nrow; ++i) {
			outVec[i] = at(i, j);
		}
	}

	void getCol(uli_t j, std::vector<T> &outVec) const {

		if (outVec.size() != this->par.nrow)
			outVec.resize(this->par.nrow);

		for (uli_t i = 0; i < this->par.nrow; ++i) {
			outVec[i] = at(i, j);
		}
	}

	void getColMiss(uli_t j, std::vector<T> &outVec) const {

		if (outVec.size() != this->par.nrow)
			outVec.resize(this->par.nrow);

		for (uli_t i = 0; i < this->par.nrow; ++i) {
			if (missingsNotAllowed()) {
				outVec[i] = at(i, j);
			} else if (this->missings.at(i, j)) {
				outVec[i] = par.missingcode;
			} else {
				outVec[i] = at(i, j);
			}

		}
	}

	bool anyMissingsInCol(uli_t j) const {
		if (missingsNotAllowed())
			return false;

		for (uli_t i = 0; i < this->par.nrow; ++i)
			if (missings.at(i, j))
				return true;

		return false;
	}

	void getColOrig(uli_t j, T *inVec) const {
		for (uli_t i = 0; i < this->par.nrow / 2; ++i)
			inVec[i] = at(i, j);
	}

	void setColSynth(uli_t j, T *outVec) {
		for (uli_t i = this->par.nrow / 2; i < this->par.nrow; ++i)
			set(i, j, outVec[i - this->par.nrow / 2]);
	}

	void getColNoMissings(uli_t j, std::vector<T> &outVec) const {

		outVec.clear();

		for (uli_t i = 0; i < this->par.nrow; ++i)
			if (at(i, j) != par.missingcode)
				outVec.push_back(at(i, j));
	}

	void getColNoMissingsExt(uli_t j, std::vector<T> &outVec) const {

		outVec.clear();

		for (uli_t i = 0; i < this->par.nrow; ++i) {
			if (missingsNotAllowed())
				outVec.push_back(par.missingcode);
			else if (!(this->missings).at(i, j))
				outVec.push_back(at(i, j));
			else
				outVec.push_back(par.missingcode);
		}
	}

	RJunglePar getDataFromCSV() {

		try {
			igzstream infile(par.filename);
			this->lom.clear();

			std::istringstream linestream;
			std::istringstream element;
			std::string line;
			std::string token;
			uli_t i;
			uli_t j;
			uli_t irows;
			uli_t jcols;
			uli_t rows;
			uli_t cols;
			double val;

			rows = this->par.nrow;
			cols = this->par.ncol;

			if (!infile)
				throw Exception(ERRORCODE_7);
			if (rows == 0 || cols == 0)
				throw Exception(ERRORCODE_25);

			if (par.pedfile_flag && par.transpose_flag)
				throw Exception(ERRORCODE_44);
			if (par.pedfile_flag && (par.delimiter != ' '))
				throw Exception(ERRORCODE_60);

			isUnsupervised = false;
			if (strcmp(par.depVarName, "") == 0)
				isUnsupervised = true;

			if (par.pedfile_flag) {
				if (isUnsupervised)
					this->par.depVarName = (char *) "PHENOTYPE";

				isUnsupervised = false;
			}

			if (isUnsupervised) {
				this->varNames.push_back("Unsupervised");
				par.depVarName = (char *) "Unsupervised";
				this->lom.push_back(lom_nominal);

				if (par.transpose_flag) {
					this->par.nrow = cols * 2;
					this->par.ncol = rows + 1;
				} else {
					this->par.nrow = rows * 2;
					this->par.ncol = cols + 1;
				}
			} else {
				if (par.transpose_flag) {
					this->par.nrow = cols;
					this->par.ncol = rows;
				} else {
					this->par.nrow = rows;
					this->par.ncol = cols;
					if (par.pedfile_flag)
						this->par.ncol -= 4;
				}
			}

			initMatrix();

			// prepare string helpers
			char delims[] = "??\n\r";
			delims[0] = par.delimiter;
			delims[1] = par.delimScale;

			StringIterator strItLine(delims);

			// read data
			i = irows = 0;
			while (!infile.eof()) { // get next line in file

				if (irows == rows) // stop reading. Enough lines.
					break;

				// read next line (what about Microsoft's compiler bug?)
				std::getline(infile, line);

				// skip line when empty
				if (line.empty())
					continue;

				// init the string object with first token
				strItLine.init(line);

				// analyze line
				j = jcols = 0;
				while (!strItLine.empty()) {

					if (jcols == cols) // stop reading. Enough columns.
						break;

					if (i == 0) { // first line

						// check preamble of PED file if available
						if (par.pedfile_flag && (jcols < 6)) {
							if (strItLine.lastDelim != par.delimiter)
								Exception(ERRORCODE_55);

							checkVarNames(strItLine.token, jcols);
						}

						// process the level of measurement
						if (!par.pedfile_flag || (jcols > 3)) {

							// save variable name
							this->varNames.push_back(std::string(strItLine.token));

							// level of measurement (LOM)
							if (strItLine.lastDelim == par.delimScale) { // LOM was specified

								// get next token anyway
								strItLine.next();

								if (this->par.allcont_flag) // set all as numeric
									this->lom.push_back(lom_numeric);
								else
									// set specified LOM
									pushBackLom(strItLine.token);

							} else { // LOM was not specified
								this->lom.push_back(lom_undef);
							}
						}
					} else { // data beyond first line

						if (par.pedfile_flag && (jcols < 4)) { // save preamble of PED file

							ped_data[irows][jcols] = strItLine.token;
						} else { // save data

							val = (strcmp(strItLine.token, "NA")) ? strtod(strItLine.token,
									NULL) : MISSINGCODE;

							if (par.transpose_flag) { // data is transposed
								setArray(jcols, irows + (isUnsupervised ? 1 : 0),
										static_cast<T> (val));
							} else { // data is not transposed
								setArray(irows, jcols + (isUnsupervised ? 1 : 0)
										- (par.pedfile_flag ? 4 : 0), static_cast<T> (val));
							}
						}

					}

					++jcols;
					++j;
					strItLine.next();
				}

				if (0 != i)
					++irows;
				++i;
			}

			if (isUnsupervised) {
				this->storeHalfCategories();
				this->createUnsupervisedData();
			}

			this->setDepVarName(std::string(par.depVarName));

			this->storeCategories();

			this->makeDepVecs();

			this->getMissings();

			return par;
		} catch (std::exception &e) {
			std::string *out = new std::string("DataFrame::getDataFromCSV:");
			out->append(e.what());
			throw Exception(out->c_str());
		}
	}

	void getColMaskedRow(uli_t i, std::vector<uli_t> &maskVec, T *outVec) const {

		for (uli_t j = 0; j < maskVec.size(); ++j) {
			outVec[j] = at(maskVec[j], i);
		}
	}

	void setArray(uli_t i, uli_t j, T val) {
		if (i >= 0 && i < par.nrow && j >= 0 && j < par.ncol) {
			set(i, j, val);
		}
	}

#ifdef __SPARSE_DATA__
	inline void set(uli_t i, uli_t j, T val) {
		m_data[i][j / 4] = (m_data[i][j / 4] & ~msk[j % 4]) | (val << ofs[j % 4]);
	}
#else
	inline void set(uli_t i, uli_t j, T val) {
		this->m_data[i][j] = val;
	}
#endif

	inline void rm(uli_t i, uli_t j) {
		set(i, j, 0);
	}
	;

	void setAll(T val) {
		// dynamic allocation of 2-D arrays
		uli_t i, j;
		try {
			for (i = 0; i < this->par.nrow; i++) {
				for (j = 0; j < this->par.ncol; j++) {
					set(i, j, val);
				}
			}
		} catch (std::exception &e) {

			throw;
		}
	}

	inline uli_t getnrow() {
		return par.nrow;
	}
	;
	inline uli_t getncol() {
		return par.ncol;
	}
	;

	std::ostream &print(std::ostream &os) const {
		uli_t i, j, k;

		if (par.pedfile_flag)
			os << "FID IID PAT MAT ";

		if (this->par.ncol > 0) {
			os << varNames[0];
			os << rjungle::strLom[(size_t) lom[0]];
		}

		for (j = 1; j < this->par.ncol; j++) {
			os << par.delimiter << varNames[j];
			os << rjungle::strLom[(size_t) lom[j]];
		}
		os << std::endl;

		for (i = 0; i < this->par.nrow; i++) {
			// ped file data
			if (par.pedfile_flag) {
				for (k = 0; k < 4; k++)
					os << ped_data[i][k] << par.delimiter;
			}

			for (j = 0; j < this->par.ncol; j++) {
				// data
				if (par.transpose_flag) {
					os << (double) at(j, i);
					if (!missingsNotAllowed())
						if (missings.at(j, i))
							os << "(NA)";
				} else {
					os << (double) at(i, j);
					if (!missingsNotAllowed())
						if (missings.at(i, j))
							os << "(NA)";
				}
				os << "\t";
			}
			os << std::endl;
		}

		os << "Class-Indexer of variables:" << std::endl;

		for (i = 0; i < varClassIndexer.size(); i++) {
			os << varNames[i] << ":\t";
			for (j = 0; j < varClassIndexer[i].size(); j++) {
				os << j << " -> ";
				os << (double) varClassIndexer[i][j];
				os << "\t";
			}
			os << std::endl;
		}

		return os;
	}

	void printCSV(std::ostream &os, std::vector<uli_t> *colMaskVec = NULL,
			bool printDepVar = true) const {

		uli_t i, j, k;
		unsigned int skip = (isUnsupervised ? 1 : 0);
		uli_t outnrow = (isUnsupervised ? (par.nrow / 2) : (par.nrow));
		std::vector<std::string> varNamesSel;
		bool wasNULL = false;

		if (par.ncol < 2)
			throw Exception(ERRORCODE_36);

		//get variable names of selection
		if (colMaskVec != NULL) {
			for (i = 0; i < colMaskVec->size(); ++i) {
				varNamesSel.push_back(varNames[colMaskVec->at(i)]);
			}
		} else {
			wasNULL = true;
			colMaskVec = new std::vector<uli_t>();

			for (i = skip; i < varNames.size(); ++i) {
				colMaskVec->push_back(i);
				varNamesSel.push_back(varNames[i]);
			}
		}

		// add dep. variable
		typename std::vector<uli_t>::iterator it;
		it = std::find(colMaskVec->begin(), colMaskVec->end(), par.depVar);
		if ((it == colMaskVec->end()) && (strcmp(par.depVarName, "Unsupervised")
				!= 0) && (printDepVar)) {
			colMaskVec->push_back(par.depVar);
			varNamesSel.push_back(varNames[par.depVar]);
		}

		// print variable names
		bool isFirst = true;

		// put ped variable names
		if (par.pedfile_flag) {
			os << "FID IID PAT MAT SEX PHENOTYPE";
			isFirst = false;
		}

		for (j = 0; j < varNamesSel.size(); j++) {

			if (!par.pedfile_flag || (varNamesSel[j] != "SEX" && varNamesSel[j]
					!= "PHENOTYPE")) {
				isFirst ? (isFirst = false) : (os << par.delimiter);

				os << varNamesSel[j];

				if (lom[j] != lom_undef)
					os << par.delimScale << rjungle::strLom[(size_t) lom[j]];

			}
		}
		os << std::endl;

		// reset SEX and PHENOTYPE
		if (par.pedfile_flag) {

			// PHENOTYPE
			it = std::find(colMaskVec->begin(), colMaskVec->end(), par.depVar);
			if (colMaskVec->end() != it) {
				colMaskVec->erase(it);
			}
			colMaskVec->insert(colMaskVec->begin(), par.depVar);

			// SEX
			it = std::find(colMaskVec->begin(), colMaskVec->end(), 0);
			if (colMaskVec->end() != it) {
				colMaskVec->erase(it);
			}
			colMaskVec->insert(colMaskVec->begin(), 0);

		}

		// print data
		double val;
		if (par.transpose_flag) {
		} else {
			for (i = 0; i < outnrow; i++) {
				// ped file prefix data
				if (par.pedfile_flag) {
					for (k = 0; k < 4; k++) {
						os << ped_data[i][k] << par.delimiter;
					}
				}

				// data
				os << (double) at(i, colMaskVec->at(0));
				for (j = 1; j < colMaskVec->size(); j++) {
					os << par.delimiter;

					val = (double) at(i, colMaskVec->at(j));
					if (val == (double) par.missingcode)
						os << "NA";
					else
						os << val;
				}

				os << std::endl;
			}
		}

		if (wasNULL)
			delete colMaskVec;
	}

	inline const std::vector<std::string> &getVarNames() const {
		return this->varNames;
	}

	void setVarNames(std::vector<std::string> &names) {
		this->varNames = names;
	}

	void makeDepVecs() {
		uli_t i, idx;
		depValVec = varClassIndexer[par.depVar];

		depIdxVec.clear();
		for (i = 0; i < this->par.nrow; ++i) {
			idx = find(at(i, par.depVar), depValVec);
			depIdxVec.push_back(idx);
		}
	}

	uli_t find(T val, std::vector<T> &vec) {
		return std::find(vec.begin(), vec.end(), val) - vec.begin();
	}

	void storeCategories() {
		uli_t i, col, idx;
		typename std::map<T, uli_t>::iterator it;

		varCategories.clear();
		varClassIndexer.clear();
		for (i = 0; i < this->par.ncol; ++i) {
			varCategories.push_back(std::vector<uli_t>());
			varClassIndexer.push_back(std::vector<T>());
		}

		for (col = 0; col < this->par.ncol; ++col) {
			if (col != this->par.depVar)
				continue;

			// save number of samples in category
			for (i = 0; i < this->par.nrow; ++i) {
				idx = find(at(i, col), varClassIndexer[col]);

				// if cate. is unknown then add it
				if (idx == varClassIndexer[col].size()) {
					varClassIndexer[col].push_back(at(i, col));
					varCategories[col].push_back(0);
				}

				varCategories[col][idx]++;
			}

		}
	}

	void storeHalfCategories() {
		uli_t i, col, idx;
		typename std::map<T, uli_t>::iterator it;

		varHalfCategories.clear();
		varClassIndexer.clear();
		for (i = 0; i < this->par.ncol; ++i) {
			varHalfCategories.push_back(std::vector<uli_t>());
			varClassIndexer.push_back(std::vector<T>());
		}

		for (col = 0; col < this->par.ncol; ++col) {
			// save number of samples in category
			for (i = 0; i < this->par.nrow / 2; ++i) {
				idx = find(at(i, col), varClassIndexer[col]);

				// if cate. is unknown then add it
				if (idx == varClassIndexer[col].size()) {
					varClassIndexer[col].push_back(at(i, col));
					varHalfCategories[col].push_back(0);
				}

				varHalfCategories[col][idx]++;
			}

		}
	}

	void printSummary() {
		uli_t j;

		for (j = 0; j < varNames.size(); ++j)
			std::cout << varNames[j] << " ";
		std::cout << std::endl;

		std::cout << "categories of data:" << std::endl;
		for (uli_t i = 0; i < getncol(); ++i) {
			std::cout << "var " << varNames[i] << ": ";

			//for (j = 0; j < varCategories[i].size(); ++j) {
			//  std::cout << varClassIndexer[i][j]
			//  << " " << varCategories[i][j] << ",";
			//}
			std::cout << std::endl;
		}
	}

	void setMissingCode(int val) {
		this->par.missingcode = val;
	}
	inline T getMissingCode() {
		return this->par.missingcode;
	}

	void setDepVar(uli_t idx) {
		this->par.depVar = idx;
	}
	inline uli_t getDepVar() {
		return this->par.depVar;
	}

	void setDepVarName(std::string depVarName) {
		std::vector<std::string>::iterator it;

		it = std::find(varNames.begin(), varNames.end(), depVarName);
		if (it == varNames.end())
			throw Exception(ERRORCODE_34);
		this->par.depVar = it - varNames.begin();
	}

	void setCol(uli_t col, T *vec) {

		for (uli_t i = 0; i < this->par.nrow; i++)
			set(i, col, vec[i]);
	}

	void setColMaskedRow(uli_t i, std::vector<uli_t> &maskVec, T* inVec) {

		for (uli_t j = 0; j < maskVec.size(); ++j) {
			set(maskVec[j], i, inVec[j]);
		}
	}

	void getCol(uli_t col, T *vec) {

		for (uli_t i = 0; i < this->par.nrow; i++)
			vec[i] = at(i, col);
	}

	std::vector<uli_t>* getIndexOfVarNames(std::vector<std::string> selNames) {
		std::vector<std::string>::iterator it;
		std::vector<std::string>::iterator pos;
		std::vector<uli_t>* colMaskVec = NULL;

		it = selNames.begin();
		if (it != selNames.end())
			colMaskVec = new std::vector<uli_t>();
		while (it != selNames.end()) {
			pos = std::find(varNames.begin(), varNames.end(), *it);
			if (pos == varNames.end())
				throw Exception(ERRORCODE_41);
			colMaskVec->push_back(pos - varNames.begin());
			++it;
		}

		return colMaskVec;
	}

	uli_t getIndexOfVarName(std::string selName) {
		std::vector<std::string>::iterator pos;

		pos = std::find(varNames.begin(), varNames.end(), selName);
		if (pos == varNames.end())
			throw Exception(ERRORCODE_47);

		return pos - varNames.begin();
	}

	void getMissings() {
		unsigned int i, j;

		if (missingsNotAllowed()) {
			// get missings
			for (i = 0; i < this->par.nrow; i++) {
				for (j = 0; j < this->par.ncol; j++) {
					if (at(i, j) == getMissingCode())
						throw Exception(ERRORCODE_37);
				}
			}
		} else {
			this->missings.init(par.nrow, par.ncol);

			// get missings
			for (i = 0; i < this->par.nrow; i++) {
				for (j = 0; j < this->par.ncol; j++) {
					if (at(i, j) == getMissingCode())
						missings.at(i, j, true);
				}
			}
		}

	}

	void div(T val) {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				set(i, j, at(i, j) / val);
			}
		}
	}

	void divDiag(T val) {
		unsigned int i;

		for (i = 0; i < this->par.nrow; i++) {
			set(i, i, at(i, i) / val);
		}
	}

	void divNoDiag(T val) {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				if (i == j)
					continue;
				set(i, j, at(i, j) / val);
			}
		}
	}

	void mult(T val) {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				set(i, j, at(i, j) * val);
			}
		}
	}

	void div(DataFrame<T> &dataFrame) {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				set(i, j, at(i, j) / dataFrame.at(i, j));
			}
		}
	}

	void pow2() {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				set(i, j, pow(at(i, j), 2));
			}
		}
	}

	void squareroot() {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				set(i, j, sqrt(at(i, j)));
			}
		}
	}

	void nanToMinusFive() {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				if (at(i, j) != at(i, j)) {//NaN
					at(i, j) = -5;
				}
			}
		}
	}

	void minusXPow2(DataFrame<T> &dataFrame) {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				set(i, j, at(i, j) - pow(dataFrame.at(i, j), 2));
			}
		}
	}

	void minus(DataFrame<T> &dataFrame) {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				set(i, j, at(i, j) - dataFrame.at(i, j));
			}
		}
	}

	void add(DataFrame<T> &dataFrame) {
		unsigned int i, j;

		for (i = 0; i < this->par.nrow; i++) {
			for (j = 0; j < this->par.ncol; j++) {
				set(i, j, at(i, j) + dataFrame.at(i, j));
			}
		}
	}

	void add(unsigned int i, unsigned int j, T val) {
		set(i, j, at(i, j) + val);
	}

	void createUnsupervisedData() {
		// create randomly sorted data
		unsigned int size = this->par.nrow / 2;
		T *colSynthSrc = new T[size];
		T *colSynthDst = new T[size];
		unsigned int i;

		// over all columns
		for (unsigned int j = 0; j < this->par.ncol; ++j) {

			// synthetic class labels
			if (j == 0) {
				for (i = 0; i < size; ++i)
					set(i, 0, 0);

				for (i = size; i < this->par.nrow; ++i)
					set(i, 0, 1);

				continue;
			}

			// save sampled content of original column data
			getColOrig(j, colSynthSrc);

			gsl_ran_sample(this->par.rng, colSynthDst, size, colSynthSrc, size,
					sizeof(T));

			setColSynth(j, colSynthDst);
		}

		delete[] colSynthDst;
		delete[] colSynthSrc;
	}

	void checkVarNames(char* token, unsigned int i) {
		switch (i) {
		case 0:
			if (strcmp(token, "FID") == 0)
				return;
			break;
		case 1:
			if (strcmp(token, "IID") == 0)
				return;
			break;
		case 2:
			if (strcmp(token, "PAT") == 0)
				return;
			break;
		case 3:
			if (strcmp(token, "MAT") == 0)
				return;
			break;
		case 4:
			if (strcmp(token, "SEX") == 0)
				return;
			break;
		case 5:
			if (strcmp(token, "PHENOTYPE") == 0)
				return;
			break;
		default:
			break;
		}
		throw Exception(ERRORCODE_44);
	}

	//private:
	virtual void initMatrix() {
		// dynamic allocation of 2-D arrays
		uli_t i;

		try {
			// remove
			if (this->m_data != NULL) {
				for (uli_t i = 0; i < this->par.nrow; i++) {
					if (this->m_data[i] != NULL)
						delete[] this->m_data[i];
				}
				delete[] this->m_data;
			}

			uli_t nrowReal = this->par.nrow;
			if (isUnsupervised)
				nrowReal /= 2;

			if (this->ped_data != NULL) {
				for (uli_t i = 0; i < nrowReal; i++) {
					if (this->ped_data[i] != NULL)
						delete[] this->ped_data[i];
				}
				delete[] this->ped_data;
			}

			// allocate
			this->m_data = new T*[this->par.nrow];
			for (i = 0; i < this->par.nrow; i++) {
#ifdef __SPARSE_DATA__
				this->m_data[i] = new T[(this->par.ncol/4) + ((this->par.ncol%4)?1:0)];
#else
				this->m_data[i] = new T[this->par.ncol];
#endif
			}

			this->ped_data = new std::string*[nrowReal];
			for (i = 0; i < nrowReal; i++) {
				this->ped_data[i] = new std::string[4];
			}
		} catch (std::exception &e) {

			throw;
		}

	}

	/*
	 * Assign LOM to variable
	 */
	void pushBackLom(char *token) {
		if (strcmp(token, rjungle::strLom[lom_nominal]) == 0)
			this->lom.push_back(lom_nominal);
		else if (strcmp(token, rjungle::strLom[lom_ordinal]) == 0)
			this->lom.push_back(lom_ordinal);
		else if (strcmp(token, rjungle::strLom[lom_numeric]) == 0)
			this->lom.push_back(lom_numeric);
		else
			throw Exception(ERRORCODE_17);
	}

	inline bool missingsNotAllowed() const {
		return ((par.imputeIt < 1) && (par.extractdata_flag < 1));
	}

#ifdef __SPARSE_DATA__
	static const int msk[];
	static const int ofs[];
#endif

	RJunglePar par;

	BitMatrix missings;
	T **m_data; // data will be analyzed
	std::string **ped_data; // data from ped file if --pedfile option was chosen

	uli_t catLimit;
	bool isUnsupervised;

	std::vector<std::string> varNames; // variable names
	std::vector<char> lom; // level of measurement
	std::vector<std::vector<uli_t> > varCategories, varHalfCategories; // unique values
	std::vector<std::vector<T> > varClassIndexer; // indexes
	std::vector<uli_t> depIdxVec;
	std::vector<T> depValVec;
};

#ifdef __SPARSE_DATA__
template <class T> const int DataFrame<T>::msk[4] = {192,48,12,3};
template <class T> const int DataFrame<T>::ofs[4] = {6,4,2,0};
#endif

template<class T>
inline std::ostream & operator<<(std::ostream &os, DataFrame<T> &dataFrame) {
	return dataFrame.print(os);
}

#endif /*DATAFRAME_H_*/

