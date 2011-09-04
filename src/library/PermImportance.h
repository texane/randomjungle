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

#ifndef PERMIMPORTANCE_H_
#define PERMIMPORTANCE_H_

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "treedefs.h"
#include "Helper.h"
#include "TImportance.h"
#include "Exception.h"
#include "SaveCollector.h"

template<class T>
class PermImportance: public TImportance<T> {
public:
	PermImportance() :
		TImportance<T> () {
	}
	;

	virtual ~PermImportance() {
	}
	;

	virtual void save() {
		this->saveCol.par = this->par;

		Helper::getRanks(this->amounts, this->saveCol.orderRow, true);
		Helper::seq<size_t>(this->ids, 0, this->amounts.size() - 1, 1);

		*this->io->outImportance2 << this->saveCol;
	}

	virtual void load() {

	}

	void setBaseData(RJunglePar *par, RJungleIO *io, RJungleGen<T> *gen,
			DataFrame<T> *data) {
		this->par = par;
		this->io = io;
		this->gen = gen;
		this->data = data;
	}

	virtual void init() {
		this->saveCol.clear();

		// register vectors
		this->saveCol.push_back(&this->iterationVec, "iteration");
		this->saveCol.push_back(&this->ids, "id");
		this->saveCol.push_back(&this->data->varNames, "varname");

		this->saveCol.push_back(&this->amounts, "sorted_by_score");
		this->saveCol.push_back(&this->amounts_raw, "raw_score");
		this->saveCol.push_back(&this->amounts_breiman, "breiman_score");
		this->saveCol.push_back(&this->amounts_liaw, "liaw_score");
		this->saveCol.push_back(&this->amounts_meng, "meng_score");

		this->saveCol.push_back(&sd_breiman, "sd_breiman");
		this->saveCol.push_back(&sd_liaw, "sd_liaw");

		this->saveCol.push_back(&ntrees, "ntrees");

		// adjust this->saveCol
		this->saveCol.repeatLast = true;
		this->saveCol.showDepVar = false;

		this->saveCol.isAvailable.assign(this->par->ncol, true);

		reset();
	}

	virtual void reset() {
		// reset vectors
		ntrees.assign(this->par->ncol, 0);

		expX_values.assign(this->par->ncol, 0);
		expX2_breiman.assign(this->par->ncol, 0);
		expX2_liaw.assign(this->par->ncol, 0);

		amounts_raw.assign(this->par->ncol, 0);
		amounts_breiman.assign(this->par->ncol, 0);
		amounts_liaw.assign(this->par->ncol, 0);
		amounts_meng.assign(this->par->ncol, 0);

		this->amounts.assign(this->par->ncol, 0);

		sd_breiman.assign(this->par->ncol, 0);
		sd_liaw.assign(this->par->ncol, 0);

	}

	static Importance<T>* newImportanceObject() {
		return new PermImportance<T> ();
	}

	virtual void getBestVariables(size_t size, std::vector<uli_t> &outVec) {
		size_t i, j, idx;
		i = j = 0;

		Helper::getRanks(this->amounts, this->saveCol.orderRow, true);

		//evil conversion uli_t -> size_t
		outVec.clear();
		size_t len = this->saveCol.orderRow.size();
		i = len;
		while (i > 0) {
			idx = this->saveCol.orderRow[i - 1];

			if ((this->par->backSel == bs_NEG && this->amounts[idx] < 0) || (j
					>= size))
				this->saveCol.isAvailable[idx] = 0;

			if ((this->saveCol.isAvailable[idx] == 1) && (j < size)) {
				outVec.push_back(idx);
				++j;
			}

			--i;
		}
	}

	void finalize() {
		//double variation, accur;
		size_t i;
		uli_t ntree;

#ifdef HAVE_MPI
		ntree = (this->par->mpiId == 0)?(this->par->ntreeMpi):(this->par->ntree);
#else
		ntree = this->par->ntree;
#endif

		for (i = 0; i < expX_values.size(); ++i) {

			// get unscaled value
			amounts_raw[i] = expX_values[i] / ntree;
			amounts_meng[i] = ntrees[i] ? (expX_values[i] / ntrees[i]) : 0;

			// standard deviation (scaling)
			sd_breiman[i] = getScaling(amounts_raw[i], expX2_breiman[i]);
			sd_liaw[i] = getScaling(amounts_raw[i], expX2_liaw[i]);

			// get all scaled values
			amounts_breiman[i] = (sd_breiman[i] > 0) ? (amounts_raw[i]
					/ sd_breiman[i]) : 0;
			amounts_liaw[i] = (sd_liaw[i] > 0) ? (amounts_raw[i] / sd_liaw[i]) : 0;
		}

		// get variation
		switch (this->par->impMeasure) {
		case im_perm_raw:
			this->amounts = amounts_raw;
			break;
		case im_perm_breiman:
			this->amounts = amounts_breiman;
			break;
		case im_perm_liaw:
			this->amounts = amounts_liaw;
			break;
		case im_perm_meng:
			this->amounts = amounts_meng;
			break;
		default:
			this->amounts = amounts_raw;
			break;
		}
	}

	double getScaling(double accur, double variation) {
		// store decrease of accuracy in vector
		// this approximation results in a more stable variance estimate
		uli_t ntree;
#ifdef HAVE_MPI
		ntree = (this->par->mpiId == 0)?(this->par->ntreeMpi):(this->par->ntree);
#else
		ntree = this->par->ntree;
#endif

		variation /= ntree;
		// next, variation equals standard error
		// (a standard error calculation with div. by k - 1 would be biased)
		variation -= accur * accur;
		variation = sqrt(variation / ntree);
		return variation;
	}

	/*
	 * make permutation importance
	 *
	 */
	void add(CmpldTree<T> *cmpldTree, DataTreeSet &oobSet, long int treeId,
			Prediction<T> &predOrg, std::vector<uli_t> *&colMaskVec) {

		T * vec;
		T * bakVec;
		std::vector<uli_t> idx, rowMaskVec;
		double curacc;
		std::vector<uli_t>::iterator itIdx;
		std::vector<CmpldTree<T> *> cmpldTreeBox;
		Prediction<T> predPerm;
		double dump, accur, variation, diffAccur;
		uli_t i, j;

		// reset tree / accuracy size
		predPerm.initOne(this->par, this->io, this->gen, this->data);

		// create variable vector
		if (colMaskVec == NULL) {
			for (uli_t i = 0; i < this->par->ncol; ++i) {
				if (i == this->par->depVar)
					continue;
				idx.push_back(i);
			}
		} else {
			idx.assign(colMaskVec->begin(), colMaskVec->end());
		}

		// get row mask
		Helper::convertOOBtoMask(oobSet, treeId, rowMaskVec);

		// no samples?
		if (rowMaskVec.size() == 0)
			return;

		// permutation vectors
		vec = new T[rowMaskVec.size()];
		bakVec = new T[rowMaskVec.size()];

		// get cutoffs in tree for use in conditional importance
		std::vector<std::vector<T> > cutoffs;

		if (this->par->condimp >= 0) {
			for (j = 0; j < idx.size(); ++j) {
				cutoffs.push_back(std::vector<T>());
				this->gen->fct.getCutoffs(*cmpldTree, idx[j], cutoffs[j]);
			}
		}

		// perform permutation importance
		i = 0;
		itIdx = idx.begin();
		bool keepRunning = true;

		if (itIdx != idx.end()) {
			do {
				// if variable in tree
				if ((cmpldTree->hasVarID(*itIdx)) && (this->saveCol.isAvailable[*itIdx]
						== 1)) {
					// permutate current OOB column

					this->data->getColMaskedRow(*itIdx, rowMaskVec, vec);
					memcpy(bakVec, vec, rowMaskVec.size() * sizeof(T));

					dump = 0;
					accur = 0;
					variation = 0;

					++ntrees[*itIdx];

					// calc. tree accuracy / variation
					// add it to global accuracy / variation
					// all calculate importance
					// get decrease of accuracy

					if (this->par->condimp >= 0) { // conditional importance

						std::vector<uli_t> groups, grid;
						grid.clear();
						grid.resize(rowMaskVec.size(), 0);

						// TODO: for nominal

						if (this->par->treeType != 1 && this->par->treeType != 3
								&& this->par->treeType != 5)
							throw Exception(ERRORCODE_62);

						// continuous
						// get conditional group using pearson's correlation coefficient
						for (j = 0; j < idx.size(); ++j) {
							if ((idx[j] == this->par->depVar) || (idx[j] == *itIdx))
								continue;

							if (fabs(Helper::getPearsonCor<T>(*this->data, rowMaskVec,
									*itIdx, idx[j])) >= this->par->condimp) {

								this->gen->fct.getCutoffs(*cmpldTree, idx[j], cutoffs[j]);

								if (cutoffs[j].size() == 0)
									continue;

								Helper::getIntervalGroup(*this->data, rowMaskVec, idx[j],
										cutoffs[j], groups);

								Helper::addToCIGrid(groups, grid);

							}
						}

						Helper::permuteInGrid<T>(this->par->rng, vec, grid);
					} else { // standard importance
						gsl_ran_shuffle(this->par->rng, vec, rowMaskVec.size(), sizeof(T));
					}

					this->data->setColMaskedRow(*itIdx, rowMaskVec, vec);

					// get new accuracy
					predPerm.addOne(cmpldTree, &oobSet, treeId);

					curacc = this->gen->fct.getAccuracyOfCmpldTrees(this->data,
							predPerm.pred, this->data->getDepVar(), oobSet, treeId, true,
							false, *this->io);

					// sum of differences
					diffAccur = predOrg.acc[treeId] - curacc;
					expX_values[*itIdx] += diffAccur;
					// sum of squares

					// calculation of imp. measure in random forests (Breiman, Cutler)
					expX2_breiman[*itIdx] += pow(diffAccur, 2);

					// due to calculation of imp. measure in R's randomForest (Liaw, Wiener)
					expX2_liaw[*itIdx] += pow(diffAccur, 2) * oobSet.getSmplSize(treeId);

					// restore column
					this->data->setColMaskedRow(*itIdx, rowMaskVec, bakVec);
				}
				// prepare for next iteration
				++i;
				++itIdx;
				if (itIdx == idx.end()) {
					keepRunning = false;
				} else {
					if (*itIdx == this->par->depVar)
						++itIdx;
					if (itIdx == idx.end())
						keepRunning = false;
				}

			} while (keepRunning);

		}

		delete[] bakVec;
		delete[] vec;
	}

	/*
	 * make permutation importance
	 *
	 */
	virtual void combineMpi(std::vector<uli_t>*& colMaskVec) {
#ifdef HAVE_MPI
		MPI_Status status;
		std::vector<double> dataVec;

		MPI_Barrier(MPI_COMM_WORLD);

		*this->io->outVerbose << "Combining permutation importance scores..." << std::endl;

		if (this->par->mpiId == 0) {
			// master

			std::vector<double> message(this->par->ncol, 0);
			size_t nslaves;

			// ntree vector
			if (colMaskVec == NULL) {
				for (nslaves = 0; nslaves < (size_t)(this->par->mpiSize - 1); nslaves++) {
					// get ntrees
					MPI_Recv(&message[0], this->par->ncol, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_TAG_NTREE, MPI_COMM_WORLD, &status);
					Helper::add(ntrees, message);

					// get expX_values
					MPI_Recv(&message[0], this->par->ncol, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_TAG_EXV, MPI_COMM_WORLD, &status);
					Helper::add(expX_values, message);

					// get expX2_breiman
					MPI_Recv(&message[0], this->par->ncol, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_TAG_EX2B, MPI_COMM_WORLD, &status);
					Helper::add(expX2_breiman, message);

					// get expX2_liaw
					MPI_Recv(&message[0], this->par->ncol, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_TAG_EX2L, MPI_COMM_WORLD, &status);
					Helper::add(expX2_liaw, message);

					// print status
					//*this->io->outVerbose << nslaves + 1 << "/" << (this->par->mpiSize - 1) << std::endl;
				}
			} else {
				message.resize(colMaskVec->size());
				*this->io->outVerbose << "Variables no. " << colMaskVec->size() << std::endl;

				for (nslaves = 0; nslaves < (size_t)(this->par->mpiSize - 1); nslaves++) {
					// get ntrees
					MPI_Recv(&message[0], colMaskVec->size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_TAG_NTREE, MPI_COMM_WORLD, &status);
					Helper::addMasked(ntrees, *colMaskVec, message);

					// get expX_values
					MPI_Recv(&message[0], colMaskVec->size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_TAG_EXV, MPI_COMM_WORLD, &status);
					Helper::addMasked(expX_values, *colMaskVec, message);

					// get expX2_breiman
					MPI_Recv(&message[0], colMaskVec->size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_TAG_EX2B, MPI_COMM_WORLD, &status);
					Helper::addMasked(expX2_breiman, *colMaskVec, message);

					// get expX2_liaw
					MPI_Recv(&message[0], colMaskVec->size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_TAG_EX2L, MPI_COMM_WORLD, &status);
					Helper::addMasked(expX2_liaw, *colMaskVec, message);

					// print status
					//*this->io->outVerbose << nslaves + 1 << "/" << (this->par->mpiSize - 1) << std::endl;
				}
			}
		} else {
			// slave

			if (colMaskVec == NULL) {
				// send ntree vector
				MPI_Send(&(ntrees[0]), this->par->ncol, MPI_DOUBLE, 0, MPI_TAG_NTREE, MPI_COMM_WORLD);

				// send expX_values vector
				MPI_Send(&(expX_values[0]), this->par->ncol, MPI_DOUBLE, 0, MPI_TAG_EXV, MPI_COMM_WORLD);

				// send expX2_breiman vector
				MPI_Send(&(expX2_breiman[0]), this->par->ncol, MPI_DOUBLE, 0, MPI_TAG_EX2B, MPI_COMM_WORLD);

				// send expX2_liaw vector
				MPI_Send(&(expX2_liaw[0]), this->par->ncol, MPI_DOUBLE, 0, MPI_TAG_EX2L, MPI_COMM_WORLD);

			} else {

				// send ntree vector
				Helper::getMaskedVec(ntrees, *colMaskVec, dataVec);
				MPI_Send(&(dataVec[0]), dataVec.size(), MPI_DOUBLE, 0, MPI_TAG_NTREE, MPI_COMM_WORLD);

				// send expX_values vector
				Helper::getMaskedVec(expX_values, *colMaskVec, dataVec);
				MPI_Send(&(dataVec[0]), dataVec.size(), MPI_DOUBLE, 0, MPI_TAG_EXV, MPI_COMM_WORLD);

				// send expX2_breiman vector
				Helper::getMaskedVec(expX2_breiman, *colMaskVec, dataVec);
				MPI_Send(&(dataVec[0]), dataVec.size(), MPI_DOUBLE, 0, MPI_TAG_EX2B, MPI_COMM_WORLD);

				// send expX2_liaw vector
				Helper::getMaskedVec(expX2_liaw, *colMaskVec, dataVec);
				MPI_Send(&(dataVec[0]), dataVec.size(), MPI_DOUBLE, 0, MPI_TAG_EX2L, MPI_COMM_WORLD);

				*this->io->outVerbose << "Variables no. " << colMaskVec->size() << std::endl;
			}

			// print status
			*this->io->outVerbose << "Sent permutation scores" << std::endl;
		}

		// sync all processes
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}

	RJungleGen<T> *gen;

	std::vector<double> ntrees; // number of trees in which variable was

	std::vector<double> expX_values; // unscaled importance
	std::vector<double> expX2_breiman; // pre-scaling factor
	std::vector<double> expX2_liaw; // pre-scaling factor

	std::vector<double> amounts_raw; // pre-scaling factor
	std::vector<double> amounts_breiman;
	std::vector<double> amounts_liaw;
	std::vector<double> amounts_meng;

	std::vector<double> sd_breiman;
	std::vector<double> sd_liaw;
};

#endif /* PERMIMPORTANCE_H_ */
