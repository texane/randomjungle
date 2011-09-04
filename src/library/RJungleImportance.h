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

#ifndef RJUNGLEIMPORTANCE_H_
#define RJUNGLEIMPORTANCE_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "RJungleGen.h"
#include "RJungleIO.h"
#include "RJunglePar.h"
#include "CmpldTree.h"
#include "Helper.h"
#include "treedefs.h"

template<class T>
class RJungleImportance {
public:
  RJungleImportance();
  virtual ~RJungleImportance();

  /*
   * make permutation importance
   *
   */
  static void makePermVarImp2(
      RJunglePar &par, RJungleIO &io, RJungleGen<T> &gen, DataFrame<T> &data,
      DataTreeSet &oobSet, std::vector<uli_t> *&colMaskVec, std::vector<
          CmpldTree<T> *> &cmpldTrees,
      std::vector<std::pair<double, uli_t> > &permVarImpPairVec) {

    T *vec;
    T *bakVec;
    std::vector<uli_t> idx;
    std::vector<std::pair<uli_t, uli_t> > numOfTrees;
    std::vector<double> accuracy;
    //std::vector<double > variation;
    std::vector<uli_t>::iterator itIdx;
    typename std::vector<CmpldTree<T> *>::iterator itCurTree;
    uli_t curTree;
    double dump, accur, variation, k;

    clock_t start = 0;
    clock_t end = 0;
    int timeSpan = 50; // maximal 50
    int timeShift = 10;
    double timeEst = 0;
    int totalTrees = 0;
    double percentRatio = 0.1;
    uli_t percentRate = 0;
    uli_t percent = 0;
    uli_t treeBufferSize;
    uli_t maxTreeBufferSize = 10;
    int i;

    // reset tree / accuracy size
    cmpldTrees.resize(par.ntree, NULL);
    accuracy.resize(par.ntree, -10);
    //variation.resize(ntree, 0);

    permVarImpPairVec.clear();

    *io.outVerbose << "Calculating permutation importance..." << std::endl;

    // create variable vector
    if (colMaskVec == NULL) {
      for (uli_t i = 0; i < par.ncol; ++i) {
        if (i == par.depVar)
          continue;
        idx.push_back(i);
      }
    } else {
      idx.assign(colMaskVec->begin(), colMaskVec->end());
    }

    // set time span
    timeSpan = 1 + std::min((uli_t) 100, (uli_t) idx.size() / 10);
    percentRate = std::max((uli_t) 1, (uli_t) floor(idx.size() * percentRatio));

    // permutation vectors
    vec = new T[par.nrow];
    bakVec = new T[par.nrow];

    // perform perm. imp.
    i = 0;
    itIdx = idx.begin();
    bool keepRunning = true;

    // create multiple threads
#pragma omp parallel num_threads(1) default(shared)
    if (itIdx != idx.end()) {
      //DataFrame data(data);
      int ntree_omp;
      int numOfThreads_omp;
      std::vector<double> accuracy_omp;
      std::vector<double> variation_omp;
      std::vector<CmpldTree<T> *> cmpldTrees_omp;
      std::vector<Tree<T, uli_t> *> trees_omp;
      int startTree_omp = 0;
      int endTree_omp = 0;
      int id_omp;
      uli_t j_omp;
      double dumpAccur_omp, dumpVariation_omp, k_omp, curacc_omp, diffAccur_omp;

      numOfThreads_omp = 1;
      id_omp = 0;

#pragma omp master
      numOfTrees.resize(numOfThreads_omp, std::make_pair(0, 0));

#pragma omp barrier

      // balance trees over all threads
      ntree_omp = (int) floor(par.ntree / numOfThreads_omp);

      //#pragma omp critical
      //if (id_omp != 0) {totalTrees = 0; ntree_omp = 0;}

#pragma omp critical
      if (id_omp != 0) {
        startTree_omp = totalTrees;
        endTree_omp = totalTrees + ntree_omp;
        numOfTrees[id_omp] = std::make_pair(startTree_omp, endTree_omp);
        totalTrees += ntree_omp;
      }

#pragma omp barrier

#pragma omp master
      {
        ntree_omp = (int) par.ntree - totalTrees;
        startTree_omp = totalTrees;
        endTree_omp = totalTrees + ntree_omp;
        numOfTrees[id_omp] = std::make_pair(startTree_omp, endTree_omp);
        treeBufferSize = 1 + std::min((uli_t) floor(par.ntree
            / numOfThreads_omp), maxTreeBufferSize);

        *io.outVerbose << numOfThreads_omp << " thread(s) permuting "
            << par.ntree << " trees" << std::endl;
        /*
         << ntree << " trees (" << ntree_omp;

         for (uli_t j = 1; j < numOfTrees.size(); ++j) {
         *io.outVerbose << "+" << numOfTrees[j].second - numOfTrees[j].first;
         }

         *io.outVerbose << ")" << std::endl;
         */
        *io.outVerbose << "Get initial accuracy" << std::flush;
      }

#pragma omp barrier

#pragma omp critical
      {
        // get original accuracy of training set
        accuracy_omp.resize(par.ntree, -10);
        variation_omp.resize(par.ntree, 0);
        cmpldTrees_omp = cmpldTrees;
      }

      for (int j = startTree_omp; j < endTree_omp; ++j) {
/*
        accuracy_omp[j] = gen.fct.getAccuracyOfCmpldTrees(
            &data, cmpldTrees_omp, gen.fct, data.getDepVar(), oobSet, j, false,
            io);
*/
      }

      // merging compiled trees
#pragma omp critical
      {
        for (int j = startTree_omp; j < endTree_omp; ++j) {
          accuracy[j] = accuracy_omp[j];
        }
        *io.outVerbose << "." << std::flush;
      }

#pragma omp barrier

      // each thread gets all info
#pragma omp critical
      {
        for (int j = 0; j < (int) accuracy.size(); ++j) {
          accuracy_omp[j] = accuracy[j];
        }
      }

#pragma omp barrier

#pragma omp master
      *io.outVerbose << std::endl << "Permuting..." << std::endl;

      // start with all sync.
      do {
        // master permutes next column
#pragma omp master
        {
          //*io.outVerbose << *itIdx << ":" << std::flush;
          // management for time estimation
          if (i < timeSpan + timeShift && i >= timeShift)
            start = clock();
          if ((uli_t) i == percent + percentRate) {
            percent = percent + percentRate;
            *io.outVerbose << "progress: " << (uli_t) (100.0 * (double) percent
                / idx.size()) << "%" << std::endl;
          }

          // permutate current column
          data.getCol(*itIdx, vec);
          memcpy(bakVec, vec, par.nrow * sizeof(T));

          curTree = 0;
          itCurTree = cmpldTrees.begin();
          dump = 0;
          accur = 0;
          variation = 0;
          k = 0;

        }

        // all are waiting for all
#pragma omp barrier

        // reset variables
        dumpAccur_omp = 0;
        dumpVariation_omp = 0;
        k_omp = 0;

        while (itCurTree != cmpldTrees.end()) {
          // read some trees in buffer
#pragma omp critical
          {
            j_omp = 0;
            startTree_omp = curTree;
            while (itCurTree != cmpldTrees.end() && j_omp < treeBufferSize) {
              itCurTree++;
              j_omp++;
              curTree++;
            }
            endTree_omp = curTree;
          }

          // calc. tree accuracy / variation
          // add it to global accuracy / variation
          // all calculate importance
          // get decrease of accuracy
          for (int j = startTree_omp; j < endTree_omp; ++j) {
            // Does tree contains the variable?
            if (!(cmpldTrees_omp[j]->hasVarID(*itIdx)))
              continue;

            gsl_ran_shuffle(par.rng, vec, par.nrow, sizeof(T));
            data.setCol(*itIdx, vec);

            // get new accuracy
/*
            curacc_omp = gen.fct.getAccuracyOfCmpldTrees(
                &data, cmpldTrees_omp, gen.fct, data.getDepVar(), oobSet, j, false,
                io);
*/

            ++k_omp;

            // sum of differences
            diffAccur_omp = accuracy_omp[j] - curacc_omp;
            dumpAccur_omp += diffAccur_omp;
            // sum of squares

            if (par.impMeasure == im_perm_breiman) {
              // calculation of imp. measure in random forests (Breiman, Cutler)
              dumpVariation_omp += pow(diffAccur_omp, 2);
            } else if (par.impMeasure == im_perm_liaw) {
              // due to calculation of imp. measure in R's randomForest (Liaw, Wiener)
              dumpVariation_omp += pow(diffAccur_omp, 2)
                  * oobSet.getSmplSize(j);
            }
          }
          //#pragma omp barrier
        }

#pragma omp critical
        if (k_omp > 0) {
          accur += dumpAccur_omp;
          variation += dumpVariation_omp;
          k += k_omp;
          //*io.outVerbose << "." << std::flush;
        }

        // end up with all sync.
#pragma omp barrier

        // master restores content of column and prepares for next iteration
#pragma omp master
        {
          if ((par.impMeasure == im_perm_breiman) || (par.impMeasure == im_perm_liaw)
              || (par.impMeasure == im_perm_raw)) {
            k = par.ntree;
          }

          if (k > 0)
            accur /= k;

          if ((par.impMeasure == im_perm_breiman) || (par.impMeasure == im_perm_liaw)) {
            // store decrease of accuracy in vector
            // this approximation results in a more stable variance estimate
            if (k > 0) {
              variation /= k;
              // next, variation equals standard error
              // (a standard error calculation with div. by k - 1 would be biased)
              variation -= accur * accur;
              variation = sqrt(variation / k);
            }
          } else if ((par.impMeasure == im_perm_raw) || (par.impMeasure
              == im_perm_meng)) {
            variation = 1;
          }

          if ((k > 0) && (variation > 0)) {
            permVarImpPairVec.push_back(std::make_pair(
                accur / variation, *itIdx));
          } else {
            // -5 was suggested in original fortran code of breiman/cutler
            permVarImpPairVec.push_back(std::make_pair(0, *itIdx));
          }

          // restore column
          data.setCol(*itIdx, bakVec);

          // time estimate
          if (i < timeSpan + timeShift && i >= timeShift) {
            end = clock();
            timeEst += (double) (end - start) / CLOCKS_PER_SEC;
          }
          if (i == timeSpan - 1 + timeShift) {
            timeEst = (idx.size() / numOfThreads_omp - 1) * timeEst / timeSpan;
            *io.outVerbose << "Permutation time estimate: ";
            Helper::printTime(timeEst, *io.outVerbose);
            *io.outVerbose << std::endl;
          }

          // prepare for next iteration
          ++i;
          ++itIdx;
          if (itIdx == idx.end()) {
            keepRunning = false;
          } else {
            if (*itIdx == par.depVar)
              ++itIdx;
            if (itIdx == idx.end())
              keepRunning = false;
          }
          //*io.outVerbose << std::endl;
        }

        // end up with all sync.
#pragma omp barrier

      } while (keepRunning);

    }

    // only master thread remains
    // sort importance values
    sort(permVarImpPairVec.begin(), permVarImpPairVec.end());

    delete[] bakVec;
    delete[] vec;
  }

  /*
   * make permutation importance
   *
   * this is fast but inaccurate, because only 1 permutation per variable
   * is made.
   */
  static void makePermVarImp3(
      RJunglePar &par, RJungleIO &io, RJungleGen<T> &gen, DataFrame<T> &data,
      DataTreeSet &oobSet, std::vector<uli_t> *&colMaskVec, std::vector<
          CmpldTree<T> *> &cmpldTrees,
      std::vector<std::pair<double, uli_t> > &permVarImpPairVec) {

    T *vec;
    T *bakVec;
    std::vector<uli_t> idx;
    std::vector<std::pair<uli_t, uli_t> > numOfTrees;
    std::vector<double> accuracy;
    //std::vector<double > variation;
    std::vector<uli_t>::iterator itIdx;
    typename std::vector<CmpldTree<T> *>::iterator itCurTree;
    uli_t curTree;
    double dump, accur, variation, k;

    clock_t start = 0;
    clock_t end = 0;
    int timeSpan = 50; // maximal 50
    int timeShift = 10;
    double timeEst = 0;
    int totalTrees = 0;
    double percentRatio = 0.1;
    uli_t percentRate = 0;
    uli_t percent = 0;
    uli_t treeBufferSize;
    uli_t maxTreeBufferSize = 10;
    int i;

    // reset tree / accuracy size
    cmpldTrees.resize(par.ntree, NULL);
    accuracy.resize(par.ntree, -10);
    //variation.resize(ntree, 0);

    permVarImpPairVec.clear();

    *io.outVerbose << "Calculating permutation importance..." << std::endl;

    // create variable vector
    if (colMaskVec == NULL) {
      for (uli_t i = 0; i < par.ncol; ++i) {
        if (i == par.depVar)
          continue;
        idx.push_back(i);
      }
    } else {
      idx.assign(colMaskVec->begin(), colMaskVec->end());
    }

    // set time span
    timeSpan = 1 + std::min((uli_t) 100, (uli_t) idx.size() / 10);
    percentRate = std::max((uli_t) 1, (uli_t) floor(idx.size() * percentRatio));

    // permutation vectors
    vec = new T[par.nrow];
    bakVec = new T[par.nrow];

    // perform perm. imp.
    i = 0;
    itIdx = idx.begin();
    bool keepRunning = true;

    // create multiple threads
#pragma omp parallel default(shared)
    if (itIdx != idx.end()) {
      //DataFrame data(data);
      int ntree_omp;
      int numOfThreads_omp;
      std::vector<double> accuracy_omp;
      std::vector<double> variation_omp;
      std::vector<CmpldTree<T> *> cmpldTrees_omp;
      std::vector<Tree<T, uli_t> *> trees_omp;
      int startTree_omp = 0;
      int endTree_omp = 0;
      int id_omp;
      uli_t j_omp;
      double dumpAccur_omp, dumpVariation_omp, k_omp, curacc_omp, diffAccur_omp;

#ifdef _OPENMP
      numOfThreads_omp = omp_get_num_threads();
      id_omp = omp_get_thread_num();
#else
      numOfThreads_omp = 1;
      id_omp = 0;
#endif

#pragma omp master
      numOfTrees.resize(numOfThreads_omp, std::make_pair(0, 0));

#pragma omp barrier

      // balance trees over all threads
      ntree_omp = (int) floor(par.ntree / numOfThreads_omp);

      //#pragma omp critical
      //if (id_omp != 0) {totalTrees = 0; ntree_omp = 0;}

#pragma omp critical
      if (id_omp != 0) {
        startTree_omp = totalTrees;
        endTree_omp = totalTrees + ntree_omp;
        numOfTrees[id_omp] = std::make_pair(startTree_omp, endTree_omp);
        totalTrees += ntree_omp;
      }

#pragma omp barrier

#pragma omp master
      {
        ntree_omp = (int) par.ntree - totalTrees;
        startTree_omp = totalTrees;
        endTree_omp = totalTrees + ntree_omp;
        numOfTrees[id_omp] = std::make_pair(startTree_omp, endTree_omp);
        treeBufferSize = 1 + std::min((uli_t) floor(par.ntree
            / numOfThreads_omp), maxTreeBufferSize);

        *io.outVerbose << numOfThreads_omp << " thread(s) permuting "
            << par.ntree << " trees" << std::endl;
        /*
         << ntree << " trees (" << ntree_omp;

         for (uli_t j = 1; j < numOfTrees.size(); ++j) {
         *io.outVerbose << "+" << numOfTrees[j].second - numOfTrees[j].first;
         }

         *io.outVerbose << ")" << std::endl;
         */
        *io.outVerbose << "Get initial accuracy" << std::flush;
      }

#pragma omp barrier

#pragma omp critical
      {
        // get original accuracy of training set
        accuracy_omp.resize(par.ntree, -10);
        variation_omp.resize(par.ntree, 0);
        cmpldTrees_omp = cmpldTrees;
      }

      for (int j = startTree_omp; j < endTree_omp; ++j) {
        accuracy_omp[j] = gen.fct.getAccuracyOfCmpldTrees(
            &data, cmpldTrees_omp, gen.fct, data.getDepVar(), oobSet, j);
      }

      // merging compiled trees
#pragma omp critical
      {
        for (int j = startTree_omp; j < endTree_omp; ++j) {
          accuracy[j] = accuracy_omp[j];
        }
        *io.outVerbose << "." << std::flush;
      }

#pragma omp barrier

      // each thread gets all info
#pragma omp critical
      {
        for (int j = 0; j < (int) accuracy.size(); ++j) {
          accuracy_omp[j] = accuracy[j];
        }
      }

#pragma omp barrier

#pragma omp master
      *io.outVerbose << std::endl << "Permuting..." << std::endl;

      // start with all sync.
      do {
        // master permutes next column
#pragma omp master
        {
          //*io.outVerbose << *itIdx << ":" << std::flush;
          // management for time estimation
          if (i < timeSpan + timeShift && i >= timeShift)
            start = clock();
          if ((uli_t) i == percent + percentRate) {
            percent = percent + percentRate;
            *io.outVerbose << "progress: " << (uli_t) (100.0 * (double) percent
                / idx.size()) << "%" << std::endl;
          }

          // permutate current column
          data.getCol(*itIdx, vec);
          memcpy(bakVec, vec, par.nrow * sizeof(T));
          gsl_ran_shuffle(par.rng, vec, par.nrow, sizeof(T));
          data.setCol(*itIdx, vec);

          curTree = 0;
          itCurTree = cmpldTrees.begin();
          dump = 0;
          accur = 0;
          variation = 0;
          k = 0;

        }

        // all are waiting for all
#pragma omp barrier

        // reset variables
        dumpAccur_omp = 0;
        dumpVariation_omp = 0;
        k_omp = 0;

        while (itCurTree != cmpldTrees.end()) {
          // read some trees in buffer
#pragma omp critical
          {
            j_omp = 0;
            startTree_omp = curTree;
            while (itCurTree != cmpldTrees.end() && j_omp < treeBufferSize) {
              itCurTree++;
              j_omp++;
              curTree++;
            }
            endTree_omp = curTree;
          }

          // calc. tree accuracy / variation
          // add it to global accuracy / variation
          // all calculate importance
          // get decrease of accuracy
          for (int j = startTree_omp; j < endTree_omp; ++j) {
            // Does tree contains the variable?
            if (!(cmpldTrees_omp[j]->hasVarID(*itIdx)))
              continue;

            // get new accuracy
            curacc_omp = gen.fct.getAccuracyOfCmpldTrees(
                &data, cmpldTrees_omp, gen.fct, data.getDepVar(), oobSet, j);

            ++k_omp;

            // sum of differences
            diffAccur_omp = accuracy_omp[j] - curacc_omp;
            dumpAccur_omp += diffAccur_omp;
            // sum of squares

            if (par.impMeasure == im_perm_breiman) {
              // calculation of imp. measure in random forests (Breiman, Cutler)
              dumpVariation_omp += pow(diffAccur_omp, 2);
            } else if (par.impMeasure == im_perm_liaw) {
              // due to calculation of imp. measure in R's randomForest (Liaw, Wiener)
              dumpVariation_omp += pow(diffAccur_omp, 2)
                  * oobSet.getSmplSize(j);
            }
          }
          //#pragma omp barrier
        }

#pragma omp critical
        if (k_omp > 0) {
          accur += dumpAccur_omp;
          variation += dumpVariation_omp;
          k += k_omp;
          //*io.outVerbose << "." << std::flush;
        }

        // end up with all sync.
#pragma omp barrier

        // master restores content of column and prepares for next iteration
#pragma omp master
        {
          if ((par.impMeasure == im_perm_breiman) || (par.impMeasure == im_perm_liaw)
              || (par.impMeasure == im_perm_raw)) {
            k = par.ntree;
          }

          if (k > 0)
            accur /= k;

          if ((par.impMeasure == im_perm_breiman) || (par.impMeasure == im_perm_liaw)) {
            // store decrease of accuracy in vector
            // this approximation results in a more stable variance estimate
            if (k > 0) {
              variation /= k;
              // next, variation equals standard error
              // (a standard error calculation with div. by k - 1 would be biased)
              variation -= accur * accur;
              variation = sqrt(variation / k);
            }
          } else if ((par.impMeasure == im_perm_raw) || (par.impMeasure
              == im_perm_meng)) {
            variation = 1;
          }

          if ((k > 0) && (variation > 0)) {
            permVarImpPairVec.push_back(std::make_pair(
                accur / variation, *itIdx));
          } else {
            // -5 was suggested in original fortran code of breiman/cutler
            permVarImpPairVec.push_back(std::make_pair(0, *itIdx));
          }

          // restore column
          data.setCol(*itIdx, bakVec);

          // time estimate
          if (i < timeSpan + timeShift && i >= timeShift) {
            end = clock();
            timeEst += (double) (end - start) / CLOCKS_PER_SEC;
          }
          if (i == timeSpan - 1 + timeShift) {
            timeEst = (idx.size() / numOfThreads_omp - 1) * timeEst / timeSpan;
            *io.outVerbose << "Permutation time estimate: ";
            Helper::printTime(timeEst, *io.outVerbose);
            *io.outVerbose << std::endl;
          }

          // prepare for next iteration
          ++i;
          ++itIdx;
          if (itIdx == idx.end()) {
            keepRunning = false;
          } else {
            if (*itIdx == par.depVar)
              ++itIdx;
            if (itIdx == idx.end())
              keepRunning = false;
          }
          //*io.outVerbose << std::endl;
        }

        // end up with all sync.
#pragma omp barrier

      } while (keepRunning);

    }

    // only master thread remains
    // sort importance values
    sort(permVarImpPairVec.begin(), permVarImpPairVec.end());

    delete[] bakVec;
    delete[] vec;
  }

  /*
   * print intrinsic importance of variables in tree nodes
   */
  static void printIVarImp(
      RJunglePar &par, RJungleIO &io, DataFrame<T> &data, std::vector<
          std::pair<double, uli_t> > &freqPairs, uli_t iteration = 0) {

    Helper::printTreePair<T>(
        freqPairs, data, par.ntree, iteration, *io.outImportance);
  }

  static void printPermVarImp(
      RJunglePar &par, RJungleIO &io, DataFrame<T> &data, std::vector<
          std::pair<double, uli_t> > &freqPairs, uli_t iteration = 0) {

    //makePermVarImp(data, NULL, this->permVarImpPairVec);
    Helper::printTreePair<T>(
        freqPairs, data, par.ntree, iteration, *io.outImportance);
  }

};

#endif /* RJUNGLEIMPORTANCE_H_ */
