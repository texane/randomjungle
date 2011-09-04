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

#ifndef RJUNGLECOMPILER_H_
#define RJUNGLECOMPILER_H_

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
class RJungleCompiler {
public:
  RJungleCompiler();
  virtual ~RJungleCompiler();

  static void compileTrees(
    RJunglePar &par, RJungleIO &io, RJungleGen<T> &gen, std::vector<Tree<T,
      uli_t> *> &trees, std::vector<CmpldTree<T> *> &cmpldTrees) {

    std::vector<std::pair<uli_t, uli_t> > numOfTrees;
    int totalTrees = 0;
    uli_t treeBufferSize;
    uli_t maxTreeBufferSize = 10;

    *io.outVerbose << "Compiling trees" << std::flush;

    // burn old compiled jungle
    typename std::vector<CmpldTree<T> *>::iterator it;
    it = cmpldTrees.begin();
    while (it != cmpldTrees.end()) {
      if (*it != NULL) {
        delete *it;
        *it = NULL;
      }
      ++it;
    }
    cmpldTrees.clear();
    cmpldTrees.resize(par.ntree, NULL);

    // compile new jungle
#pragma omp parallel default(shared)
    {
      int ntree_omp;
      int numOfThreads_omp;
      std::vector<CmpldTree<T> *> cmpldTrees_omp;
      std::vector<Tree<T, uli_t> *> trees_omp;
      int startTree_omp = 0;
      int endTree_omp = 0;
      int id_omp;

#ifdef _OPENMP
      numOfThreads_omp = omp_get_num_threads();
      id_omp = omp_get_thread_num();
#else
      numOfThreads_omp = 1;
      id_omp = 0;
#endif

#pragma omp master
      {
        numOfTrees.resize(numOfThreads_omp, std::make_pair(0, 0));
      }

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
        treeBufferSize = 1 + std::min(
          par.ntree / numOfThreads_omp, maxTreeBufferSize);
      }

#pragma omp critical
      {
        cmpldTrees_omp.resize(par.ntree, NULL);
        trees_omp.assign(trees.begin(), trees.end());
      }

      for (int j = startTree_omp; j < endTree_omp; ++j) {
#pragma omp critical
        cmpldTrees_omp[j] = gen.fct.cmplTree(*trees_omp[j]);
      }

#pragma omp critical
      {
        for (int j = startTree_omp; j < endTree_omp; ++j) {
          cmpldTrees[j] = cmpldTrees_omp[j];
        }
        *io.outVerbose << "." << std::flush;
      }
    }

    *io.outVerbose << std::endl;
  }

};

#endif /* RJUNGLECOMPILER_H_ */
