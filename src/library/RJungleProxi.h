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

#ifndef RJUNGLEPROXI_H_
#define RJUNGLEPROXI_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include <ctime>

#include "RJunglePar.h"
#include "RJungleIO.h"
#include "RJungleGen.h"
#include "DataFrame.h"
#include "DataTreeSet.h"
#include "CmpldTree.h"
#include "Helper.h"
#include "treedefs.h"

template<class T>
class RJungleProxi {
public:
  RJungleProxi();
  virtual ~RJungleProxi();

  // finalize proximities
  static void finalize(Proximities<T> &proxi) {
    size_t i, j;

    // normalize proximities
    for (i = 0; i < proxi.data->par.nrow; ++i) {
      for (j = i + 1; j < proxi.data->par.nrow; ++j) {
        proxi.set(j, i, proxi.at(i, j));
      }

      proxi.set(i, i, 1);
    }
  }

  static void updateProximitiesCmpld(
      DataFrame<T> *data, CmpldTree<T> *cmpldTree, BuildinGenFct<T> &genFct,
      DataTreeSet &oobSet, std::vector<double> &oobpair, Proximities<T> &proxi,
      uli_t treeId) {

    T* classMeVec;

    uli_t i, j;
    uli_t nrow = data->getnrow();
    uli_t inf = std::numeric_limits<uli_t>::infinity();
    std::vector<T> classRes;
    std::vector<T> classVec;
    std::vector<T> colVec;

    std::vector<uli_t> smplTermNode(nrow, inf);

    // get terminal nodes
    for (j = 0; j < nrow; ++j) {
      //if (j >= oobSet.getNsmpl())
      //  throw Exception(ERRORCODE_31);

      if (data->at(j, data->getDepVar()) == data->getMissingCode())
        continue;

      data->getRow(j, classMeVec);

      //if (oobSet.at(j, treeId))
      smplTermNode[j] = genFct.classifyIDcmpldTree(classMeVec, cmpldTree);
    }

    // get proxi
    for (i = 0; i < nrow; ++i)
      for (j = i + 1; j < nrow; ++j) {
        //if ((smplTermNode[i] != inf) && (oobSet.at(i, treeId) && oobSet.at(j, treeId))) {
        //if (oobSet.at(i, treeId) && oobSet.at(j, treeId)) {

        // proxi. estimation using OOB
        if (oobSet.at(i, treeId) ^ oobSet.at(j, treeId)) {
          ++oobpair[i - 1 + j * (nrow - 1)];
          if (smplTermNode[i] == smplTermNode[j])
            proxi.add(i, j, 1);
        }
      }
  }

  static void normalizeProximitiesCmpld(
      DataFrame<T> *data, std::vector<double> &oobpair, Proximities<T> &proxi) {

    uli_t i, j;
    uli_t nrow = data->getnrow();
    double temp;
    i = j = 0;
    temp = 0;

    // normalize proximities
    for (i = 0; i < nrow; ++i) {
      for (j = i + 1; j < nrow; ++j) {
        temp = oobpair[i - 1 + j * (nrow - 1)];
        if (temp > 0)
          proxi.div(i, j, temp);
      }
    }
  }

  /*
   * Variable proximities of classification trees (CART)
   */

  static void initVarProx(DataFrame<T> &data, Proximities<T> &varProxi) {
    varProxi.data = &data;
    varProxi.setDim(data.par.ncol, data.par.ncol);
    varProxi.initMatrix();
    varProxi.reset();
  }

  // get variable proximities
  // make proximity 2nd order with top variables
  static void updateVarProximitiesCmpld(
      RJungleGen<T> &gen, DataFrame<T> &data, CmpldTree<T> *&cmpldTree,
      std::vector<uli_t> *&colMask, Proximities<T> &varProxi) {

    // add proximity values to data frame
    switch (data.par.varproximities) {
    case vp_SIMPLE:
      updateVarProximitiesCmpld1(&data, cmpldTree, gen.fct, colMask, varProxi);
      break;
    case vp_EXTEND1:
      updateVarProximitiesCmpld2(&data, cmpldTree, gen.fct, colMask, varProxi);
      break;
    default:
      ;
    }
  }

  static void finalizeVarProx(RJungleIO &io, DataFrame<T> &data, std::vector<
      uli_t> *colMask, Proximities<T> &varProxi) {

    std::vector<std::string> varNames;
    uli_t i;

    //*io.outVerbose << "Calculating variable proximity..." << std::endl;

    // init more vriables

    size_t ncol;
    if (colMask == NULL) {
      ncol = data.par.ncol;

      for (i = 0; i < ncol; ++i) {
        varNames.push_back(data.varNames[i]);
      }

    } else {
      ncol = colMask->size();

      for (i = 0; i < ncol; ++i) {
        varNames.push_back(data.varNames[colMask->at(i)]);
      }
    }

    // normalize proximities
    varProxi.div(data.par.ntree);
    varProxi.setDiag(1);

    for (i = 0; i < ncol; ++i) {
      if (i == data.par.depVar)
        continue;
      if (i > 0)
        *io.outVarProximity << data.par.delimiter;
      *io.outVarProximity << varNames[i];
    }
    *io.outVarProximity << std::endl;
    varProxi.printCSV(*io.outVarProximity);

  }

  /*
   * Are two variables in the same tree? Then add proximity by one.
   */
  static void updateVarProximitiesCmpld1(
      DataFrame<T> *data, CmpldTree<T> *cmpldTree, BuildinGenFct<T> &genFct,
      std::vector<uli_t> *colMask, Proximities<T> &proxi) {

    std::vector<uli_t> vars;
    size_t i, j;

    // get variables in tree
    if (colMask == NULL) {
      for (i = 0; i < data->par.ncol; ++i) {
        // Does tree contains the variable?
        if (cmpldTree->hasVarID(i))
          vars.push_back(i);
      }
    } else {
      for (i = 0; i < colMask->size(); ++i) {
        // Does tree contains the variable?
        if (cmpldTree->hasVarID(colMask->at(i)))
          vars.push_back(i);
      }
    }

    // update proximities
    for (i = 0; i < vars.size(); ++i) {
      for (j = i + 1; j < vars.size(); ++j) {
        proxi.add(vars[i], vars[j], 1);
        proxi.add(vars[j], vars[i], 1);
      }
    }
  }

  /*
   * Is variable i kid of variable j? Then add proximity by one.
   */
  static void updateVarProximitiesCmpld2(
      DataFrame<T> *data, CmpldTree<T> *cmpldTree, BuildinGenFct<T> &genFct,
      std::vector<uli_t> *colMask, Proximities<T> &proxi) {

    if (genFct.addVarProx == NULL)
      throw Exception(ERRORCODE_50);

    genFct.addVarProx(*data, proxi, cmpldTree, colMask);
  }

};

#endif /* RJUNGLEPROXI_H_ */
