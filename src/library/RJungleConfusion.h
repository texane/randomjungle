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

#ifndef RJUNGLECONFUSION_H_
#define RJUNGLECONFUSION_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include <ctime>

#include "RJungleIO.h"
#include "RJungleGen.h"
#include "DataFrame.h"
#include "DataTreeSet.h"
#include "CmpldTree.h"
#include "Helper.h"
#include "treedefs.h"

template<class T>
class RJungleConfusion {
public:
  RJungleConfusion();
  virtual ~RJungleConfusion();

  static void printConfMatNew(RJungleIO &io, RJungleGen<T> &gen,
      DataFrame<T> &data, T* pred, DataTreeSet &oobSet, uli_t iteration = 0,
      uli_t numOfVars = 0) {

    *io.outVerbose << "Writing accuracy information..." << std::endl;
    *io.outVerbose << "Calculating confusion matrix..." << std::endl;
    *io.outConfusion << "Iteration: " << iteration << std::endl;
    *io.outConfusion << "Number of variables: " << numOfVars << std::endl;

    /*
     *io.outConfusion << "Training set: " << std::endl;
     */

    // print header
    if (io.isClassTree(data.par) && data.par.votes_flag) {
      *io.outVotes << "margin";
      for (size_t idx = 0; idx < data.varClassIndexer[data.par.depVar].size(); ++idx) {
        *io.outVotes << data.par.delimiter
            << data.varClassIndexer[data.par.depVar][idx];
      }
      *io.outVotes << std::endl;
    }

    if (io.isClassTree(data.par) && (iteration == 0)) {
      // header of second confusion file
      *io.outConfusion2 << "iteration";
      *io.outConfusion2 << data.par.delimiter << "numOfVars";
      for (size_t j = 0; j < data.depValVec.size(); ++j)
        *io.outConfusion2 << data.par.delimiter << (double) data.depValVec[j];
      *io.outConfusion2 << data.par.delimiter << "error";
      *io.outConfusion2 << std::endl;
    }

    if (io.isClassTree(data.par))
      *io.outConfusion2 << iteration << data.par.delimiter << numOfVars
          << data.par.delimiter;

    *io.outConfusion << "Test/OOB set: " << std::endl;

    gen.fct.getAccuracyOfCmpldTrees(&data, pred, data.getDepVar(), oobSet, -1,
        false, true, io);

  }

};

#endif /* RJUNGLECONFUSION_H_ */
