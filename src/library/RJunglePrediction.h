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

#ifndef RJUNGLEPREDICTION_H_
#define RJUNGLEPREDICTION_H_

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
#include "Prediction.h"


template<class T>
class RJunglePrediction {
public:
  RJunglePrediction();
  virtual ~RJunglePrediction();

  static void showPredictionCmpld(
      RJungleIO &io, DataFrame<T> &data, Prediction<T> &pred,
      DataTreeSet &dataSet, uli_t iteration = 0) {

    std::vector<T> classRes;
    std::vector<T> classVec;
    std::vector<T> colVec;
    uli_t j, m;
    T medianVote;

    *io.outVerbose << "Writing prediction ..." << std::endl;

    for (j = 0; j < pred.nrow; ++j) {
      if (j >= dataSet.getNsmpl())
        throw Exception(ERRORCODE_31);

      // forest
      classRes.clear();
      for (m = 0; m < dataSet.getNtree(); ++m)
        if (dataSet.at(j, m))
          classRes.push_back(pred.pred[j + data.getnrow() * m]);

      medianVote = Helper::getMostFreqProp<T>(
          data.par, classRes, data.getMissingCode(), &data, false, NULL, j);

      if (medianVote != data.getMissingCode())
        *io.outPrediction << (double) medianVote << std::endl;
      else
        *io.outPrediction << "NA" << std::endl;
    }

  }

};

#endif /* RJUNGLEPREDICTION_H_ */
