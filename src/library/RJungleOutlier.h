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

#ifndef RJUNGLEOUTLIER_H_
#define RJUNGLEOUTLIER_H_

#include <math.h>

#include "Helper.h"
#include "RJungleIO.h"
#include "DataFrame.h"
#include "DataTreeSet.h"
#include "Proximities.h"

class RJungleOutlier {
public:

	RJungleOutlier() {}

	virtual ~RJungleOutlier() {}

  // get outliers
  template<class T>
  static void getOutlier(
    RJungleIO &io, DataFrame<T> &data, Proximities<T> &proxi) {

    size_t i, j;
    std::vector<double> outlierScore(data.getnrow(), 0);

    *io.outVerbose << "Getting outlier score..." << std::endl;

    // get P^hat(i)
    for (i = 0; i < data.getnrow(); ++i)
      for (j = 0; j < data.getnrow(); ++j) {
        // what if i==j ?
        if ((data.at(i, data.par.depVar) != data.at(j, data.par.depVar)))
          continue;

        outlierScore[i] += pow(proxi.at(i, j), 2);
      }

    // get Score
    for (i = 0; i < data.getnrow(); ++i) {
      // outlier definition described on breiman's website
      // and copycat in paper "Using random forest for reliable classification and
      // cost-sensitive learning for medical diagnosis" by Yang et al.
      outlierScore[i] = data.varCategories[data.par.depVar][data.depIdxVec[i]]
        / ((outlierScore[i] > 0 ? outlierScore[i] : 1));

      // outlier definition actually used in breiman's and liaw's code
      //outlierScore[i] = data.getnrow() / ((outlierScore[i] > 0
      //  ? outlierScore[i] : 1));
    }

    // additional normalization in liaw's code
    if (data.par.outlier == 2) {
      std::vector<double> tempVec;
      double temp1, temp2;
      temp1 = temp2 = 0;

      // iterate over all groups
      for (j = 0; j < data.varCategories[data.par.depVar].size(); ++j) {
        tempVec.clear();

        // get outlier scores of group j
        for (i = 0; i < data.getnrow(); ++i)
          if (data.depIdxVec[i] == j)
            tempVec.push_back(outlierScore[i]);

        temp1 = Helper::getMedian<double>(tempVec);
        temp2 = Helper::getMad<double>(tempVec);

        // normalize scores
        if (temp2 == 0)
          temp2 = 1;

        for (i = 0; i < data.getnrow(); ++i)
          if (data.depIdxVec[i] == j)
            outlierScore[i] = (outlierScore[i] - temp1) / temp2;
      }
    }

    // write to file
    *io.outOutlier << "sample_ID" << data.par.delimiter << "outlier_score"
      << std::endl;

    for (i = 0; i < data.getnrow(); ++i) {
      *io.outOutlier << i << data.par.delimiter << outlierScore[i] << std::endl;
    }
  }
};

#endif /* RJUNGLEOUTLIER_H_ */
