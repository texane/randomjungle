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

#ifndef RJUNGLEPROTO_H_
#define RJUNGLEPROTO_H_

class ProtCount {
public:
  ProtCount(size_t numOfSamples, uli_t row) :
    numOfSamples(numOfSamples), row(row) {
  }
  ;

  ProtCount() :
    numOfSamples(0), row(0) {
  }
  ;

  void loadUp(ProtCount &prot) {
    if (prot.numOfSamples > numOfSamples) {
      numOfSamples = prot.numOfSamples;
      row = prot.row;
    }
  }

  size_t numOfSamples;
  uli_t row;
};

template<class T>
class RJungleProto {
public:
  RJungleProto();
  virtual ~RJungleProto();

  static void getPrototypes(
    RJungleIO &io, RJungleGen<T> &gen, DataFrame<T> &data,
    std::vector<uli_t> *colMaskVec, Proximities<T> &proxi) {

    size_t i, j, k, numOfNeighbours;
    std::vector<size_t> ranks;
    std::vector<double> vec;
    std::vector<ProtCount> protCounts, protCountsNext;

    numOfNeighbours = data.par.prototypes;

    // init protCounts
    for (i = 0; i < data.varClassIndexer[data.par.depVar].size(); ++i) {
      protCounts.push_back(ProtCount(0, 0));
    }

    // for each sample
    for (i = 0; i < data.par.nrow; ++i) {
      // get ranks
      vec = proxi.getRow(i);
      ranks = Helper::getRanks<double>(vec, false);

      // prune ranks
      ranks.erase(ranks.begin() + numOfNeighbours, ranks.end());

      // count number of samples of each class
      countNumOfSamples(ranks, i, data, protCountsNext);

      // store results
      for (j = 0; j < protCounts.size(); ++j)
        protCounts[j].loadUp(protCountsNext[j]);
    }

    // print header
    size_t num = (colMaskVec != NULL) ? colMaskVec->size() : data.par.ncol;
    *io.outPrototypes << "class";
    for (size_t j = 0; j < num; ++j) {
      k = (colMaskVec != NULL) ? colMaskVec->at(j) : j;
      if (k == data.par.depVar)
        continue;
      *io.outPrototypes << data.par.delimiter << data.varNames[k];
    }
    *io.outPrototypes << std::endl;

    // get neighbours
    for (j = 0; j < protCounts.size(); ++j) {
      i = protCounts[j].row;
      *io.outPrototypes << data.depValVec[j];

      // get ranks
      vec = proxi.getRow(i);
      ranks = Helper::getRanks<double>(vec, false);

      // prune ranks
      ranks.erase(ranks.begin() + numOfNeighbours, ranks.end());

      // get x data
      for (i = 0; i < num; ++i) {
        k = (colMaskVec != NULL) ? colMaskVec->at(i) : i;
        if (k == data.par.depVar)
          continue;

        data.getRowMaskedColNoMissing3(k, ranks, vec);
        *io.outPrototypes << data.par.delimiter << Helper::getMedianX2<double>(
          vec);
      }
      *io.outPrototypes << std::endl;
    }
  }

  static void countNumOfSamples(
    std::vector<size_t> &ranks, uli_t row, DataFrame<T> &data, std::vector<
      ProtCount> &protCounts) {

    size_t i;
    protCounts.clear();

    // init protCounts
    for (i = 0; i < data.varClassIndexer[data.par.depVar].size(); ++i) {
      protCounts.push_back(ProtCount(0, row));
    }

    // load up protCounts
    for (i = 0; i < ranks.size(); ++i) {
      protCounts[data.depIdxVec[ranks[i]]].numOfSamples++;
    }
  }

};

#endif /* RJUNGLEPROTO_H_ */
