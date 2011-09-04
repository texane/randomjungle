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

#ifndef RJUNGLEIMPUTE_H_
#define RJUNGLEIMPUTE_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include <ctime>

#include "config.h"
#include "RJunglePar.h"
#include "RJungleIO.h"
#include "RJungleGen.h"
#include "RJungleHelper.h"
#include "Proximities.h"
#include "DataFrame.h"
#include "DataTreeSet.h"
#include "Helper.h"
#include "treedefs.h"


template <class T >
class RJungleImpute {
public:
  RJungleImpute();
  virtual ~RJungleImpute();

  static void imputeData(
      RJunglePar &par,
      RJungleIO &io,
      DataFrame<T > &data,
      std::vector<uli_t > *&colMaskVec,
      Proximities<T > &proxi
    ) {

    std::vector<std::pair<uli_t, uli_t> > numOfTrees;
    typename std::vector<CmpldTree<T > *>::iterator itCurTree;
    std::vector<double > ldVec;
    std::vector<T > row, col, colMiss;
    unsigned int i, j;
    std::pair<T, double > res, res2;
    double val;
    double val2;
    clock_t start = 0;
    clock_t end = 0;
    unsigned int timeSpan = 2;
    unsigned int timeShift = 1;
    double timeEst = 0;


    //*io.outVerbose << "data" << data << std::endl;
    *io.outVerbose << "Imputing missings values..." << std::endl;

    // replace missings
    for(j = 0; j < par.ncol; j++) {

      if (j == timeShift) start = clock();

      //data.getColNoMissingsExt(j, col);
      if(!data.par.gwa_flag && (data.lom[j] != lom_nominal))
        data.getColMiss(j, colMiss);
      else
        data.getCol(j, col);

      // in GWA case: get variables (SNPs) in high LD with current SNP
      if (data.par.gwa_flag && (j != data.par.depVar))
        Helper::getSNPsInLD(data, j, colMaskVec, ldVec);

      for(i = 0; i < par.nrow; i++)
        if (data.missings.at(i, j)) {
          if(data.par.gwa_flag) {
            data.getRow(i, row);
            res = Helper::getMostFreqLD<T >(
              i,
              row,
              par.missingcode,
              ldVec);
            val = res.first;
            //std::cout << "ld val " << val << std::endl;
            //std::cout << "ld rel freq " << res.second << std::endl;

            res2 = Helper::getMostFreqProx2<T >(
              i,
              col,
              par.missingcode,
              proxi);
            val2 = res2.first;
            //std::cout << "proxi val " << val2 << std::endl;
            //std::cout << "proxi rel freq " << res2.second << std::endl;

            val = (val * res.second + val2 * res2.second);
            val /= res.second + res2.second;
            #ifndef HAVE_ROUND
            data.setArray(i, j, (T)round(val));
            #else
            data.setArray(i, j, (T)floor(val));
            #endif
          } else
            if ((data.lom[j] != lom_numeric) && !Helper::isGuessingLomNumeric(data.par)) {
              val = Helper::getMostFreqProx<T >(
                i,
                colMiss,
                par.missingcode,
                proxi);
              data.setArray(i, j, (T)(val));
              #ifndef HAVE_ROUND
              //data.setArray(i, j, (T)round(val));
              #else
              //data.setArray(i, j, (T)floor(val));
              #endif
            } else {
              val = Helper::getMeanProx<T >(
                i,
                colMiss,
                par.missingcode,
                proxi);
              val2 = Helper::getMeanProxX<T >(
                i,
                colMiss,
                par.missingcode,
                proxi);
              if (val2 != 0) {
                data.setArray(i, j, (T)(val / val2));
              } else {
                data.setArray(i, j, (T)(Helper::getMean2<T>(colMiss, par.missingcode)));
              }

              #ifndef HAVE_ROUND
              //data.setArray(i, j, (T)round(val / val2));
              #else
              //data.setArray(i, j, (T)floor(val / val2));
              #endif
            }
        }

      if (j == (timeSpan + timeShift)) {
        end = clock();
        timeEst = (double)(end - start) / CLOCKS_PER_SEC;
        timeEst = (par.nrow - j) * timeEst / timeSpan;
        *io.outVerbose << "Time estimate:";
        Helper::printTime(timeEst, *io.outVerbose);
        *io.outVerbose << std::endl;
      }


    }

    //std::cout << proxi << std::endl;

  }

  // impute missings with median values
  static void imputeTrivially(DataFrame<T > &data) {
    unsigned int i, j;
    std::vector<T > col;
    T val;

    for (j = 0; j < data.getncol(); ++j) {

      data.getColNoMissings(j, col);

      if (data.lom[j] != lom_numeric)
        val = Helper::getMostFreqX(col, data.getMissingCode());
      else
        val = Helper::getMedianX2(col);

      for (i = 0; i < data.getnrow(); ++i)
        if (data.at(i, j) == data.getMissingCode())
          data.setArray(i, j, val);
    }
  }

  static void imputeTriviallyGWA(RJunglePar &par, DataFrame<T > &data) {
    unsigned int i, j;
    std::vector<T > col;
    typename std::map<T, unsigned int > classMap;
    typename std::map<T, unsigned int >::iterator itMap;
    double p, pop;
    T val;

    for (j = 0; j < data.getncol(); ++j) {

      data.getColNoMissings(j, col);

      getClassesAndRate(col, classMap);

      pop = classMap[2] + classMap[1] + classMap[0];
      p   = classMap[2] + classMap[1] / 2;
      p  /= pop;

      for (i = 0; i < data.getnrow(); ++i)
        if (data.at(i, j) == data.getMissingCode()) {
          val = 0;
          if (gsl_rng_uniform(par.rng) < p) val++;
          if (gsl_rng_uniform(par.rng) < p) val++;
          data.setArray(i, j, val);
        }
    }
  }

  static void imputeSNPs(
      RJunglePar &par,
      RJungleIO &io,
      RJungleGen<T > &gen,
      DataFrame<T > &data,
      std::vector<uli_t > *colMaskFullVec
    ) {

    typename std::vector<CmpldTree<T > *>::iterator itCurTree;
    std::vector<std::pair<uli_t, uli_t> > numOfTrees;
    std::vector<double > ldVecFull, ldVec;
    std::vector<uli_t > *colMaskVec;
    std::pair<T, double > res, res2;
    unsigned int timeShift = 1;
    unsigned int timeSpan = 2;
    std::vector<T > row, col;
    unsigned int i, j, k, idx;
    Proximities<T > proxi;
    double timeEst = 0;
    clock_t start = 0;
    clock_t end = 0;
    double val, val2;
    bool wasNull = false;

    if (colMaskFullVec == NULL) {
      wasNull = true;

      colMaskFullVec = new std::vector<uli_t >();

      for(idx = 0; idx < par.ncol; idx++) {
        if (par.depVar != idx) colMaskFullVec->push_back(idx);
      }
    }

    Helper::printVec(*colMaskFullVec);

    colMaskVec = new std::vector<uli_t >();

    //Helper::printVec<uli_t >(*colMaskFullVec);
    //*io.outVerbose << "data" << data << std::endl;
    *io.outVerbose << "Imputing missings SNP values..." << std::endl;

    RJungleImpute<T >::imputeTrivially(data);

    for (k = 0; k < par.imputeIt; ++k) {
      // replace missings
      for (idx = 0; idx < colMaskFullVec->size(); ++idx) {

        if (idx == timeShift) start = clock();

        j = colMaskFullVec->at(idx);
        if (j == data.par.depVar) continue;

        // skip if no missings in column
        if (!data.anyMissingsInCol(j)) continue;

        // get data of current column
        data.getCol(j, col);

        // get SNPs in LD with current SNP
        Helper::getSNPsInLD(data, j, colMaskFullVec, ldVecFull);
        //std::cout << "j " << j << std::endl;

        // get best SNPs (r^2 > par.cutoffHighLD) but at least the best 4 SNPs
        Helper::getSNPsInHighLD(data, ldVecFull, ldVec, colMaskFullVec, colMaskVec);
        //Helper::printVec(*colMaskVec);

        // grow jungle now!
        RJungleGrow<T >::growForImputation(par, io, gen, data, colMaskVec, proxi);

        for(i = 0; i < par.nrow; i++) {
          if (data.missings.at(i, j)) {

            // get data of current row
            data.getRow(i, row);

            // guess value via LD structure
            res = Helper::getMostFreqLD<T >(i, row, par.missingcode, ldVec);
            val = res.first;

            // guess value via samples structure
            res2 = Helper::getMostFreqProx2<T >(i, col, par.missingcode, proxi);
            val2 = res2.first;

            /*
            std::cout << "i " << i << " j " << j << std::endl;
            std::cout << "ld val " << val << std::endl;
            std::cout << "ld rel freq " << res.second << std::endl;
            std::cout << "proxi val " << val2 << std::endl;
            std::cout << "proxi rel freq " << res2.second << std::endl;
            */

            val = (val * res.second + val2 * res2.second);
            val /= res.second + res2.second;

            #ifndef HAVE_ROUND
            data.setArray(i, j, round(val));
            #else
            data.setArray(i, j, floor(val));
            #endif
          }
        }

        if (idx == (timeSpan + timeShift)) {
          end = clock();
          timeEst = (double)(end - start) / CLOCKS_PER_SEC;
          timeEst = (par.nrow - idx) * timeEst / timeSpan;
          *io.outVerbose << "Time estimate:";
          Helper::printTime(timeEst, *io.outVerbose);
          *io.outVerbose << std::flush;
        }
      }
    }

    //data.printCSV(*io.outImputedData, colMaskFullVec);

    Helper::printVec(ldVec);
    Helper::printVec(ldVecFull);


    // clean up memory
    if (wasNull) {
      delete colMaskFullVec;
      colMaskFullVec = NULL;
    }



    delete colMaskVec;
  }

};

#endif /* RJUNGLEIMPUTE_H_ */
