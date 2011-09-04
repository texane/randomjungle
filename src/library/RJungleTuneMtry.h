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

#ifndef RJUNGLETUNEMTRY_H_
#define RJUNGLETUNEMTRY_H_

#include "Prediction.h"

template<class T>
class RJungleTuneMtry {
public:
  RJungleTuneMtry();
  virtual ~RJungleTuneMtry();

  /*
   * Tune the mtry parameter
   */
  static void tuneMtry(RJunglePar &_par, RJungleIO &_io, RJungleGen<T> &_gen,
      DataFrame<T> &data, std::vector<uli_t> *colMaskVec) {

    // pre conditions
    if (_par.ntree == 0)
      throw Exception(ERRORCODE_32);

    // declare

    //binder
    RJungleBinder<T> binder;
    binder.par = _par;
    binder.io = _io;
    binder.gen = _gen;

    RJunglePar &par = binder.par;
    RJungleIO &io = binder.io;
    RJungleGen<T> &gen = binder.gen;
    DataTreeSet &oobSet = binder.oobSet;
    std::vector<CmpldTree<T> *> &cmpldTrees = binder.cmpldTrees;
    Importance<T> *&intrinsicImportance = binder.intrinsicImportance;
    PermImportance<T> *&permImportance = binder.permImportance;
    Prediction<T> &pred = binder.pred;

    double pe;
    uli_t i, mtry;
    double theta;
    std::vector<uli_t> mtrys;

    // define
    mtry = par.ncol;
    theta = par.tunemtry;
    pe = 1.0;
		uli_t myTreeSize = par.ntree;
#ifdef HAVE_MPI
			if (par.mpiId == 0)
				myTreeSize = par.ntreeMpi;
#endif

    // importance
    intrinsicImportance = gen.fct.newImportanceObject();
    intrinsicImportance->setBaseData(&par, &io, &data);

    permImportance = new PermImportance<T> ();
    permImportance->setBaseData(&par, &io, &gen, &data);

    // create vector containing different mtry values
    Helper::getSlicedVec(par.ncol - 1, mtrys, theta);

    *io.outTuneMtry << "mtry" << par.delimiter << "OOB_prediction_error"
        << std::endl;

    // mtry loop
    for (i = 0; i < mtrys.size(); ++i) {
			// reset
			oobSet.init(par.nrow, myTreeSize);

			Helper::clearPtrVec<CmpldTree<T> *>(cmpldTrees);
			cmpldTrees.resize(myTreeSize, NULL);

      // importance
      intrinsicImportance->init();
      permImportance->init();

      // grow forest
      par.mtry = mtrys[i];

      *io.outVerbose << "Testing mtry: " << par.mtry << std::endl;

      // grow forest
      RJungleGrow<T>::growLocal(data, binder, colMaskVec, i);

      // output results
      *io.outConfusion << "mtry: " << par.mtry << std::endl;
      *io.outConfusion << "Test/OOB set: " << std::endl;

      // get oob error
      pe = 1 - gen.fct.getAccuracyOfCmpldTrees(&data, pred.pred,
          data.getDepVar(), oobSet, -1, false, false, io);

      RJungleConfusion<T>::printConfMatNew(io, gen, data, pred.pred, oobSet, i,
          par.ncol);

      *io.outConfusion << std::endl;

      *io.outTuneMtry << par.mtry << par.delimiter << pe << std::endl;

      *io.outVerbose << std::endl;

      // importance
      if (par.impMeasure > im_intrinsic) {
        permImportance->setIteration(i);
        permImportance->save();
      }
      intrinsicImportance->setIteration(i);
      intrinsicImportance->save();
    }

    // Burning down the forest
    Helper::clearPtrVec<CmpldTree<T> *>(cmpldTrees);
  }

};

#endif /* RJUNGLETUNEMTRY_H_ */
