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

#ifndef PREDICTION_H_
#define PREDICTION_H_

#include <vector>

#include "treedefs.h"
#include "DataFrame.h"
#include "RJunglePar.h"
#include "RJungleIO.h"
#include "RJungleGen.h"

template<class T>
class Prediction {
public:
  Prediction() :
    pred(NULL), acc(NULL), nrow(0), ncol(0) {
  }
  ;

  virtual ~Prediction() {
    removePred();
  }
  ;

  void removePred() {
    // remove
    if (this->pred != NULL)
      delete[] this->pred;

    if (this->acc != NULL)
      delete[] this->acc;

    this->pred = NULL;
    this->acc = NULL;
  }

  void init(RJunglePar *par, RJungleIO *io, RJungleGen<T> *gen,
      DataFrame<T> *data) {

    setBaseData(par, io, data, gen);

    uli_t ntree;

    // master should contain all data in the end
#ifdef HAVE_MPI
    ntree = (par->mpiId == 0)?(par->ntreeMpi):(par->ntree);
#else
    ntree = par->ntree;
#endif

    initMatrix(data->getnrow(), ntree);

  }

  virtual void initMatrix(size_t nrow, size_t ntree) {
    // remove
    removePred();

    this->nrow = nrow;
    this->ncol = ntree;

    // allocate
    this->pred = new T[nrow * ntree];
    this->acc = new double[ntree];
  }

  void initOne(RJunglePar *par, RJungleIO *io, RJungleGen<T> *gen,
      DataFrame<T> *data) {

    setBaseData(par, io, data, gen);
    initMatrix(1, data->getnrow());
  }

  void setBaseData(RJunglePar *par, RJungleIO *io, DataFrame<T> *data,
      RJungleGen<T> *gen) {

    this->par = par;
    this->io = io;
    this->gen = gen;
    this->data = data;
  }

  void add(CmpldTree<T> *cmpldTree, DataTreeSet *oobSet, long int treeId) {
    uli_t j;
    T *classMeVec;

    // classify samples
    for (j = 0; j < data->getnrow(); ++j) {
      // get sample

      data->getRow(j, classMeVec);

      if (data->at(j, par->depVar) == data->getMissingCode())
        continue;

      // get predictiony
      pred[j + nrow * treeId] = (oobSet->at(j, treeId)) ? gen->fct.classifyCmpldTree(
          classMeVec, cmpldTree) : data->getMissingCode();
    }

    acc[treeId] = gen->fct.getAccuracyOfCmpldTrees(this->data, pred,
        this->data->getDepVar(), *oobSet, treeId, false, false, *this->io);
  }

  void addOne(CmpldTree<T> *cmpldTree, DataTreeSet *oobSet, long int treeId) {
    uli_t j;
    T *classMeVec;

    // classify samples
    for (j = 0; j < data->getnrow(); ++j) {

      // get sample
      data->getRow(j, classMeVec);

      if (data->at(j, par->depVar) == data->getMissingCode())
        continue;

      // get prediction
      pred[j] = (oobSet->at(j, treeId)) ? gen->fct.classifyCmpldTree(
          classMeVec, cmpldTree) : data->getMissingCode();
    }
  }

  inline T at(size_t sample, size_t tree) {
    return pred[sample + nrow * tree];
  }

  void combineMpi() {
#ifdef HAVE_MPI
    // get data type
    MPI_Datatype dataType;

#ifdef __SPARSE_DATA__
    dataType = MPI_CHAR;
#else
    switch (par->memMode) {
      case 0:
      dataType = MPI_DOUBLE;
      break;
      case 1:
      dataType = MPI_FLOAT;
      break;
      case 2:
      dataType = MPI_CHAR;
      break;
    }
#endif

    MPI_Barrier(MPI_COMM_WORLD);

    *this->io->outVerbose << "Combining internal prediction..." << std::endl;

    if (par->mpiId == 0) {
      // master
      int slaveId;
      int length;
      int totalLength = par->ntree * par->nrow; // do not overwrite master's data

      MPI_Status status;
      
      // receive data from slaves
      for (slaveId = 1; slaveId < par->mpiSize; slaveId++) {
        // receive length
        MPI_Recv(&length, 1, MPI_INT, slaveId, MPI_TAG_PREDICTION_LEN, MPI_COMM_WORLD, &status);

        // receive pred data
        MPI_Recv(pred + totalLength, length, dataType, slaveId, MPI_TAG_PREDICTION_PRED, MPI_COMM_WORLD, &status);

        // print status
        //*this->io->outVerbose << slaveId << "/" << (par->mpiSize - 1) << std::endl;

        totalLength += length;
      }
      *this->io->outVerbose << "Received internal prediction" << std::endl;

    } else {
      //slave

      int length = nrow * ncol;

      MPI_Send(&length, 1, MPI_INT, 0, MPI_TAG_PREDICTION_LEN, MPI_COMM_WORLD);

      MPI_Send(pred, nrow * ncol, dataType, 0, MPI_TAG_PREDICTION_PRED, MPI_COMM_WORLD);

      // print status
      *this->io->outVerbose << "Sent internal prediction" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  T* pred;
  double* acc;
  size_t nrow, ncol;

  RJunglePar *par;
  RJungleIO *io;
  RJungleGen<T> *gen;
  DataFrame<T> *data;

};

#endif /* PREDICTION_H_ */
