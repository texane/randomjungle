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

#ifndef TIMPORTANCE_H_
#define TIMPORTANCE_H_

#include "treedefs.h"
#include "Helper.h"
#include "Importance.h"
#include "SaveCollector.h"
#include "Exception.h"

template<class T>
class TImportance: public Importance<T> {
public:
  TImportance() :
    Importance<T> () {
  }
  ;
  virtual ~TImportance() {
  }
  ;

  virtual void save() {
    saveCol.par = this->par;

    Helper::getRanks(this->amounts, saveCol.orderRow, true);
    Helper::seq<size_t>(ids, 0, this->amounts.size() - 1, 1);

    *this->io->outImportance << saveCol;
  }

  virtual inline void add(uli_t varID, double val) {
    amounts[varID] += val;
  }

  virtual void add(Importance<T> *impx) {
    TImportance<T> *imp = (TImportance<T>*) impx;
    for (uli_t i = 0; i < this->amounts.size(); ++i)
      this->amounts[i] += imp->amounts[i];
  }

  virtual void init() {
    this->saveCol.clear();

    // register vectors
    this->saveCol.push_back(&this->iterationVec, "iteration");
    this->saveCol.push_back(&this->ids, "id");
    this->saveCol.push_back(&this->data->varNames, "varname");

    this->saveCol.push_back(&this->amounts, "gini_index");

    // adjust this->saveCol
    this->saveCol.repeatLast = true;
    this->saveCol.showDepVar = false;
    this->saveCol.isAvailable.assign(this->par->ncol, 1);

    reset();
  }

  virtual void reset() {
    // reset vectors
    amounts.assign(this->par->ncol, 0);
  }

  static Importance<T>* newImportanceObject() {
    return new TImportance<T> ();
  }

  virtual void getBestVariables(size_t size, std::vector<uli_t> &outVec) {
    size_t i, j, idx;
    i = j = 0;

    Helper::getRanks(this->amounts, this->saveCol.orderRow, true);

    //evil conversion uli_t -> size_t
    outVec.clear();
    size_t len = this->saveCol.orderRow.size();
    while (i < len) {
      idx = this->saveCol.orderRow[i];

      if (j >= size)
        this->saveCol.isAvailable[idx] = 0;

      if ((this->saveCol.isAvailable[idx] == 1) && (j < size)) {
        outVec.push_back(idx);
        ++j;
      }

      ++i;
    }
  }

  void setIteration(uli_t it) {
    iterationVec.assign(1, it);

    if (it == 0)
      this->enableHeader();
    else
      this->disableHeader();
  }

  inline void enableHeader() {
    saveCol.showHeader = true;
  }

  inline void disableHeader() {
    saveCol.showHeader = false;
  }

  std::vector<std::pair<double, uli_t> > freqPairs;
  std::vector<double> amounts;

  std::vector<uli_t> iterationVec;
  std::vector<size_t> ids; // variable ids

  SaveCollector saveCol;
};

#endif /* TIMPORTANCE_H_ */
