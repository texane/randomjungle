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

#ifndef RJUNGLEBINDER_H_
#define RJUNGLEBINDER_H_

#include "Prediction.h"
#include "Importance.h"
#include "TImportance.h"
#include "PermImportance.h"
#include "Helper.h"
#include "CmplFct.h"
#include "Proximities.h"

template<class T>
class RJungleBinder {
public:
  RJungleBinder() :
    intrinsicImportance(NULL), permImportance(NULL) {
  }
  ;

  virtual ~RJungleBinder() {
    if (trees.size() != 0) {
      typename std::vector<Tree<T, uli_t> *>::iterator it;
      it = trees.begin();
      while (it != trees.end()) {
        if (*it != NULL)
          delete *it;
        ++it;
      }
      trees.clear();
    }

    if (cmpldTrees.size() != 0) {
      typename std::vector<CmpldTree<T> *>::iterator it;
      it = cmpldTrees.begin();
      while (it != cmpldTrees.end()) {
        if (*it != NULL)
          delete *it;
        ++it;
      }
      cmpldTrees.clear();
    }

    if (this->intrinsicImportance != NULL)
      delete intrinsicImportance;

    if (this->permImportance != NULL)
      delete permImportance;
  }

  RJunglePar par;
  RJungleIO io;
  RJungleGen<T> gen;

  std::vector<std::pair<double, uli_t> > freqPairs;
  std::vector<std::pair<double, std::pair<uli_t, uli_t> > > twoWayInteraction;

  // cache variables
  std::vector<Tree<T, uli_t> *> trees;
  std::vector<CmpldTree<T> *> cmpldTrees;
  Importance<T> *intrinsicImportance;
  PermImportance<T> *permImportance;
  Prediction<T> pred;

  // out of bag data set
  DataTreeSet oobSet;

  // sample proximities
  Proximities<T> proxi;

  // variable proximity
  Proximities<T> varProxi;
};

#endif /* RJUNGLEBINDER_H_ */
