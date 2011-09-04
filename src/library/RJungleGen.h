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

#ifndef RJUNGLEGEN_H_
#define RJUNGLEGEN_H_

/*
 * Includes
 */

#include "CmplFct.h"
#include "BuildinGenerators.h"
#include "DataFrame.h"
#include "RJunglePar.h"
#include "ClassAtom.h"
#include "CmpldTree.h"
#include "DataTreeSet.h"
#include "Tree.h"


/*
 * Def.
 */

template <class T >
class RJungleGen{
public:
  RJungleGen() {
  	// call macros to declare important generator functions
  	wasInitialized = false;
  }

  ~RJungleGen() {}

  void init(RJunglePar &par, DataFrame<T > &data) {
    wasInitialized = true;

    // define building functions
    if (par.treeType == tt_CART) {
      rjungle::Generator::
      CART::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CARTcontCont) {
      rjungle::Generator::
      CARTcontCont::assignFunctions<T >(&this->fct);
    } else if (par.treeType == tt_CARTcontCate) {
      rjungle::Generator::
      CARTcontCate::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_CARTcateCate) {
      rjungle::Generator::
      CARTcateCate::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_LOTUS) {
      rjungle::Generator::
      LOTUS::assignFunctions<T >(&this->fct);
  	} else if (par.treeType == tt_CARTsrt) {
      rjungle::Generator::
      CARTsrt::assignFunctions<T >(&this->fct);
  	} else {
      rjungle::Generator::
      CART::assignFunctions<T >(&this->fct);
  	}

    // check usability of trees regarding data
    if (fct.checkUsability != NULL) (*fct.checkUsability)(data);
  }

  BuildinGenFct<T > fct;

  bool wasInitialized;
};

#endif /* RJUNGLEGEN_H_ */
