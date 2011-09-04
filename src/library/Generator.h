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

#ifndef GENERATOR_H_
#define GENERATOR_H_

/*
 * Includes
 */

#include <iostream>
#include <cstring>
#include <vector>
#include <string>
#include <limits>
#include <ctime>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "Helper.h"
#include "CmplFct.h"
#include "treedefs.h"

/*
 * Def.
 */

template <class T >
class Generator {
public:
	// const. / dest.
	Generator() {
	}
	virtual ~Generator() {
	}

	// find the best variable in each node
virtual static GETBESTVARFCT(T) {

}

virtual static GETNODEPERFFCT(T) {

}

virtual static GETTERMRESFCT(T) {

}

virtual static CMPLTREEFCT(T) {

}

virtual static CLASSIFYCMPLDTREEFCT(T) {

}
}

#endif /* GENERATOR_H */
