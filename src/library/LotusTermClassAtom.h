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

#ifndef LOTUSTERMCLASSATOM_H_
#define LOTUSTERMCLASSATOM_H_


#include <vector>
#include <limits>
#include "TermClassAtom.h"
#include "treedefs.h"

template <class T, class C >
class LotusTermClassAtom : public TermClassAtom<T, C > {
public:
	LotusTermClassAtom(): val(T()) { };
	LotusTermClassAtom(const LotusTermClassAtom<T, C > &termCA) : val(termCA.val) { };
	LotusTermClassAtom(T val) : val(val) { };

	virtual ~LotusTermClassAtom(void) {
	};

	/// waste
	T val;

	/// x predictor variable
	uli_t varID;

	/// Betas of logistic regression
  std::vector<double > betas;
};

#endif /*LOTUSTERMCLASSATOM_H_*/
