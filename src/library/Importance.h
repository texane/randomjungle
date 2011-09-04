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

#ifndef IMPORTANCE_H_
#define IMPORTANCE_H_

#include <vector>

#include "treedefs.h"
#include "RJunglePar.h"
#include "RJungleIO.h"

template<class T>
class Importance {
public:
	Importance() {
	}
	;

	virtual ~Importance() {
	}
	;

	virtual void printHeader() {
	}
	virtual void print() {
	}
	virtual void save() {
	}

	void setBaseData(RJunglePar *par, RJungleIO *io, DataFrame<T> *data) {
		this->par = par;
		this->io = io;
		this->data = data;
	}

	void setBaseDataFrom(Importance<T> &imp) {
		this->par = imp.par;
		this->io = imp.io;
		this->data = imp.data;
	}

	virtual void add(Importance<T> *imp) {
	}

	virtual void init() {
	}
	virtual void reset() {
	}

	virtual void setIteration(uli_t iteration) {
		this->iteration = iteration;
	}
	;
	virtual uli_t getIteration() {
		return iteration;
	}
	;

	virtual void getBestVariables(size_t size, std::vector<uli_t> &outVec) {
	}
	;
	virtual void combineMpi() {
	}
	;

	uli_t iteration;
	RJunglePar *par;
	RJungleIO *io;
	DataFrame<T> *data;

};

#endif /* IMPORTANCE_H_ */
