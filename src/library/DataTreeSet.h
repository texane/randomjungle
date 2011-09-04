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

#ifndef DATATREESET_H_
#define DATATREESET_H_

/*
 * Includes
 */

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <boost/dynamic_bitset.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#include <string>
#endif

#include "BitMatrix.h"

/*
 * DataTreeSet
 */

class DataTreeSet: public BitMatrix {
public:
	DataTreeSet() :
		BitMatrix() {
	}
	;

	DataTreeSet(uli_t nrow, uli_t ncol) :
		BitMatrix(nrow, ncol) {
	}
	;

	/**
	 * Copy constructor.
	 * @param dataFrame
	 */
	DataTreeSet(const DataTreeSet &set) {
		this->dep = set.dep;
		this->nrow = set.nrow;
		this->ncol = set.ncol;
	}

	inline uli_t getNsmpl() {
		return getNrow();
	}
	inline uli_t getNtree() {
		return getNcol();
	}
	inline uli_t getSmplSize(unsigned int j) {
		return getSizeOfCol(j);
	}

	virtual void combineMpi(RJunglePar &par) {
#ifdef HAVE_MPI
		MPI_Barrier(MPI_COMM_WORLD);

		if (par.mpiId == 0) {
			// master
			MPI_Status status;
			int nslaves;
			int slaveId;
			std::vector<std::string > strBitMatrix(par.mpiSize, std::string());
			BitMatrix bitMatrix;

			// add master's data
			boost::to_string(this->dep, strBitMatrix[0]);

			// receive bitmatrixes
			for (nslaves = 0; nslaves < par.mpiSize - 1; nslaves++) {
				// receive id
				MPI_Recv(&slaveId, 1, MPI_INT, MPI_ANY_SOURCE, MPI_TAG_DATATREESET_SLAVEID, MPI_COMM_WORLD, &status);

				// receive bitmatrix
				strBitMatrix[(size_t)(slaveId)] = bitMatrix.receiveStrMpi(slaveId);
			}

			// combine bitmatrixes
			std::string str;
			for (nslaves = 0; nslaves < (int)strBitMatrix.size(); nslaves++)
			str.append(strBitMatrix[(size_t)nslaves]);

			dep = boost::dynamic_bitset<>(str);
		} else {
			// slave

			// send id
			MPI_Send(&par.mpiId, 1, MPI_INT, 0, MPI_TAG_DATATREESET_SLAVEID, MPI_COMM_WORLD);

			// send bitmatrix
			this->sendMpi();
		}

		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
};

inline std::ostream & operator<<(std::ostream &os,
		const DataTreeSet &dataTreeSet) {
	return dataTreeSet.print(os);
}

#endif /*DATATREESET_H_*/

