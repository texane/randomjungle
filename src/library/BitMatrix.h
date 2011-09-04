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

#ifndef BITMATRIX_H_
#define BITMATRIX_H_

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

#ifndef NULL
#define NULL 0
#endif

/** 
 * \brief Representation of a bit matrix that can store bits exclusively.
 */
class BitMatrix {
public:
	/** 
	 * \brief Constructor that sets the data dimensions to 0.
	 */
	BitMatrix() :
		nrow(0), ncol(0) {
	}
	;

	/** 
	 * \brief Destructor
	 */
	virtual ~BitMatrix() {
	}
	;

	/** 
	 * \brief Constructor that sets the data dimensions.
	 * 
	 * @param nrow Number of rows.
	 * @param ncol Number of columns.
	 */
	BitMatrix(uli_t nrow, uli_t ncol) {
		init(nrow, ncol);
	}
	;

	/** 
	 * \brief Initializes the bit matrix with "false".
	 * Dimensions have to be given and old matrix data
	 * will be deleted.
	 * 
	 * @param nrow.
	 * @param ncol.
	 */
	void init(uli_t nrow, uli_t ncol) {
		this->nrow = nrow;
		this->ncol = ncol;
		dep.reset();
		dep.resize(nrow * ncol, false);
	}

	/** 
	 * \brief Sets all cells in bit matrix to "true"
	 */
	void setAll() {
		dep.set();
	}

	/** 
	 * \brief Sets all cells in bit matrix to "false"
	 */
	void setNone() {
		dep.reset();
	}

	/** 
	 * \brief Adds the given data to bit matrix (OR operator)
	 * 
	 * @param dataSet A bit matrix that should be added.
	 */
	inline void operator|=(const BitMatrix &dataSet) {
		dep |= dataSet.dep;
	}

	/** 
	 * \brief Returns logical value at a specific position
	 * in bit matrix.
	 * 
	 * @param row Row position in bit matrix.
	 * @param col Column position in bit matrix.
	 * 
	 * @return Logical value at that specific position
	 */
	inline bool at(uli_t row, uli_t col) const {
		return dep[row + col * nrow];
	}

	/** 
	 * \brief Sets logical value at a specific position.
	 * 
	 * @param row Row position in bit matrix.
	 * @param col Column position in bit matrix.
	 * @param val Value to be set.
	 */
	inline void at(uli_t row, uli_t col, bool val) {
		dep[row + col * nrow] = val;
	}

	/** 
	 * \brief Prints the content of bit matrix to output
	 * stream.
	 * 
	 * @param os Output stream.
	 * 
	 * @return Output stream
	 */
	std::ostream &print(std::ostream &os) const {
		uli_t i, j;

		for (i = 0; i < this->nrow; i++) {
			for (j = 0; j < this->ncol; j++) {
				os << this->dep[i + j * nrow] << "\t";
			}
			os << std::endl;
		}
		return os;
	}

	/** 
	 * \brief Returns the sum of row i.
	 * 
	 * @param i Index of row.
	 * 
	 * @return Sum of row i
	 */
	uli_t getSizeOfRow(unsigned int i) {
		uli_t size = 0;

		for (unsigned int j = 0; j < this->ncol; j++) {
			if (dep[i + j * nrow])
				++size;
		}

		return size;
	}

	/** 
	 * \brief Returns the sum of column j.
	 * 
	 * @param j Index of column.
	 * 
	 * @return Sum of column j.
	 */
	uli_t getSizeOfCol(unsigned int j) {
		uli_t size = 0;

		for (unsigned int i = 0; i < this->nrow; i++) {
			if (dep[i + j * nrow])
				++size;
		}

		return size;
	}

	/** 
	 * \brief Returns the row number of bit matrix.
	 * 
	 * @return Number of rows.
	 */
	inline uli_t getNrow() {
		return nrow;
	}

	/** 
	 * \brief Returns the column number of bit matrix.
	 * 
	 * @return Number of columns.
	 */
	inline uli_t getNcol() {
		return ncol;
	}

	/** 
	 * \brief Returns the total number of "true"s in bit matrix.
	 * 
	 * @return Total number of "true"
	 */
	inline uli_t count() {
		return dep.count();
	}

	/** 
	 * \brief Sends bit matrix to master process.
	 */
	virtual void sendMpi() {
#ifdef HAVE_MPI

		// slave
		std::string message;

		// send bitset
		MPI_Send(&nrow, 1, MPI_UNSIGNED, 0, MPI_TAG_BITMATRIX_NROW, MPI_COMM_WORLD);
		MPI_Send(&ncol, 1, MPI_UNSIGNED, 0, MPI_TAG_BITMATRIX_NCOL, MPI_COMM_WORLD);

		// convert bitset to string
		boost::to_string(dep, message);

		// send string and string's zero byte
		MPI_Send((char *)message.c_str(), message.size() + 1, MPI_CHAR, 0, MPI_TAG_BITMATRIX_DEP, MPI_COMM_WORLD);
#endif
	}

	/** 
	 * \brief Receives bit matrix from a slave. 
	 * 
	 * @param mpiId Slave id.
	 */
	virtual void receiveMpi(int mpiId) {
#ifdef HAVE_MPI
		MPI_Status status;

		// master
		char* message;

		// receive bitset
		MPI_Recv(&nrow, 1, MPI_UNSIGNED, mpiId, MPI_TAG_BITMATRIX_NROW, MPI_COMM_WORLD, &status);
		MPI_Recv(&ncol, 1, MPI_UNSIGNED, mpiId, MPI_TAG_BITMATRIX_NCOL, MPI_COMM_WORLD, &status);

		// string length comprises the zero byte
		message = new char[nrow * ncol + 1];

		// receive dep bitset
		MPI_Recv(message, nrow * ncol + 1, MPI_CHAR, mpiId, MPI_TAG_BITMATRIX_DEP, MPI_COMM_WORLD, &status);
		dep = boost::dynamic_bitset<>(std::string(message));

		delete[] message;
#endif
	}

	/** 
	 * \brief Receives bit matrix from a slave. 
	 * Matrix is coded as a string.
	 * 
	 * @param mpiId Slave id.
	 */
	virtual std::string receiveStrMpi(int mpiId) {
		std::string str;
#ifdef HAVE_MPI
		MPI_Status status;

		// master
		char* message;

		// receive bitset
		MPI_Recv(&nrow, 1, MPI_UNSIGNED, mpiId, MPI_TAG_BITMATRIX_NROW, MPI_COMM_WORLD, &status);
		MPI_Recv(&ncol, 1, MPI_UNSIGNED, mpiId, MPI_TAG_BITMATRIX_NCOL, MPI_COMM_WORLD, &status);

		// string length comprises the zero byte
		message = new char[nrow * ncol + 1];

		// receive dep bitset
		MPI_Recv(message, nrow * ncol + 1, MPI_CHAR, mpiId, MPI_TAG_BITMATRIX_DEP, MPI_COMM_WORLD, &status);
		str = std::string(message);

		delete[] message;

#endif
		return str;
	}

	/// bit matrix containing data-col dependencies
	boost::dynamic_bitset<> dep;

	/// Row number of bit matrix
	uli_t nrow;

	/// Column number of bit matrix
	uli_t ncol;
};

/** 
 * \brief Prints a bit matrix to stream.
 * 
 * @param os Output stream.
 * @param bitMatrix Bit matrix to be printed.
 * 
 * @return Output stream.
 */
inline std::ostream & operator<<(std::ostream &os, const BitMatrix &bitMatrix) {
	return bitMatrix.print(os);
}

#endif /*BITMATRIX_H_*/

