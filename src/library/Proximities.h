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

#ifndef PROXIMITIES_H_
#define PROXIMITIES_H_

/*
 * Includes
 */

#include <vector>
#include <limits>

#include "DataFrame.h"
#include "Exception.h"
#include "BitMatrix.h"
#include "treedefs.h"

/*
 * Proximities
 */

template<class T>
class Proximities {
public:

  Proximities() :
    nrow(0), ncol(0), proxMatrix(NULL), data(NULL) {
  }
  ;

  Proximities(uli_t nrow, uli_t ncol) :
    nrow(nrow), ncol(ncol), proxMatrix(NULL), data(NULL) {

    initMatrix();
  }
  ;

  Proximities(uli_t nrow, uli_t ncol, DataFrame<T> *data) :
    nrow(nrow), ncol(ncol), proxMatrix(NULL), data(data) {

    initMatrix();
  }
  ;

  Proximities(DataFrame<T> *data) :
    proxMatrix(NULL), data(data) {

    if (data == NULL)
      throw Exception(ERRORCODE_35);

    nrow = data->getnrow();
    ncol = nrow;

    initMatrix();
    reset();
  }
  ;

  void initWithData(DataFrame<T> *data) {
    this->data = data;
    if (data == NULL)
      throw Exception(ERRORCODE_35);

    nrow = data->getnrow();
    ncol = nrow;

    initMatrix();
    reset();
  }

  virtual ~Proximities() {
    uli_t i;
    if (this->proxMatrix != NULL) {
      for (i = 0; i < this->nrow; i++) {
        if (this->proxMatrix[i] != NULL)
          delete[] this->proxMatrix[i];
      }
      delete[] this->proxMatrix;
    }
  }

  // sample proximities

  void initMatrix() {
    // dynamic allocation of 2-D arrays
    uli_t i;
    try {
      if (this->proxMatrix != NULL) {
        for (i = 0; i < this->nrow; i++) {
          if (this->proxMatrix[i] != NULL)
            delete[] this->proxMatrix[i];
        }
        delete[] this->proxMatrix;
      }

      this->proxMatrix = new double*[this->nrow];
      for (i = 0; i < this->nrow; i++) {
        this->proxMatrix[i] = new double[this->ncol];
      }
    } catch (std::exception &e) {
      throw ;
    }
  }

  inline uli_t getnrow() {return nrow;};
  inline uli_t getncol() {return ncol;};

  std::vector<double> getRow(uli_t j) {
    std::vector<double> out;
    uli_t i;

    for (i = 0; i < ncol; ++i) {
      out.push_back(proxMatrix[j][i]);
    }

    return out;
  }

  std::vector<double> getCol(uli_t j) {
    std::vector<double> out;
    uli_t i;

    for (i = 0; i < nrow; ++i) {
      out.push_back(proxMatrix[i][j]);
    }

    return out;
  }

  std::ostream &print(std::ostream &os) const {
    uli_t i, j;

    os << "Proximity Matrix:" << std::endl;

    for(i = 0; i < this->nrow; i++) {
      for(j = 0; j < this->ncol; j++) {
        os << this->proxMatrix[i][j] << "\t";
      }
      os << std::endl;
    }

    //os << "Data Matrix:" << std::endl;

    //os << *data << std::endl;

    return os;
  }

  void printCSV(std::ostream &os) const {
    uli_t i, j;
    uli_t nrowReal = this->nrow;

    if (data->isUnsupervised) nrowReal /= 2;

    for(i = 0; i < nrowReal; i++) {
      if (i == data->par.depVar) continue;
      for(j = 0; j < nrowReal; j++) {
        if (j == data->par.depVar) continue;
        if (j> 0) os << data->par.delimiter;
        os << this->proxMatrix[i][j];
      }
      os << std::endl;
    }
  }

  void reset() {
    uli_t i, j;

    for(i = 0; i < this->nrow; i++)
    for(j = 0; j < this->ncol; j++)
    this->proxMatrix[i][j] = 0.0;
  }

  //fast, without security checks
  inline double at(uli_t i, uli_t j) const {return proxMatrix[i][j];}
  //inline const T* const operator[](uli_t i) {return proxMatrix[i];};

  inline void add(size_t i, size_t j, double val) {
    proxMatrix[i][j] += val;
  }

  void add(Proximities<T> &proxi) {
    uli_t i, j;

    for(i = 0; i < this->nrow; i++)
    for(j = 0; j < this->ncol; j++)
    this->proxMatrix[i][j] += proxi.proxMatrix[i][j];
  }

  void div(double val) {
    uli_t i, j;

    for(i = 0; i < this->nrow; i++)
    for(j = 0; j < this->ncol; j++)
    this->proxMatrix[i][j] /= val;
  }

  inline void div(size_t i, size_t j, double val) {
    proxMatrix[i][j] /= val;
  }

  inline void set(size_t i, size_t j, double val) {
    proxMatrix[i][j] = val;
  }

  void setDiag(double val) {
    uli_t i;

    for(i = 0; i < this->nrow; i++) {
      if (i < this->ncol) this->proxMatrix[i][i] = val;
    }
  }

  inline void setDim(size_t i, size_t j) {
    nrow = i;
    ncol = j;
  }

  /*
   * Combine proximities
   */
  virtual void combineMpi() {
#ifdef HAVE_MPI
    MPI_Status status;
    std::vector<double> dataVec;
		size_t nslaves, i;

    MPI_Barrier(MPI_COMM_WORLD);

		//std::cout << "Combining proximities..." << this->data->par.mpiId << std::endl;

		if (this->data->par.mpiId == 0) {
			// master

			std::vector<double> message(ncol, 0);

			// ntree vector
			message.resize(ncol);
			for (nslaves = 0; nslaves < (size_t)(this->data->par.mpiSize - 1); nslaves++) {
				// copy matrix row by row
				for (i = 0; i < nrow; i++) {
					MPI_Recv(&message[0], ncol, MPI_DOUBLE, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status);
					addVecs(proxMatrix[i], message);
				}

				// print status
				//std::cout << "Received proximities" << this->data->par.mpiId << std::endl;
			}

		} else {
			// slave

			for (i = 0; i < nrow; i++) {
				MPI_Send(proxMatrix[i], ncol, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
			}

			// print status
			//std::cout << "Sent proximities" << this->data->par.mpiId << std::endl;
		}

		// sync all processes
		MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

	void addVecs(double* &to, std::vector<double> &from) {
		for (unsigned int i = 0; i < from.size(); ++i)
			to[i] += from[i];
	}


  uli_t nrow;
  uli_t ncol;
  double **proxMatrix;
  DataFrame<T> *data;
};

template<class T>
inline std::ostream & operator<<(std::ostream &os, Proximities<T> &proximities) {
  return proximities.print(os);
}

#endif /*PROXIMITIES_H_*/

