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

#ifndef SAVECOLLECTOR_H_
#define SAVECOLLECTOR_H_

#include <iostream>
#include <vector>
#include <string>

#include "treedefs.h"
#include "RJungleIO.h"
#include "RJunglePar.h"

class SaveCollector {
public:
  SaveCollector();
  virtual ~SaveCollector();

  void push_back(std::vector<double> *, std::string);
  void push_back(std::vector<int> *, std::string);
  void push_back(std::vector<unsigned int> *, std::string);
  void push_back(std::vector<unsigned long int> *, std::string);
  void push_back(std::vector<std::string> *, std::string);

  void clear();

  std::ostream& print(std::ostream&) const;

  /*
   * Clac. Idx Max
   */
  template<class T>
  static size_t getIdxMax(const std::vector<std::vector<T> *> &dataVec) {
    size_t idxMax, size;
    idxMax = size = 0;

    for (size_t i = 0; i < dataVec.size(); ++i) {
      size = dataVec[i]->size();
      if (size > idxMax)
        idxMax = size;
    }

    return idxMax;
  }

  // print it if there
  template<class T>
  static void printElement(
    std::ostream& os, const std::vector<T> &vec, size_t pos, bool repeatLast) {
    if (pos < vec.size()) {
      os << vec[pos];
    } else {
      if (repeatLast && (vec.size() > 0)) {
        os << vec[vec.size() - 1];
      } else {
        os << "NA";
      }
    }
  }

  bool showHeader; // should the header be shown
  bool showDepVar; // should dependent variable be shown
  bool repeatLast; // repeat last element in vector if vector is shorter than
  // maximal length

  std::vector<size_t> orderCol;
  std::vector<size_t> orderRow;
  std::vector<std::string> colNames;
  std::vector<char> isAvailable;

  std::vector<std::vector<double> *> doubleVec;
  std::vector<std::vector<int> *> intVec;
  std::vector<std::vector<unsigned int> *> size_tVec;
  std::vector<std::vector<unsigned long int> *> uli_tVec;
  std::vector<std::vector<std::string> *> stringVec;

  RJunglePar *par;
};

// ostream
inline std::ostream& operator<<(
  std::ostream& os, const SaveCollector &saveCollector) {
  return saveCollector.print(os);
}

#endif /* SAVECOLLECTOR_H_ */
