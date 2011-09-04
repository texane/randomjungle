/*
 * SaveCollector.cpp
 *
 *  Created on: 31.07.2009
 *      Author: schwarz
 */

#include "SaveCollector.h"
#include "Helper.h"

SaveCollector::SaveCollector() {
  par = NULL;
  showHeader = false;
  showDepVar = true;
  repeatLast = false;
}

SaveCollector::~SaveCollector() {}

void SaveCollector::push_back(std::vector<double> * vec, std::string name) {
  orderCol.push_back(sc_double);
  colNames.push_back(name);
  doubleVec.push_back(vec);
}

void SaveCollector::push_back(std::vector<int> * vec, std::string name) {
  orderCol.push_back(sc_int);
  colNames.push_back(name);
  intVec.push_back(vec);
}

void SaveCollector::push_back(std::vector<unsigned int> * vec, std::string name) {
  orderCol.push_back(sc_size_t);
  colNames.push_back(name);
  size_tVec.push_back(vec);
}

void SaveCollector::push_back(std::vector<unsigned long int> * vec, std::string name) {
  orderCol.push_back(sc_uli_t);
  colNames.push_back(name);
  uli_tVec.push_back(vec);
}

void SaveCollector::push_back(std::vector<std::string> * vec, std::string name) {
  orderCol.push_back(sc_string);
  colNames.push_back(name);
  stringVec.push_back(vec);
}

void SaveCollector::clear() {
  orderCol.clear();
  doubleVec.clear();
  intVec.clear();
  size_tVec.clear();
  uli_tVec.clear();
  stringVec.clear();
}

std::ostream& SaveCollector::print(std::ostream& os) const {
  size_t i, j, row;
  char delimiter = (par == NULL) ? ' '
      : par->delimiter;

  // indexes
  std::vector<size_t> idx;
  idx.resize(sc_size, 0);

  // max index

  size_t idxMax, sizeMax;
  idxMax = sizeMax = 0;

  // set row order
  if (orderRow.size() > 0) {
    idxMax = orderRow.size();
  } else {
    for (j = 0; j < sc_size; ++j) {
      switch (j) {
      case sc_double:
        sizeMax = SaveCollector::getIdxMax(doubleVec);
        break;
      case sc_int:
        sizeMax = SaveCollector::getIdxMax(intVec);
        break;
      case sc_size_t:
        sizeMax = SaveCollector::getIdxMax(size_tVec);
        break;
      case sc_uli_t:
        sizeMax = SaveCollector::getIdxMax(uli_tVec);
        break;
      case sc_string:
        sizeMax = SaveCollector::getIdxMax(stringVec);
        break;
      }
      if (sizeMax > idxMax)
        idxMax = sizeMax;
    }
  }

  // print col names
  if (showHeader) {
    for (j = 0; j < colNames.size(); ++j) {
      if (j > 0)
        os << delimiter;

      os << colNames[j];
    }
    os << std::endl;
  }

  bool showRow;
  // print data
  for (j = 0; j < idxMax; ++j) {
    // reset
    idx.assign(sc_size, 0);

    row = (orderRow.size() > 0) ? orderRow[j]
        : j;

    // show row ?
    showRow = true;
    if (par != NULL) {
      if (!showDepVar && (par->depVar == row))
        showRow = false;

      if (row < isAvailable.size())
        if (isAvailable[row] == 0)
          showRow = false;
    }

    // show row !
    if (showRow) {
      for (i = 0; i < orderCol.size(); ++i) {
        if (i > 0)
          os << delimiter;

        switch (orderCol[i]) {
        case sc_double:
          SaveCollector::printElement(
              os, *doubleVec[idx[orderCol[i]]], row, repeatLast);
          break;
        case sc_int:
          SaveCollector::printElement(
              os, *intVec[idx[orderCol[i]]], row, repeatLast);
          break;
        case sc_size_t:
          SaveCollector::printElement(
              os, *size_tVec[idx[orderCol[i]]], row, repeatLast);
          break;
        case sc_uli_t:
          SaveCollector::printElement(
              os, *uli_tVec[idx[orderCol[i]]], row, repeatLast);
          break;
        case sc_string:
          SaveCollector::printElement(
              os, *stringVec[idx[orderCol[i]]], row, repeatLast);
          break;
        }

        ++idx[orderCol[i]];
      }
      os << std::endl;
    }
  }

  // hand on ostream
  return os;
}
