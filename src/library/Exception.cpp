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

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "Exception.h"

Exception::Exception(std::string str): str(str), status(0) {
  std::cerr << "ERROR:" << std::endl << str << std::endl;
}

const char* Exception::what() const throw() {
  std::string strret("");

  if (status != 0) {
    char result[100];
    sprintf(result, "sta:%d ", status);
    strret.append(result);
  }

  return strret.c_str();
}
