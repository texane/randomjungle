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

#ifndef JTREETWEAKS_H_
#define JTREETWEAKS_H_

typedef unsigned long int uli_t;

#ifdef __cplusplus
extern "C" {
#endif

  // UNUSED STRUCT ! SEE RJunglePar.h
  typedef struct JTreeTweaks {
    uli_t maxTreeDepth;
    uli_t targetPartitionSize;
  } JTreeTweaks;

#ifdef __cplusplus
}
#endif

#endif /* JTREETWEAKS_H_ */
