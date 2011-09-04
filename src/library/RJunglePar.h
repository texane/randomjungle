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

#ifndef RANDOMJUNGLEPAR_H_
#define RANDOMJUNGLEPAR_H_

#include <gsl/gsl_rng.h>

/*
 * Def.
 */
#ifndef HAVE__BOOL
#define bool int
#endif

//typedef unsigned long int uli_t;

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct RJunglePar {
    char *filename;
    char delimiter;
    unsigned int treeType;
    uli_t ntree;
    uli_t ntreeMpi; // Total number of trees when MPI is used
    uli_t mtry;
    unsigned int depVar;
    char *depVarName;
    uli_t nrow;
    uli_t ncol;
    unsigned int varNamesRow;
    unsigned int depVarCol;
    char *outprefix;
    uli_t skipRow;
    uli_t skipCol;
    int missingcode;
    unsigned int impMeasure;
    unsigned int backSel;
    uli_t numOfImpVar;
    bool downsampling_flag;
    bool verbose_flag;
    unsigned int memMode;
    unsigned int saveJungleType;
    char *predict;
    int varproximities;
    bool summary_flag;
    bool testlib_flag;
    char *colSelection;
    unsigned int imputeIt;
    bool gwa_flag;
    bool allcont_flag;
    bool transpose_flag;
    bool sampleproximities_flag;
    bool weightsim_flag;
    bool extractdata_flag;
    unsigned int seed;
    int nthreads;
    bool pedfile_flag;
    double tunemtry;
    int outlier;
    bool votes_flag;
    int prototypes;
    int mpi;
    int mpiId;
    int mpiSize;
    double condimp;
    bool permresponse_flag;

    char delimScale;
    double cutoffHighLD;

    gsl_rng *rng;
    char *version;

    // JTree tweaks
    uli_t maxTreeDepth;
    uli_t targetPartitionSize;
  } RJunglePar;

#ifdef __cplusplus
}
#endif

#endif /*RANDOMJUNGLEPAR_H_*/
