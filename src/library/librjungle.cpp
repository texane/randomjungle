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

/*
 * Includes
 */

#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "config.h"
#include "treedefs.h"
#include "librjungle.h"
#include "Exception.h"
#include "gzstream.h"
#include "CmplFct.h"
#include "RJungleIO.h"
#include "RJungleGen.h"
#include "RJungleCtrl.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

/*
 * Source
 */

//void __cdecl randomJungle(RJunglePar par){
void randomJungle(RJunglePar par) {
#ifdef HAVE_MPI
  // test for permutation importance
  if (par.impMeasure < 2) {
    throw Exception(ERRORCODE_58);
  }

  // variables only for mpi
  std::stringstream sstr;
  std::string str(par.outprefix);

  // init mpi
  MPI_Comm_size(MPI_COMM_WORLD, &par.mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &par.mpiId);

  // test for ntree vs mpi size
  if (par.mpiSize - 1 > par.ntree) {
    throw Exception(ERRORCODE_59);
  }
#endif

  //set number of threads that will be used
#ifdef _OPENMP
#ifndef HAVE_MPI
  par.mpiId = 0;
#endif
  if ((par.nthreads < omp_get_max_threads()) && (par.nthreads != 0)) {
    // omp_set_dynamic defined different ???
#ifdef __WINDOWS__
    //omp_set_dynamic(0);
#else
    //omp_set_dynamic(1);
#endif
    omp_set_num_threads(par.nthreads);
  }
  omp_set_num_threads(par.nthreads);
#endif

#ifdef HAVE_MPI
  // reset parameters of different open mpi processes
  // rename output file
  if (par.mpiId != 0) {
    sstr << par.mpiId;
    str.append("_mpi_id_");
    str.append(sstr.str());
    par.outprefix = (char *) str.c_str();
  }
  // modify seed

  // modify number of tree(s)
  uli_t ntreeBak = par.ntree;
  par.ntreeMpi = par.ntree;
  par.ntree = (uli_t) floor(par.ntree / par.mpiSize);

  // allot remaining tree(s)
  int remainder = ntreeBak % par.mpiSize;
  if (par.mpiId < remainder)
    par.ntree++;
#endif

  // choose function regarding memory mode
#ifdef __SPARSE_DATA__
  par.memMode = 2;
	RJungleCtrl<char >::autoBuild(par);
#else
  switch (par.memMode) {
    case 0:
      RJungleCtrl<double>::autoBuild(par);
      break;
    case 1:
      RJungleCtrl<float>::autoBuild(par);
      break;
    case 2:
      RJungleCtrl<char>::autoBuild(par);
      break;
    default:
      throw Exception(ERRORCODE_39);
  }
#endif

}

RJunglePar __cdecl initRJunglePar() {
  RJunglePar par;
#ifdef HAVE_LIBZ
  par.filename = (char *) "input.dat.gz";
  par.predict = (char *) ""; //par.predict = (char *)"rjungle.jungle.xml.gz";
#else
  par.filename = (char *)"input.dat";
  par.predict = (char *)""; //par.predict = (char *)"rjungle.jungle.xml";
#endif
  par.depVarName = (char *) "";
  par.colSelection = (char *) "";
  par.outprefix = (char *) "rjungle";
  par.delimiter = ' ';
  par.delimScale = '@';
  par.treeType = 1; // type of base classifier
  par.ntree = 500;
  par.mtry = 0; // var sert size at each node while growing one tree [if = 0 then mtry = sqrt(samples)]
  par.depVar = 0; // which var is the dep one in the data
  par.nrow = 0; // number of rows (samples) in data [if = 0 then read all rows from data file]
  par.ncol = 0; // number of columns (variables) in data [if = 0 then read all columns  from data file]
  par.varNamesRow = 0; // at what position are the var names in the file
  par.depVarCol = 0; // at what position is the dep var in the data file
  par.skipRow = 0; // how many rows should be skipped while reading the file
  par.skipCol = 0; // how many columns should be skipped
  par.missingcode = MISSINGCODE;
  par.impMeasure = 1; // imp. measure
  par.backSel = 0; // back. sel.
  par.numOfImpVar = 100;
  par.imputeIt = 0;
  par.verbose_flag = false;
  par.downsampling_flag = false;
#ifdef __SPARSE_DATA__
  par.memMode = 3; //sparse
#else
  par.memMode = 0;
#endif
  par.saveJungleType = 0;
  par.varproximities = 0;
  par.sampleproximities_flag = false;
  par.summary_flag = false;
  par.testlib_flag = false;
  par.weightsim_flag = false;
  par.extractdata_flag = false;
  par.gwa_flag = false;
  par.allcont_flag = false;
  par.transpose_flag = false;
  par.pedfile_flag = false;
  par.seed = 1;
  par.nthreads = 0;
  par.rng = NULL;
  par.cutoffHighLD = 0.8;
  par.version = (char *) PACKAGE_VERSION;
  par.tunemtry = 0;
  par.outlier = 0;
  par.votes_flag = false;
  par.prototypes = 0;
  par.mpi = 0;
  par.mpiId = -1;
  par.mpiSize = 0;
  par.condimp = -1.0;
  par.permresponse_flag = 0;

  // JTree tweaks
  par.maxTreeDepth = 15000;
  par.targetPartitionSize = 1;

  return par;
}
;
