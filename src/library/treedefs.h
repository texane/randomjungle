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

#ifndef TREEDEFS
#define TREEDEFS

#define MPI_TAG_NTREE 1001
#define MPI_TAG_EXV 1003
#define MPI_TAG_EX2B 1004
#define MPI_TAG_EX2L 1005
#define MPI_TAG_BITMATRIX_NROW 1006
#define MPI_TAG_BITMATRIX_NCOL 1007
#define MPI_TAG_BITMATRIX_DEP 1008
#define MPI_TAG_DATATREESET_SLAVEID 1009
#define MPI_TAG_PREDICTION_LEN 1010
#define MPI_TAG_PREDICTION_PRED 1011
#define MPI_TAG_PREDICTION_ACC 1012


// you can use double, float, char, uli_t, int, ...
#define MYDATATYPE float

#ifdef __SPARSE_DATA__
#define MISSINGCODE 3
#else
#define MISSINGCODE -99
#endif

#define SAVECOLLECTORTYPE
enum SaveCollectorType {
  sc_double = 0,
  sc_int = 1,
  sc_size_t = 2,
  sc_uli_t = 3,
  sc_string = 4,
  sc_size = 5
};

#define LVLOFMEASUREMENT
enum LvlOfMeasurement {
  lom_undef = 0,
  lom_nominal = 1,
  lom_ordinal = 2,
  lom_numeric = 3
};

#define CLASSATOMTYPE
enum ClassAtomType {
  ca_BASE = 0,
  ca_TERM = 3,
  ca_T2 = 4,
  ca_TM = 5,
  ca_S = 6,
  ca_IAM = 7
};

#define TREETYPE
enum TreeType {
  // CART
	tt_CART = 1,	        // Class. trees: y nominal, x numeric
	tt_CARTcateCate = 2,	  // Class. trees: y nominal, x nominal
	tt_CARTcontCont = 3,       // Reg.   trees: y numeric, x numeric
  tt_CARTcontCate = 4, // Reg.   trees: y numeric, x nominal

  // CART extension
	tt_CARTsrt = 5,       // Class. trees: y nominal, x numeric (faster for many different values in x)

  // Experimental
  tt_LOTUS = 205, // LOTUS (Tree presentation with CART)
};

#define VARTYPE
enum VarType {
  vt_CONT = 1, // contineous
  vt_CATE = 2  // categorical
};


#define IMPORTANCEMEASURE
enum ImportanceMeasure {
	im_intrinsic = 1,	// Intrinsic imp. (i.e. GINI-Index)
	im_perm_breiman = 2, // Permutation imp. Breiman
  im_perm_liaw = 3, // Permutation imp. Liaw
  im_perm_raw = 4, // Permutation imp. Raw / no normalization
  im_perm_meng = 5, // Permutation imp. Meng
};

#define BACKWARDELIMINATION
enum BackwardElimination {
  bs_50P = 1,// backward elimination. grow RJ with all vars, calc GINI-Index imp, take top50% vars,
                      // then grow tree with them and so on. until numOfImpVar vars remain.
                      // Discard 50% at each step
  bs_33P = 2,// backward elimination. grow RJ with all vars, calc GINI-Index imp, take top50% vars,
                      // then grow tree with them and so on. until numOfImpVar vars remain.
                      // Discard 25% at each step
  bs_NEG = 3,
  bs_DIAZURIATE = 4
};


#define VARPROXI
enum VarProxi {
  vp_SIMPLE = 1,
  vp_EXTEND1 = 2
};


#define GSL_RNG_SEED 15051980
#define UZI unsigned short int
//typedef unsigned long int uli_t;
//typedef size_t uli_t; // == size_t
typedef unsigned int uli_t;

#ifndef NULL
#define NULL 0
#endif

#define TERMINALNODE 3

// detect OS
#if !defined( __WINDOWS__ ) && ( defined( _Windows ) || defined( _WINDOWS ) )
  #define __WINDOWS__
#endif /* !__WINDOWS__ && ( _Windows || _WINDOWS ) */
#if !defined( __WIN32__ ) && ( defined( WIN32 ) || defined( _WIN32 ) )
  #ifndef __WINDOWS__
    #define __WINDOWS__
  #endif /* __WINDOWS__ */
  #define __WIN32__
#endif /* !__WIN32__ && ( WIN32 || _WIN32 ) */
#if defined( __WINDOWS__ ) && !defined( __WIN32__ )
  #define __WIN16__
#endif /* __WINDOWS__ && !__WIN32__ */


// library for windows
#ifdef __WINDOWS__
  // calling convention, should be used for functions called by dll-client
	#ifdef __BUILDING_DLL__
		#define __lrj __declspec(dllexport)
	#else
		#define __lrj __declspec(dllimport)
	#endif
#endif

// library for ...
#ifndef __WINDOWS__
	#define __lrj
	#define __cdecl
#endif

#endif //TREEDEFS

