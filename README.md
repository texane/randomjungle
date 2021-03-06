THE BASIC Random Jungle (RJ) README
===================================
See "INSTALL" for help on installation


INTRODUCTION
---------------

This directory contains the source code, tests, autobuild files and autobuild tools for RJ. RJ is a generalized implementation of Random Forests(tm) (RF) by Leo Breiman and Adele Cutler. In genetics, it is an well established method for analysing genetic data (i.e. GWA data).

RJ is free software distributed under a GNU-style copyleft.

RF is a powerful machine learning method. [Most interesting features](http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm):

* Variable selection: 	`Estimates the importance of variables`
* Efficiency:			`It runs efficiently on large data`
* Missing values: 		`It has an effective method for estimating missing data`
* Classifier:			`It creates classifier for further analyses`
* General. err.:		`It generates an internal unbiased estimate of the generalization error as the building progresses`
* Proximities:			`It computes proximities between pairs of cases that can be used in clustering, locating outliers, or (by scaling) give interesting views of the data. Can be used as the distance matrix for Multidimensional Scaling (MDS)`

HISTORY
----------
RJ was initially written by Daniel F. Schwarz of the Institute of Medical Biometry and Statistics of University of Lubeck in Germany. He started in the beginning of 2008.

GOALS
--------
Create a machine learning tool for analysing data as an equivalent to alternative statistical methods.

DEPS
-------
UBUNTU (packages):

* gcc
* g++
* libboost-dev
* zlib1g-dev
* libgsl0-dev
* libxml2-dev

MPI:
----

* mpi-default-dev
* mpi-default-bin

Linker libraries:

* xml2
* gsl
* gslcblas
* z 


Make documentation:

* texinfo
* texlive


Sincerely,
Daniel


[Random Forests Website](http://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm)
