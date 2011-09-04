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

#ifndef RJUNGLEHELPER_H_
#define RJUNGLEHELPER_H_

/*
 * Includes
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>

#include "RJunglePar.h"
#include "RJungleIO.h"
#include "RJungleGen.h"
#include "DataFrame.h"
#include "Helper.h"
#include "treedefs.h"

template<class T>
class RJungleHelper {
public:
  RJungleHelper();
  virtual ~RJungleHelper();

  static void printXml(RJungleIO &io, std::vector<Tree<T, uli_t> *> &trees) {

    typename std::vector<Tree<T, uli_t> *>::iterator itTrees;
    uli_t i;

    *io.outVerbose
        << "Warning: Not implemented yet. (Writing to jungle type ID = 1)"
        << std::endl;

    *io.outXmlJungle << "<jungle type=\"exhaustive\" id=\"0\" size=\""
        << trees.size() << "\">" << std::endl;
    i = 0;
    /*
     itTrees = trees.begin();
     while (itTrees != trees.end()) {
     *io.outXmlJungle << "<tree id=\"" << i << "\">" << std::endl;
     (*itTrees)->printXml(*io.outXmlJungle);
     *io.outXmlJungle << "</tree>" << std::endl;

     ++i;
     ++itTrees;
     }
     */
    *io.outXmlJungle << "</jungle>" << std::endl;

  }

  static void printHeader(RJunglePar &par, RJungleIO &io, time_t &start) {

    *(io.outVerbose) << "Start: " << ctime(&start)
        << "+---------------------+-----------------+-------------------+"
        << std::endl 
				<< "|    Random Jungle    |" << std::setw(12) << par.version
				<< std::setw(5) << " "
        << "|        2010       |" << std::endl
        << "+---------------------+-----------------+-------------------+"
        << std::endl
        << "|    (C) 2008-2010 Daniel F Schwarz et al., GNU GPL, v3     |"
        << std::endl
        << "|     https://sourceforge.net/projects/rjungle/support      |"
        << std::endl
        << "+-----------------------------------------------------------+"
        << std::endl << std::endl << "Output to: " << par.outprefix << ".*"
        << std::endl << "Loading data... " << std::endl;

  }

  static void printFooter(
      RJunglePar &par, RJungleIO &io, time_t &start, time_t &end) {
    *io.outVerbose << "Elapsed time: " << difftime(end, start) << " sec" << std::endl;
	  *io.outVerbose << "Finished: " << ctime(&end) << std::endl;
  }

  static void printRJunglePar(RJunglePar &par, std::ostream &os) {
    os << "file: " << par.filename << std::endl << "delimiter: "
        << par.delimiter << std::endl << "treetype: " << par.treeType
        << std::endl << "ntree: " << par.ntree << std::endl << "mtry: "
        << par.mtry << std::endl << "depvar: " << par.depVar << std::endl
        << "depvarname: " << par.depVarName << std::endl << "nrow: "
        << par.nrow << std::endl << "ncol: " << par.ncol << std::endl
        << "varnamesrow: " << par.varNamesRow << std::endl << "depvarcol: "
        << par.depVarCol << std::endl << "outprefix: " << par.outprefix
        << std::endl << "skiprow: " << par.skipRow << std::endl << "skipcol: "
        << par.skipCol << std::endl << "missingcode: " << par.missingcode
        << std::endl << "impmeasure: " << par.impMeasure << std::endl
        << "backsel: " << par.backSel << std::endl << "nimpvar: "
        << par.numOfImpVar << std::endl << "downsampling: "
        << par.downsampling_flag << std::endl << "verbose: "
        << par.verbose_flag << std::endl << "memmode: " << par.memMode
        << std::endl << "write: " << par.saveJungleType << std::endl
        << "predict: " << par.predict << std::endl << "varproximities: "
        << par.varproximities << std::endl << "summary: " << par.summary_flag
        << std::endl << "testlib: " << par.testlib_flag << std::endl
        << std::endl << "colselection: "
        << par.colSelection << std::endl << "impute: " << par.imputeIt
        << std::endl << "gwa: " << par.gwa_flag << std::endl << "impcont: "
        << par.allcont_flag << std::endl << "transpose: " << par.transpose_flag
        << std::endl << "sampleproximities: " << par.sampleproximities_flag
        << std::endl << "weightsim: " << par.weightsim_flag << std::endl
        << "extractdata: " << par.extractdata_flag << std::endl
        << std::endl << "seeed: " << par.seed << std::endl
        << "nthreads: " << par.nthreads << std::endl << "pedfile: "
        << par.pedfile_flag << std::endl << "maxTreeDepth: "
        << par.maxTreeDepth << std::endl << "targetPartitionSize: "
        << par.targetPartitionSize << std::endl << std::endl << std::endl;
  }

  // refresh xml entries
  static void printXmlHeader(RJunglePar &par, std::ostream &os) {
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl
        << "<!-- created by randomJungle " << par.version << " "
        << "under GPL v2 or later. (c) Daniel F. Schwarz 2008-2009."
        << "All rights reserved. -->" << std::endl << "<randomjungle>"
        << std::endl << "<options>" << std::endl << "<option id=\"file\">"
        << par.filename << "</option>" << std::endl
        << "<option id=\"delimiter\">" << par.delimiter << "</option>"
        << std::endl << "<option id=\"treetype\">" << par.treeType
        << "</option>" << std::endl << "<option id=\"ntree\">" << par.ntree
        << "</option>" << std::endl << "<option id=\"mtry\">" << par.mtry
        << "</option>" << std::endl << "<option id=\"depvar\">" << par.depVar
        << "</option>" << std::endl << "<option id=\"depvarname\">"
        << par.depVarName << "</option>" << std::endl << "<option id=\"nrow\">"
        << par.nrow << "</option>" << std::endl << "<option id=\"ncol\">"
        << par.ncol << "</option>" << std::endl
        << "<option id=\"varnamesrow\">" << par.varNamesRow << "</option>"
        << std::endl << "<option id=\"depvarcol\">" << par.depVarCol
        << "</option>" << std::endl << "<option id=\"outprefix\">"
        << par.outprefix << "</option>" << std::endl
        << "<option id=\"skiprow\">" << par.skipRow << "</option>" << std::endl
        << "<option id=\"skipcol\">" << par.skipCol << "</option>" << std::endl
        << "<option id=\"missingcode\">" << (int) par.missingcode
        << "</option>" << std::endl << "<option id=\"impmeasure\">"
        << par.impMeasure << "</option>" << std::endl
        << "<option id=\"backsel\">" << par.backSel << "</option>" << std::endl
        << "<option id=\"nimpvar\">" << par.numOfImpVar << "</option>"
        << std::endl << "<option id=\"downsampling\">" << par.downsampling_flag
        << "</option>" << std::endl << "<option id=\"verbose\">"
        << par.verbose_flag << "</option>" << std::endl
        << "<option id=\"memMode\">" << par.memMode << "</option>" << std::endl
        << "<option id=\"write\">" << par.saveJungleType << "</option>"
        << std::endl << "<option id=\"predict\">" << par.predict << "</option>"
        << std::endl << "<option id=\"varproximities\">" << par.varproximities
        << "</option>" << std::endl << "<option id=\"summary\">"
        << par.summary_flag << "</option>" << std::endl
        << "<option id=\"testlib\">" << par.testlib_flag << "</option>"
        << std::endl << "<option id=\"colselection\">" << par.colSelection
        << "</option>" << std::endl << "<option id=\"impute\">" << par.imputeIt
        << "</option>" << std::endl << "<option id=\"gwa\">" << par.gwa_flag
        << "</option>" << std::endl << "<option id=\"impcont\">"
        << par.allcont_flag << "</option>" << std::endl
        << "<option id=\"transpose\">" << par.transpose_flag << "</option>"
        << std::endl << "<option id=\"sampleproximities\">"
        << par.sampleproximities_flag << "</option>" << std::endl
        << "<option id=\"weightsim\">" << par.weightsim_flag << "</option>"
        << std::endl << "<option id=\"extractdata\">" << par.extractdata_flag
        << "</option>" << std::endl << "<option id=\"yaimp\">"
        << par.extractdata_flag << "</option>" << std::endl
        << "<option id=\"seeed\">" << par.seed << "</option>" << std::endl
        << "<option id=\"nthreads\">" << par.nthreads << "</option>"
        << std::endl << "<option id=\"maxtreedepth\">"
        << par.maxTreeDepth << "</option>" << std::endl
        << "<option id=\"targetpartitionsize\">" << par.targetPartitionSize
        << "</option>" << std::endl << "<option id=\"pedfile\">"
        << par.pedfile_flag << "</option>" << std::endl << "</options>"
        << std::endl << "<jungle type=\"raw\" id=\"0\" size=\"" << par.ntree
        << "\">" << std::endl;
  }

  static void printXmlFooter(RJungleIO &io) {
    *io.outXmlJungle << "</jungle>" << std::endl;
    *io.outXmlJungle << "</randomjungle>" << std::endl;
  }

  static void printXmlRaw(RJungleIO &io, CmpldTree<T>* &cmpldTree, uli_t treeId) {
    *io.outXmlJungle << "<tree id=\"" << treeId << "\">" << std::endl;
    cmpldTree->printXml(*io.outXmlJungle);
    *io.outXmlJungle << "</tree>" << std::endl;

  }

  static void summary(RJungleIO &io, CmpldTree<T> &cmpldTree) {

    //typename std::vector<Tree<T, uli_t> *>::iterator it;
    //double nmLeafs = 0.0;
    //double nmNodes = 0.0;

		// number of nodes 
    *io.outSummary << cmpldTree.branches.size() << std::endl;

		/*
    *io.outSummary << "Number of leafs:" << std::endl;
    it = trees.begin();
    while (it != trees.end()) {
      nmLeafs += (*it)->root->getNumOfLeafs();
      *io.outSummary << (*it)->root->getNumOfLeafs() << std::endl;
      ++it;
    }

    *io.outSummary << "Sum of leafs:" << std::endl;
    it = trees.begin();
    while (it != trees.end()) {
      *io.outSummary << (*it)->root->getLeafSum() << std::endl;
      ++it;
    }

    *io.outSummary << "Number of nodes:" << std::endl;
    it = trees.begin();
    while (it != trees.end()) {
      nmNodes += (*it)->root->getNumOfNodes();
      *io.outSummary << (*it)->root->getNumOfNodes() << std::endl;
      ++it;
    }
		*/

    //*io.outSummary << "Mean number of leafs:" << nmLeafs/trees.size() << std::endl;
    //*io.outSummary << "Mean number of nodes:" << nmNodes/trees.size() << std::endl;
  }

  static void tuneMtry(RJunglePar &par, std::vector<uli_t> *colMaskVec) {
    if (par.mtry == 0) {
      par.mtry = (uli_t) sqrt((double) ((colMaskVec) ? (colMaskVec->size())
          : (par.ncol)));
    }
  }

  static void tuneNtree(RJunglePar &par, std::vector<uli_t> *colMaskVec) {
    double probNtree = 0.999;

    // get new tree size
    if (par.ntree == 0) {
      par.ntree = Helper::getTreeSizeProb(probNtree, par.mtry, (colMaskVec
          == NULL) ? (par.ncol)
          : (colMaskVec->size()));
    }
  }

  static uli_t setMtryClassification(uli_t ncol) {
    return (uli_t)floor(sqrt((double)ncol));
  }

  static uli_t setMtryRegression(uli_t ncol) {
    return (uli_t)floor(((double)ncol) / 3.0);
  }

};

#endif /* RJUNGLEHELPER_H_ */
