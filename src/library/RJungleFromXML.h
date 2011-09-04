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

#ifndef RJUNGLEFROMXML_H_
#define RJUNGLEFROMXML_H_

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <libxml/xmlreader.h>

#include "init.h"
#include "Exception.h"
#include "gzstream.h"
#include "RJunglePar.h"
#include "CmpldTree.h"
#include "treedefs.h"

template<class T>
class RJungleFromXML {
public:
  RJungleFromXML() {
  }
  ;
  RJungleFromXML(RJunglePar &par) : par(par) {
  }
  ;
  virtual ~RJungleFromXML() {
  }
  ;

  //#ifdef LIBXML_READER_ENABLED
  /**
   * Processes each node in XML file and creates a whole tree successively.
   * @param cmpldTree Pointer to location where the grown tree is stored
   */
  void processNode() {

    // declare parser variables
    const xmlChar *name, *value, *varName, *optionName;
    int type;

    // get type of node 1=opened and 15=closed
    type = xmlTextReaderNodeType(reader);

    // get name of node
    name = xmlTextReaderConstName(reader);


    // if name empty set it to --
    if (name == NULL) name = BAD_CAST "--";

    // get name of variable
    varName = xmlTextReaderGetAttribute(reader, (xmlChar *) "name");

    // if name empty set it to --
    if (varName == NULL) varName = BAD_CAST "--";

    // get name of attribute
    optionName = xmlTextReaderGetAttribute(reader, (xmlChar *) "id");

    // if name empty set it to --
    if (optionName == NULL) optionName = BAD_CAST "--";

    // node is closed
    if (type == 15) {
      if (strcmp((char *) name, (char *) "options") == 0) {

        // if preable is over then stop parsing for now
        stopParsing = true;

        return;

      } else if (strcmp((char *) name, (char *) "tree") == 0) {

        // create tree

        // create empty CmpldTree
        cmpldTree = new CmpldTree<T> (
            varSize, valueSize, valueWidth, branchSize);

        // assign cached scalars
        cmpldTree->varID = varID;
        cmpldTree->classes = classes;
        cmpldTree->indexes = indexes;
        cmpldTree->termID = termID;

        // assign cached vectors
        size_t i, j;

        cmpldTree->values = std::vector<std::vector<std::vector<T> > >(
            values.size(), std::vector<std::vector<T> >(
                valueWidth, std::vector<T>(valueSize, 0)));
        for (i = 0; i < values.size(); ++i) {
          for (j = 0; j < values[i].size(); ++j) {
            cmpldTree->values[i][j] = values[i][j];
          }
        }

        for (i = 0; i < extra.size(); ++i) {
					cmpldTree->extra.push_back(std::vector<double >());
          for (j = 0; j < extra[i].size(); ++j) {
            cmpldTree->extra.at(i).push_back(extra[i][j]);
          }
        }

        for (i = 0; i < branches.size(); ++i) {
          cmpldTree->branches.push_back(branches[i]);
        }

        return;
      }
    }

    // if a variable is declared then cache it
    if ((type == 1) && (strcmp((char *) name, (char *) "variable") == 0))
      curVarName = (char *) varName;
    else if ((type == 15)
        && (strcmp((char *) varName, curVarName.c_str()) == 0))
      curVarName = "";

    // if a option is declared then cache it
    if ((type == 1) && (strcmp((char *) name, (char *) "option") == 0))
      curOptionName = (char *) optionName;
    else if ((type == 15) && (strcmp((char *) varName, curOptionName.c_str())
        == 0))
      curOptionName = "";

    /*
     std::cout << curVarName << " " << std::endl;
     std::cout << curOptionName << " " << std::endl;


     printf("%d %d %s %s %d %d",
     xmlTextReaderDepth(reader),
     type,
     name,
     varName,
     xmlTextReaderIsEmptyElement(reader),
     xmlTextReaderHasValue(reader));
     */

    // get value
    value = xmlTextReaderConstValue(reader);

    // store content of scalars and vectors
    if (value != NULL) {
      if (type == 3) {
        std::string str((char *) value);
        std::istringstream istr(str);

        if (curOptionName != "") {
          if (curOptionName == "treetype")
            istr >> this->par.treeType;
          if (curOptionName == "ntree")
            istr >> this->par.ntree;
        }

        if (curVarName != "") {
          if (curVarName == "values")
            istr >> values;
          if (curVarName == "branches")
            istr >> branches;
          //if (curVarName == "indexes") istr >> indexes;
          if (curVarName == "varID")
            istr >> varID;
          //if (curVarName == "classes") istr >> classes;
          if (curVarName == "extra")
            istr >> extra;
          if (curVarName == "varSize")
            istr >> varSize;
          if (curVarName == "valueSize")
            istr >> valueSize;
          if (curVarName == "valueWidth")
            istr >> valueWidth;
          if (curVarName == "branchSize")
            istr >> branchSize;
          if (curVarName == "termID") {
            if (str == std::string("MAX"))
              termID = std::numeric_limits<uli_t>::max();
            else
              istr >> termID;
          }
        }
      }
    }

  }

  /**
   * Parse the RJungle XML file for next tree and store it.
   * @return A CmpldTree object
   */
  CmpldTree<T>* getNextTree() {

    cmpldTree = NULL; // tree which should be returned

    if (reader != NULL) { // reader was started

      // parse file and search for a tree node after node
      do {

        // Get next node
        ret = xmlTextReaderRead(reader);

        // process node and receive tree
        if (ret != 0) processNode();

      } while ((ret == 1) && (cmpldTree == NULL));

      if ((ret != 0) || (cmpldTree == NULL)) { // failed to parse
        //fprintf(stderr, "%s : failed to parse\n", data->par.predict);
      }
    } else { // no reader was initialized
      //fprintf(stderr, "Unable to open %s\n", data->par.predict);
    }

    return cmpldTree;
  }

  /**
   * Parse the RJungle XML file for parameters and store it.
   * @return RJungle Parameter
   */
  RJunglePar getPar() {

    if (reader != NULL) { // reader was started

      // reset parsing flag
      stopParsing = false;

      // parse file and search for a tree node after node
      do {

        // Get next node
        ret = xmlTextReaderRead(reader);

        // process node and receive par
        if (ret != 0) processNode();

      } while ((ret == 1) && (!stopParsing));

      if (ret != 0) { // failed to parse
        //fprintf(stderr, "%s : failed to parse\n", data->par.predict);
      }
    } else { // no reader was initialized
      //fprintf(stderr, "Unable to open %s\n", data->par.predict);
    }

    return par;
  }

  /**
   * Initialize XML reader and open file.
   * @param Data a pointer to the input data
   */
  void startXmlReader(DataFrame<T>* data) {

    // save data
    this->data = data;

    // initialize XML reader
    reader = xmlReaderForFile(data->par.predict, NULL, 0);

    if (reader == NULL)
      Exception(ERRORCODE_61);

    // test the version of libxml2
    LIBXML_TEST_VERSION;
  }

  /**
   * Stop XML reading and clean up the memory.
   */
  void stopXmlReader() {

    // free memory of XML reader
    xmlFreeTextReader(reader);

    //Cleanup function for the XML library.
    xmlCleanupParser();

    //this is to debug memory for regression tests
    xmlMemoryDump();

  }

	RJunglePar getParFromXmlFile() {

		// start XML reader
		startXmlReader(data);

		// fetch parameteres from XML file
		return getPar();
	}

  // xml reader variables
  xmlTextReaderPtr reader;
  int ret;

  // RJungle data
  DataFrame<T> *data;

  // RJungle parameter
  RJunglePar par;

  // Variable for stoping parsing
  bool stopParsing;

  // tree pointer for caching issues
  CmpldTree<T>* cmpldTree;

  // chached variables, will be the content of a compiled tree (cmpldTree)
  std::string curVarName;
  std::string curOptionName;

  std::vector<std::vector<uli_t> > varID;
  std::vector<std::vector<std::vector<T> > > values;
  std::vector<std::vector<uli_t> > branches;
  std::vector<T> classes;
  std::vector<uli_t> indexes;
	std::vector<std::vector<double> > extra;

  uli_t varSize;
  uli_t valueSize;
  uli_t valueWidth;
  uli_t branchSize;
  uli_t termID;

};

#endif /* RJUNGLEFROMXML_H_ */
