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

#ifndef ERRORCODES_H_
#define ERRORCODES_H_

#define ERRORCODE_1 "Invalid index number (perhaps to big)"
#define ERRORCODE_2 "Writing a duplication of a pointer to a vector (same adresses)"
#define ERRORCODE_3 "Writing a wrong sized variable to a dataframe"
#define ERRORCODE_4 "Writing a variable to a dataframe, which name is already in use (same name)"
#define ERRORCODE_5 "No data served (NULL-Pointer)"
#define ERRORCODE_6 "Outcome vector has not got one or two classes."
#define ERRORCODE_7 "File not found."
#define ERRORCODE_8 "Error while parsing csv file."
#define ERRORCODE_9 "No fitting function giving in constructor of tree"
#define ERRORCODE_10 "No tree served (NULL-Pointer)"
#define ERRORCODE_11 "Empty tree served (root node == NULL-Pointer)"
#define ERRORCODE_12 "Vector has not got two classes."
#define ERRORCODE_13 "Dependent variable has to be binary"
#define ERRORCODE_14 "Classifier points to a children node that does not exist."
#define ERRORCODE_15 "Classifier points to a children node that is a NULL pointer."
#define ERRORCODE_16 "Infile has less columns than wanted."
#define ERRORCODE_17 "Unknown level of measurement."
#define ERRORCODE_18 "Dependency variable index is bigger than the number of columns."
#define ERRORCODE_19 "Can not find specified dependency variable in the data. Index too big."
#define ERRORCODE_20 "Can not find classes 0,1 in the dependency variable."
#define ERRORCODE_21 "Can not open file for writing."
#define ERRORCODE_22 "Has got allocated memory already. Array pointer not NULL."
#define ERRORCODE_23 "Array is empty. Array pointer is NULL."
#define ERRORCODE_24 "Error while calculating with mpfr."
#define ERRORCODE_25 "Can not find file. Or number of rows and/or columns is equals to 0. Or more column names than columns in data."
#define ERRORCODE_26 "Child is NULL"
#define ERRORCODE_27 "No Children there"
#define ERRORCODE_28 "In makeNode recursion depth >= 15000"
#define ERRORCODE_31 "Do not know this person"
#define ERRORCODE_32 "Tree size is zero"
#define ERRORCODE_33 "The mtry < 2"
#define ERRORCODE_34 "Unknown dependent variable name"
#define ERRORCODE_35 "DataFrame is NULL"
#define ERRORCODE_36 "Too less columns for output. (<2)"
#define ERRORCODE_37 "Missings in Data"
#define ERRORCODE_38 "The mtry is bigger than nmber of columns (ncol < mtry)"
#define ERRORCODE_39 "Unknown memmode"
#define ERRORCODE_41 "Can not find column specified in column selection file (-C parameter)"
#define ERRORCODE_42 "Can not find category"
#define ERRORCODE_43 "Variables in data are wrongly named or orderes (FID, IID, PAT, MAT, SEX, PHENOTYPE)"
#define ERRORCODE_44 "Is not a valid ped format or can not use options together: --pedfile and --transpose"
#define ERRORCODE_45 "Option -D has to been PHENOTYPE or leave it empty"
#define ERRORCODE_46 "The delimeter has to be a SPACE character in --pedfile mode"
#define ERRORCODE_47 "Can not find variable with that name"
#define ERRORCODE_48 "Data has not got 2 groups to classify on"
#define ERRORCODE_49 "Minimal partition size is too small (choose option -l 20 or higher)"
#define ERRORCODE_50 "This type of tree has got no function for updating proximities"
#define ERRORCODE_51 "Do not use backward elimination and tune mtry together"
#define ERRORCODE_52 "Do not use backward elimination and sample proximities together"
#define ERRORCODE_53 "Do not use backward elimination and variable proximities together"
#define ERRORCODE_54 "A data cell in the input file is empty"
#define ERRORCODE_55 "Invalid PED format. Needed: FID IID PAT MAT SEX PHENOTYPE"
#define ERRORCODE_56 "Can not fetch next tree from file. RJungle XML file corrupt? Or can not read from XML file. (Not there?)"
#define ERRORCODE_57 "Some rows differs in number of columns (or a space at the end of a column?)"
#define ERRORCODE_58 "Perform permutation importance when using MPI"
#define ERRORCODE_59 "Number of processes is greater than ntree. That is not allowed."
#define ERRORCODE_60 "Unknown format (--write)"
#define ERRORCODE_61 "Can not read from XML file"
#define ERRORCODE_62 "Conditional importrance is available in the following tree types (more will be added in future): 1, 3, 5"
#define ERRORCODE_63 "CPU does not support SSE"
#define ERRORCODE_64 "Lotus can be used in <<double or float>>-memory mode exclusively."

#endif /* ERRORCODES_H_ */
