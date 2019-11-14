/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                    Definition of class FFI for Free Format Input

  ==============================================================================*/

#ifndef __FFI_H
#define __FFI_H

#include <stdlib.h>
#include <ctype.h>

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::ostream;
using std::endl;

#include <fstream>
#include <iomanip>
using std::setw;
using std::ios;
using std::setprecision;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "OFELI_Config.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


/*! \file FFI.h
 *  \brief Definition file for class FFI.
 */

/*! \class FFI
 *  \ingroup IO
 *  \brief To handle input files using free format format.
 *
 * FFI is a format (or protocol) for data files. This protocol is intended
 * to simplify the comprehension of file contents. For example, mesh and field files use the
 * FFI standard:
 *
 * \li Any string following a \c \# or \c \% sign is interpreted as a comment.
 * \li A string between a \c \# and a \c ! in the first line can be interpreted as a
 * ``magic number", \e i.e., a string used to inform about the nature of the
 * file contents. \e e.g., the string \c \#MESH! informs that the current file
 * is a Mesh data (<a href="../fileformats/MDFFormat.htm">MDF</a>) file.
 * \li Data can be given separated by any number of blanks or a carriage return.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define LENGTH_CMD   30
#define LINE_LENGTH 120
#define NB_CMD       12
#define FFI_ASCII     0
#define FFI_BINARY    1

int _is_number(const string &);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

class FFI
{

 public:

/// \brief Constructor using standard input
    FFI();

/// \brief Constructor using file name
/// @param [in] file Input file name
    FFI(const string &file);

/** \brief Constructor using file name and seeking a "magic number"
 *  @param [in] file Input file name
 *  @param [in] ident Identification string
 */
    FFI(const string &file, const string &ident);

/** \brief Constructor using a list to keywords to identify
 *  @param [in] s Vector containing keywords
 */
    FFI(const vector<string>& s);

/// Copy constructor
    FFI(const FFI &ff);

/// Destructor
    ~FFI();

/// Open file.
/// @param file [in] Input file name
    void open(const string &file);

/** Open file and seek a "magic number" \a ident.
 *  @param [in] file Input file name
 *  @param [in] ident Identification string
 */
    void open(const string &file, const string &ident);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void open(const string &file, const string &ident, int c);
    void Open(const string &file) { open(file); }
    void Open(const string &file, const string &ident) { open(file,ident); }
    void Open(const string &file, const string &ident, int c) { open(file,ident,c); } 
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// Close file
    void Close(void) { _is->close(); }

/// Read and return an integer number
/// @param [in] msg Message to output (for interactive usage)
    int getI(const string &msg="0");

/// Read and return a float number
/// @param [in] msg Message to output (for interactive usage)
    real_t getD(const string &msg="0");

/// Set a list of keywords to scan later
/// @param [in] s Array of keywords to initialize
    void setKeywords(const vector<string>& s);

/// Read and return a keyword in the list kw
/// @param [in] msg Message to output (for interactive usage)
    int getKW(const string &msg="0");

/// Read and return a string
/// @param [in] msg Message to output (for interactive usage)
    string getS(const string &msg="0");

/// Read and return an algebraic expression
/// @param [in] msg Message to output (for interactive usage)
    string getE(const string &msg="0");

/// Operator =
/// @param [in] ff FFI instance to copy
    FFI & operator=(const FFI &ff);

/// Return read identification string in file
    string getIdent() const { return _ident; }

 private:

   size_t              _msg, _string_nb;
   string              _input_file, _token, _ident, _prompt;
   char                _buffer[121];
   bool                _in, _non_fatal, _comment, _eol;
   vector<string>      _kw;
   std::ifstream       *_is;
   std::istringstream  *_iss;
   int get_token();
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif
