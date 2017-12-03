/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

              Definition of class 'IOField' to manage XML Field files

  ==============================================================================*/

#ifndef __IO_FIELD_H
#define __IO_FIELD_H

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using std::ostream;
using std::endl;
using std::clog;
using std::ifstream;
using std::ofstream;
using std::setw;
using std::ios;
using std::setprecision;

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/DMatrix.h"
#include "linear_algebra/DSMatrix.h"
#include "io/XMLParser.h"

#ifdef USE_PETSC
#include "linear_algebra/petsc/PETScVect.h"
#endif


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file IOField.h
 *  \brief Definition file for class IOField.
 */

/** \defgroup IO Input/Output
 *  \brief Input/Output utility classes
 */

/*! \class IOField
 *  \ingroup IO
 * \brief Enables working with files in the XML Format.
 * \details This class has methods to store vectors in files and read from files.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class IOField : public XMLParser
{

 public:

/// Enumerated values for file access type
    enum AccessType {
       IN  = 1,
       OUT = 2
    };

/// Default constructor
    IOField();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    IOField(const string& file,
            char*         access);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Constructor using file name.
 *  @param [in] file File name.
 *  @param [in] access Access code. This number is to be chosen among two enumerated values:
 *  <ul>
 *    <li> <tt>IOField::IN</tt> to read the file
 *    <li> <tt>IOField::OUT</tt> to write on it
 *  </ul>
 *  @param [in] compact Flag to choose a compact storage or not [Default: <tt>true</tt>]
 */
    IOField(const string& file,
            AccessType    access,
            bool          compact=true);

/** \brief Constructor using file name, mesh file and mesh.
 *  @param [in] mesh_file File containing mesh
 *  @param [in] file File that contains field stored or to store
 *  @param [in] ms Mesh instance
 *  @param [in] access Access code. This number is to be chosen among two enumerated values:
 *  <ul>
 *    <li> <tt>IOField::IN</tt> to read the file
 *    <li> <tt>IOField::OUT</tt> to write on it
 *  </ul>
 *  @param [in] compact Flag to choose a compact storage or not [Default: <tt>true</tt>]
 */
    IOField(const string& mesh_file,
            const string& file,
            Mesh&         ms,
            AccessType    access,
            bool          compact=true);

/** \brief Constructor using file name and mesh.
 *  @param [in] file File that contains field stored or to store
 *  @param [in] ms Mesh instance
 *  @param [in] access Access code. This number is to be chosen among two enumerated values:
 *  <ul>
 *    <li> <tt>IOField::IN</tt> to read the file
 *    <li> <tt>IOField::OUT</tt> to write on it
 *  </ul>
 *  @param [in] compact Flag to choose a compact storage or not [Default: <tt>true</tt>]
 */
    IOField(const string& file,
            Mesh&         ms,
            AccessType    access,
                  bool    compact=true);

/** \brief Constructor using file name and field name.
 *  @param [in] file File that contains field stored or to store
 *  @param [in] access Access code. This number is to be chosen among two enumerated values:
 *  <ul>
 *    <li> <tt>IOField::IN</tt> to read the file
 *    <li> <tt>IOField::OUT</tt> to write on it
 *  </ul>
 *  @param [in] name Seek a specific field with given \a name
 */
    IOField(const string& file,
            AccessType    access,
            const string& name);

/// \brief Destructor
    ~IOField();

/// \brief Set mesh file
/// @param [in] file Mesh file
    void setMeshFile(const string& file);

/// \brief Open file
/// \details Case where file name has been previously given (in the constructor).
    void open();

/** \brief Open file
 *  @param [in] file File name.
 *  @param [in] access Access code. This number is to be chosen among two enumerated values:
 *  <ul>
 *    <li> <tt>IOField::IN</tt> to read the file
 *    <li> <tt>IOField::OUT</tt> to write on it
 *  </ul>
 */
    void open(const string& file,
              AccessType    access);

/// \brief Close file.
    void close();

/// \brief Store mesh in file.
    void put(Mesh& ms);

/// \brief Store Vect instance <tt>v</tt> in file
/// @param [in] v Vect instance to store
    void put(const Vect<real_t>& v);

#ifdef USE_PETSC
/// \brief Store PETScVect instance <tt>v</tt> in file
/// @param [in] v PETScVect instance to store
    void put(const PETScVect<real_t>& v);
#endif

/// \brief Get Vect <tt>v</tt> instance from file.
/// \details First time step is read from the <tt>XML</tt> file.
    real_t get(Vect<real_t>& v);

/** \brief Get Vect <tt>v</tt> instance from file if the field has the given name.
 *  \details First time step is read from the <tt>XML</tt> file.
 *  @param [in,out] v Vect instance
 *  @param [in] name Name to seek in the XML file
 */
    int get(Vect<real_t>& v,
            const string& name);

/** \brief Get DMatrix <tt>A</tt> instance from file if the field has the given name.
 *  \details First time step is read from the <tt>XML</tt> file.
 *  @param [in,out] A DMatrix instance
 *  @param [in] name Name to seek in the XML file
 */
    int get(DMatrix<real_t>& A,
            const string&    name);

/** \brief Get DSMatrix <tt>A</tt> instance from file if the field has the given name.
 *  \details First time step is read from the <tt>XML</tt> file.
 *  @param [in,out] A DSMatrix instance
 *  @param [in] name Name to seek in the XML file
 */
    int get(DSMatrix<real_t>& A,
            const string&     name);

/** \brief Get Vect <tt>v</tt> instance from file corresponding to a specific time value.
 *  \details The sought vector corresponding to the time value is read from the <tt>XML</tt> file.
 *  @param [in,out] v Vector instance
 *  @param [in] t Time value
 */
    int get(Vect<real_t>& v,
            real_t        t);

/** \brief Save field vectors in a file using \b GMSH format
 *  \details This member function enables avoiding the use of <tt>cfield</tt>. It
 *  must be used once all field vectors have been stored in output file. It 
 *  closes this file and copies its contents to a \b GMSH file.
 *  @param [in] output_file Output file name where to store using \b GMSH format
 *  @param [in] mesh_file File containing mesh data
 */
    void saveGMSH(string output_file,
                  string mesh_file);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int get(Mesh&                    ms,
            vector<vector<real_t> >& v,
            string&                  name);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 private:

   ofstream  *_of;
   string    _field_name;
   int       _state;
   bool      _field_opened, _compact, _no_mesh_file;
   Mesh      *_theMesh;
};

} /* namespace OFELI */

#endif
