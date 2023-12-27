/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

         Prototypes for functions to save fields in various format files

  ==============================================================================*/


#ifndef __SAVE_FIELD_H
#define __SAVE_FIELD_H

#include "OFELI_Config.h"
#include "mesh/Grid.h"
#include "linear_algebra/Matrix.h"

#ifdef USE_PETSC
#include "linear_algebra/petsc/PETScVect.h"
#endif

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file saveField.h
 *  \brief Prototypes for functions to save mesh in various file formats.
 */

template<class T_> class Vect;

/** \fn void saveField(Vect<real_t> &v, string output_file, int opt)
 *  \ingroup Util
 *  \brief Save a vector to an output file in a given file format
 *  \details Case where the vector contains mesh information
 *  @param [in] v Vect instance to save
 *  @param [in] output_file Output file where to save the vector
 *  @param [in] opt Option to choose file format to save. This is to be chosen
 *  among enumerated values: <tt>GMSH</tt>, <tt>GNUPLOT</tt>, <tt>MATLAB</tt>, 
 *  <tt>TECPLOT</tt>, <tt>VTK</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
   void saveField(Vect<real_t>& v,
                  string        output_file,
                  int           opt);

/*! \file saveField.h
 *  \brief Prototypes for functions to save mesh in various file formats.
 */

/** \fn void saveField(const Vect<real_t> &v, const Mesh &mesh, string output_file, int opt)
 *  \ingroup Util
 *  \brief Save a vector to an output file in a given file format.
 *  \details Case where the vector does not contain mesh information
 *  @param [in] v Vect instance to save
 *  @param [in] mesh Mesh instance
 *  @param [in] output_file Output file where to save the vector
 *  @param [in] opt Option to choose file format to save. This is to be chosen
 *  among enumerated values: <tt>GMSH</tt>, <tt>GNUPLOT</tt>, <tt>MATLAB</tt>,
 *  <tt>TECPLOT</tt>, <tt>VTK</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
   void saveField(const Vect<real_t>& v,
                  const Mesh&         mesh,
                  string              output_file,
                  int                 opt);

#ifdef USE_PETSC

/** \fn void saveField(PETScVect<real_t> &v, string output_file, int opt)
 *  \ingroup Util
 *  \brief Save a PETSc vector to an output file in a given file format
 *  \details Case where the vector does not contain mesh information
 *  @param [in] v PETScVect instance to save
 *  @param [in] output_file Output file where to save the vector
 *  @param [in] opt Option to choose file format to save. This is to be chosen
 *  among enumerated values: <tt>GMSH</tt>, <tt>GNUPLOT</tt>, <tt>MATLAB</tt>,
 *  <tt>TECPLOT</tt>, <tt>VTK</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
   void saveField(PETScVect<real_t>& v,
                  string             output_file,
                  int                opt);

/** \fn void saveField(PETScVect<real_t> &v, const Mesh& mesh, string output_file, int opt)
 *  \ingroup Util
 *  \brief Save a PETSc vector to an output file in a given file format
 *  \details Case where the vector does not contain mesh information
 *  @param [in] v PETScVect instance to save
 *  @param [in] mesh Mesh instance
 *  @param [in] output_file Output file where to save the vector
 *  @param [in] opt Option to choose file format to save. This is to be chosen
 *  among enumerated values: <tt>GMSH</tt>, <tt>GNUPLOT</tt>, <tt>MATLAB</tt>,
 *  <tt>TECPLOT</tt>, <tt>VTK</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
   void saveField(PETScVect<real_t>& v,
                  const Mesh&        mesh,
                  string             output_file,
                  int                opt);
#endif

/** \fn void saveField(Vect<real_t> &v, const Grid &g, string output_file, int opt=VTK)
 *  \ingroup Util
 *  \brief Save a vector to an output file in a given file format,
 *  for a structured grid data.
 *  @param [in] v Vect instance to save
 *  @param [in] g Grid instance
 *  @param [in] output_file Output file where to save the vector
 *  @param [in] opt Option to choose file format to save. This is to be chosen
 *  among enumerated values: <tt>GMSH</tt>, <tt>VTK</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveField(Vect<real_t>& v,
                   const Grid&   g,
                   string        output_file,
                   int           opt);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void saveFormat(Mesh&  mesh,
                    string input_file,
                    string output_file,
                    int    format,
                    int    f=1);

    void getfields(string                   file,
                   Mesh&                    ms,
                   size_t&                  nb_dof,
                   vector<real_t>&          t,
                   vector<vector<real_t> >& u,
                   string&                  name
                  );
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \fn void saveGnuplot(string input_file, string output_file, string mesh_file, int f)
 *  \ingroup Util
 *  \brief Save a vector to an input <a href="http://www.gnuplot.info/">Gnuplot</a> file.
 *
 *  \details Gnuplot is a command-line driven program for producing 2D and 3D plots. It is
 *  under the GNU General Public License. Available information can be found in the site:\n
 *  http://www.gnuplot.info/
 *  @param [in] input_file Input file (OFELI XML file containing a field).
 *  @param [in] output_file Output file (gnuplot format file)
 *  @param [in] mesh_file File containing mesh data
 *  @param [in] f Field is stored each \c f time step [Default: <tt>1</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveGnuplot(string input_file,
                     string output_file,
                     string mesh_file,
                     int    f=1);

/** \fn void saveGnuplot(Mesh& mesh, string input_file, string output_file, int f)
 *  \ingroup Util
 *  \brief Save a vector to an input <a href="http://www.gnuplot.info/">Gnuplot</a> file.
 *
 *  \details Gnuplot is a command-line driven program for producing 2D and 3D plots. It is
 *  under the GNU General Public License. Available information can be found in the site:\n
 *  http://www.gnuplot.info/
 *  @param [in] mesh Reference to Mesh instance
 *  @param [in] input_file Input file (OFELI XML file containing a field).
 *  @param [in] output_file Output file (gnuplot format file)
 *  @param [in] f Field is stored each \c f time step [Default: <tt>1</tt>]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveGnuplot(Mesh&  mesh,
                     string input_file,
                     string output_file,
                     int    f=1);

/**  \fn void saveTecplot(string input_file, string output_file, string mesh_file, int f)
 *  \ingroup Util
 *  \brief Save a vector to an output file to an input <a href="http://www.tecplot.com">Tecplot</a> file.
 *
 *  \details Tecplot is high quality post graphical commercial processing program developed
 *  by \b Amtec. Available information can be found in the site: http://www.tecplot.com
 *  @param [in] input_file Input file (OFELI XML file containing a field).
 *  @param [in] output_file Output file (gnuplot format file)
 *  @param [in] mesh_file File containing mesh data
 *  @param [in] f Field is stored each \c f time step [Default: <tt>1</tt>]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveTecplot(string input_file,
                     string output_file,
                     string mesh_file,
                     int    f=1);

/**  \fn void saveTecplot(Mesh& mesh, string input_file, string output_file, int f)
 *  \ingroup Util
 *  \brief Save a vector to an output file to an input <a href="http://www.tecplot.com">Tecplot</a> file.
 *
 *  \details Tecplot is high quality post graphical commercial processing program developed
 *  by \b Amtec. Available information can be found in the site: http://www.tecplot.com
 *  @param [in] mesh Reference to Mesh instance
 *  @param [in] input_file Input file (OFELI XML file containing a field).
 *  @param [in] output_file Output file (gnuplot format file)
 *  @param [in] f Field is stored each \c f time step [Default: <tt>1</tt>]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveTecplot(Mesh&  mesh,
                     string input_file,
                     string output_file,
                     int    f=1);

/** \fn saveVTK(string input_file, string output_file, string mesh_file, int f)
 *
 *  \ingroup Util
 *  \brief Save a vector to an output <a href="http://public.kitware.com/VTK/">VTK</a> file.
 *
 *  \details The Visualization ToolKit (VTK) is an open source, freely available software
 *  system for 3D computer graphics. Available information can be found in the site:\n
 *  http://public.kitware.com/VTK/
 *  @param [in] input_file Input file (OFELI XML file containing a field).
 *  @param [in] output_file Output file (VTK format file)
 *  @param [in] mesh_file File containing mesh data
 *  @param [in] f Field is stored each \c f time step [Default: <tt>1</tt>]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveVTK(string input_file,
                 string output_file,
                 string mesh_file,
                 int    f=1);

/** \fn saveVTK(Mesh& mesh, string input_file, string output_file, int f)
 *
 *  \ingroup Util
 *  \brief Save a vector to an output <a href="http://public.kitware.com/VTK/">VTK</a> file.
 *
 *  \details The Visualization ToolKit (VTK) is an open source, freely available software
 *  system for 3D computer graphics. Available information can be found in the site:\n
 *  http://public.kitware.com/VTK/
 *  @param [in] mesh Reference to Mesh instance
 *  @param [in] input_file Input file (OFELI XML file containing a field).
 *  @param [in] output_file Output file (VTK format file)
 *  @param [in] f Field is stored each \c f time step [Default: <tt>1</tt>]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveVTK(Mesh&  mesh,
                 string input_file,
                 string output_file,
                 int    f=1);

/**  \fn void saveGmsh(string input_file, string output_file, string mesh_file, int f)
 *
 *  \ingroup Util
 *  \brief Save a vector to an output <a href="http://www.geuz.org/gmsh/">Gmsh</a> file.
 *  \details  <a href="http://www.geuz.org/gmsh/">Gmsh</a> is a free mesh generator and 
 *  postprocessor that can be downloaded from the site:\n
 *  http://www.geuz.org/gmsh/
 *  @param [in] input_file Input file (OFELI XML file containing a field).
 *  @param [in] output_file Output file (Gmsh format file)
 *  @param [in] mesh_file File containing mesh data
 *  @param [in] f Field is stored each \c f time step [Default: <tt>1</tt>]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveGmsh(string input_file,
                  string output_file,
                  string mesh_file,
                  int    f=1);

/**  \fn void saveGmsh(Mesh& mesh, string input_file, string output_file, int f)
 *
 *  \ingroup Util
 *  \brief Save a vector to an output <a href="http://www.geuz.org/gmsh/">Gmsh</a> file.
 *  \details  <a href="http://www.geuz.org/gmsh/">Gmsh</a> is a free mesh generator and 
 *  postprocessor that can be downloaded from the site:\n
 *  http://www.geuz.org/gmsh/
 *  @param [in] mesh Reference to Mesh instance
 *  @param [in] input_file Input file (OFELI XML file containing a field).
 *  @param [in] output_file Output file (Gmsh format file)
 *  @param [in] f Field is stored each \c f time step [Default: <tt>1</tt>]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    void saveGmsh(Mesh&  mesh,
                  string input_file,
                  string output_file,
                  int    f=1);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void saveMatrix(const Matrix<real_t>* A,
                    string                file);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
