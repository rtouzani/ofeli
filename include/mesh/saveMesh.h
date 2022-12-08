/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

           Prototypes for functions to save mesh in various format files

  ==============================================================================*/


#ifndef __SAVE_MESH_H
#define __SAVE_MESH_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::endl;

#include <iomanip>
using std::setw;
using std::ios;
using std::setprecision;


#include "OFELI_Config.h"
#include "mesh/MeshUtil.h"
#include "util/util.h"
#include "linear_algebra/DMatrix.h"


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file saveMesh.h
 *  \brief Prototypes for functions to save mesh in various file formats.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*
struct El  {
  size_t label, type, region, nb_nodes, node[30];
  int code[6], cc;
};

struct Nd {
   size_t label;
   int code[6], cc;
   real_t x[3];
};
*/
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \fn void saveMesh(const string &file, const Mesh &mesh, ExternalFileFormat form)
 *  \ingroup Util
 *  \brief This function saves mesh data a file for a given external format
 *  @param [in] file File where to store mesh
 *  @param [in] mesh Mesh instance to save
 *  @param [in] form Format of the mesh file. This one can be chosen among the enumerated
 *  values:
 *  <ul>
 *     <li><tt>GMSH</tt>: %Mesh generator and graphical postprocessor \b Gmsh:
 *         <tt>http://www.geuz.org/gmsh/</tt>
 *     <li><tt>GNUPLOT</tt>: Well known graphics software:
 *         <tt>http://www.gnuplot.info/</tt>
 *     <li><tt>MATLAB</tt>: Matlab file:
 *         <tt>http://www.mathworks.com/products/matlab/</tt>
 *     <li><tt>TECPLOT</tt>: Commercial graphics software:
 *         <tt>http://www.tecplot.com</tt>
 *     <li><tt>VTK</tt>: Graphics format for the free postprocessor \b ParaView: 
 *         <tt>http://public.kitware.com/VTK/</tt>
 *  </ul>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void saveMesh(const string&      file,
              const Mesh&        mesh,
              ExternalFileFormat form);


/** \fn void saveGmsh(const string &file, const Mesh &mesh)
 *  \ingroup Util
 *  \brief This function outputs a Mesh instance in a file in
 *  <a href="http://www.geuz.org/gmsh/">Gmsh</a> format.
 *
 *  \note \b Gmsh is a free mesh generator that can be downloaded from the site: http://www.geuz.org/gmsh/
 *  @param [in] file Output file in \b Gmsh format.
 *  @param [in] mesh Mesh instance to save.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void saveGmsh(const string& file,
              const Mesh&   mesh);

/** \fn void saveGnuplot(const string &file, const Mesh &mesh)
 *  \ingroup Util
 *  \brief This function outputs a Mesh instance in a file in <a href="http://www.gnuplot.info/">Gmsh</a> format.
 *
 *  \note \b Gnuplot is a command-line driven program for producing 2D and 3D plots. It is
 *  under the GNU General Public License. Available information can be found in the site:\n
 *  <tt>http://www.gnuplot.info/</tt>
 *  @param [out] file Output file in \b Gnuplot format.
 *  @param [in] mesh Mesh instance to save.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void saveGnuplot(const string& file,
                 const Mesh&   mesh);

/** \fn void saveMatlab(const string &file, const Mesh &mesh)
 *  \ingroup Util
 *  \brief This function outputs a Mesh instance in a file in <a href="http://www.mathworks.com/products/matlab/">Matlab</a> format.
 *
 *  \note \b Matlab is a language of scientific computing including visualization. It is
 *  developed by <a href="http://www.mathworks.com">MathWorks</a>. Available information can be found in the site:\n
 *  <tt>http://www.mathworks.com/products/matlab/</tt>
 *  @param [out] file Output file in \b Matlab format.
 *  @param [in] mesh Mesh instance to save.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void saveMatlab(const string& file,
                const Mesh&   mesh);

/** \fn void saveTecplot(const string &file, const Mesh &mesh)
 *  \ingroup Util
 *  \brief This function outputs a Mesh instance in a file in <a href="http://www.tecplot.com">Tecplot</a> format.
 *
 *  \note \b Tecplot is high quality post graphical commercial processing program developed
 *  by <a href="http://www.amtec.com">Amtec</a>. Available information can be found in the site:\n
 *  <tt>http://www.tecplot.com</tt>
 *  @param [out] file Output file in \b Tecplot format.
 *  @param [in] mesh Mesh instance to save.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void saveTecplot(const string& file,
                 const Mesh&   mesh);

/** \fn void saveVTK(const string& file, const Mesh& mesh)
 *  \ingroup Util
 *  \brief This function outputs a Mesh instance in a file in
 *  <a href="http://public.kitware.com/VTK/">VTK</a> format.
 *
 *  \note The Visualization ToolKit (VTK) is an open source, freely available software
 *  system for 3D computer graphics. Available information can be found in the site:\n
 *  <tt>http://public.kitware.com/VTK/</tt>
 *  @param [out] file Output file in \b VTK format.
 *  @param [in] mesh Mesh instance to save.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void saveVTK(const string& file,
             const Mesh&   mesh);

/** \fn void saveBamg(const string& file, Mesh& mesh)
 *  \ingroup Util
 *  \brief This function outputs a Mesh instance in a file in
 *  <a href="http://raweb.inria.fr/rapportsactivite/RA2002/gamma/uid25.html">Bamg</a> format.
 * @param [in] file Name of a file written in the Bamg format.\n
 * \note \b Bamg is a 2-D mesh generator. It allows to construct adapted meshes from a 
 * given metric. It was developed at INRIA, France. Available information can be found 
 * in the site:\n
 * <tt>http://raweb.inria.fr/rapportsactivite/RA2002/gamma/uid25.html</tt>
 * @param [in] mesh Mesh instance.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void saveBamg(const string& file,
              const Mesh&   mesh);
              
/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
