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

                 Functions to read mesh files in various formats
                          and construct Mesh instance

  ==============================================================================*/


#ifndef __GET_MESH_H
#define __GET_MESH_H

#include <string>
using std::string;

#include <vector>
#include "OFELI_Config.h"
#include "linear_algebra/Point.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file getMesh.h
 *  \brief Definition file for mesh conversion functions.
 */

class Mesh;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
void dof_code(size_t nb_dof, int mark, int* code);
struct Entity { size_t nb; int type; std::vector<int> l; };
struct El { size_t dim, n, region, nb_nodes, node[30]; int shape, code, cc; };
struct Nd {
   size_t n;
   int code[6], cc;
   Point<real_t> x;
   Nd() : cc(0) { }
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \fn void getMesh(string file, ExternalFileFormat form, Mesh& mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a mesh file stored in
 *  an external file format.
 *  @param [in] file Input mesh file name.\n
 *  @param [in] form Format of the mesh file. This one can be chosen among the enumerated
 *  values:
 *  <ul>
 *     <li><tt>GMSH</tt>: %Mesh generator \b Gmsh, see site:\n
 *         <tt>http://www.geuz.org/gmsh/</tt>
 *     <li><tt>MATLAB</tt>: Matlab file, see site:\n
 *         <tt>http://www.mathworks.com/products/matlab/</tt>
 *     <li><tt>EASYMESH</tt>: \b Easymesh is a 2-D mesh generator, see site:\n
 *         <tt>http://web.mit.edu/easymesh_v1.4/www/easymesh.html</tt>
 *     <li><tt>GAMBIT</tt>: \b Gambit is a mesh generator associated to \b Fluent
 *         <tt>http://www.stanford.edu/class/me469b/gambit_download.html</tt>
 *     <li><tt>BAMG</tt>: %Mesh generator <tt>Bamg</tt>, see site:\n
 *         <tt>http://raweb.inria.fr/rapportsactivite/RA2002/gamma/uid25.html</tt>
 *     <li><tt>NETGEN</tt>: \b Netgen is a 3-D mesh generator, see site:\n
 *         <tt>http://www.hpfem.jku.at/netgen/</tt>
 *     <li><tt>TETGEN</tt>: \b Tetgen is a 3-D mesh generator, see site:\n
 *         <tt>http://tetgen.berlios.de/</tt>
 *     <li><tt>TRIANGLE_FF</tt>: \b Triangle is a 2-D mesh generator, see site:\n  
 *         <tt>http://www.cs.cmu.edu/~quake/triangle.html</tt>
 *  </ul>
 *  @param [out] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *  is not provided, in general, by mesh generators. Its default value here is 1.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getMesh(string             file,
             ExternalFileFormat form,
             Mesh&              mesh,
             size_t             nb_dof=1);

/** \fn void getBamg(string file, Mesh& mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a mesh file stored in
 *  <a href="www-rocq1.inria.fr/gamma/cdrom/www/bamg/eng.htm">Bamg</a> format.
 *
 *  @param [in] file Name of a file written in the Bamg format.\n
 *  \note \b Bamg is a 2-D mesh generator. It allows to construct adapted meshes from a 
 *   given metric. It was developed at INRIA, France. Available information can be found 
 *  in the site:\n
 *  <tt>http://raweb.inria.fr/rapportsactivite/RA2002/gamma/uid25.html</tt>
 *  @param [out] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *  is not provided, in general, by mesh generators. Its default value here is 1.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getBamg(string file,
             Mesh&  mesh,
             size_t nb_dof=1);

/** \fn void getEasymesh(string file, Mesh &mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a mesh file stored in
 *  <a href="http://www-dinma.univ.trieste.it/nirftc/research/easymesh/Default.htm">Easymesh</a> format.
 *
 *  @param [in] file Name of a file (without extension) written in \b Easymesh format.
 *  Actually, the function Easymesh2MDF attempts to read mesh data from files <tt>file.e</tt>,
 *  <tt>file.n</tt> and <tt>file.s</tt> produced by \b Easymesh. 
 *  
 *  \note \b Easymesh is a free program that generates 2-D, unstructured, Delaunay and 
 *  constrained Delaunay triangulations in general domains. It can be downloaded from the site:\n
 *  http://www-dinma.univ.trieste.it/nirftc/research/easymesh/Default.htm
 *  @param [in] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *  is not provided, in general, by mesh generators. Its default value here is 1.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getEasymesh(string file,
                 Mesh&  mesh,
                 size_t nb_dof=1);

/** \fn void getGambit(string file, Mesh &mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a mesh file stored in \b Gambit neutral format.
 *
 *  \note \b Gambit is a commercial mesh generator associated to the CFD code
 *  <a href="http://www.fluent.com/software/gambit/">Fluent</a>.
 *  Informations about \b Gambit can be found in the site:\n
 *  http://www.fluent.com/software/gambit/
 *  @param [in] file Name of a file written in the \b Gambit neutral format.
 *  @param [out] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *  is not provided, in general, by mesh generators. Its default value here is 1.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getGambit(string file,
               Mesh&  mesh,
               size_t nb_dof=1);

/** \fn void getGmsh(string file, Mesh &mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a mesh file stored in <a href="http://www.geuz.org/gmsh/">Gmsh</a> format.
 *
 *  \note \b Gmsh is a free mesh generator that can be downloaded from the site:\n
 *  http://www.geuz.org/gmsh/
 *  @param [in] file Name of a file written in the \b Gmsh format.
 *  @param [out] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *  is not provided, in general, by mesh generators. Its default value here is 1.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getGmsh(string file,
             Mesh&  mesh,
             size_t nb_dof=1);

/** \fn void getMatlab(string file, Mesh &mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a Matlab mesh data.
 *
 *  \note \b Matlab is a language of scientific computing including visualization. It is
 *  developed by <a href="http://www.mathworks.com">MathWorks</a>. Available information can be found in the site:\n
 *  http://www.mathworks.com/products/matlab/
 *  @param [in] file Name of a file created by Matlab by executing the script file <tt>Matlab2OFELI.m</tt>
 *  @param [out] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *  is not provided, in general, by mesh generators. Its default value here is 1.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getMatlab(string file,
               Mesh&  mesh,
               size_t nb_dof=1);

/** \fn void getNetgen(string file, Mesh &mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a mesh file stored in <a href="http://www.hpfem.jku.at/netgen/">Netgen</a> format.
 *
 *  \note \b Netgen is a tetrahedral mesh generator that can be downloaded from the site:\n
 *  http://www.hpfem.jku.at/netgen/
 *  @param [in] file Name of a file written in the Netgen format.
 *  @param [out] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *   is not provided, in general, by mesh generators. [ default = 1 ]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getNetgen(string file,
               Mesh&  mesh,
               size_t nb_dof=1);

/** \fn void getTetgen(string file, Mesh &mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a mesh file stored in
 *  <a href="http://tetgen.berlios.de/">Tetgen</a> format.
 *
 *  \note \b Tetgen is a free three-dimensional mesh generator that can be downloaded in the site:\n
 *  http://tetgen.berlios.de/
 *  @param [in] file Name of a file written in the \b Tetgen format.
 *  @param [out] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *  is not provided, in general, by mesh generators. Its default value here is 1.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getTetgen(string file,
               Mesh&  mesh,
               size_t nb_dof=1);

/** \fn void getTriangle(string file, Mesh &mesh, size_t nb_dof=1)
 *  \ingroup Util
 *  \brief Construct an instance of class Mesh from a mesh file stored in
 *  <a href="http://people.scs.fsu.edu/~burkardt/c_src/triangle/triangle.html">Triangle</a> format.
 *
 *  \note \b TRIANGLE is a C program that can generate meshes, Delaunay triangulations and Voronoi 
 *  diagrams for 2D pointsets that can be downloaded in the site:\n
 *  http://people.scs.fsu.edu/~burkardt/c_src/triangle/triangle.html/
 *  @param [in] file Name of a file written in the \b Tetgen format.
 *  @param [out] mesh Mesh instance created by the function.
 *  @param [in] nb_dof Number of degrees of freedom for each node. This information
 *  is not provided, in general, by mesh generators. Its default value here is 1.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void getTriangle(string file,
                 Mesh&  mesh,
                 size_t nb_dof=1);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
