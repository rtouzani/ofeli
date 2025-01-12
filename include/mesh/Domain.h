/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

               Definition of class Domain for domain definition

  ==============================================================================*/

#ifndef __DOMAIN_H
#define __DOMAIN_H

#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

#include "OFELI_Config.h"
#include "linear_algebra/Point.h"
#include "mesh/Mesh.h"
#include "io/FFI.h"
#include "io/Fct.h"

using std::string;
using std::vector;
using std::ofstream;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
int bamg1(string emfile, string outfile);
int main_bamg(string emfile, string outfile);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Domain.h
 *  \brief Definition file for class Domain.
 */

/** \class Domain
 *  \ingroup Mesh
 *  \brief To store and treat finite element geometric information.
 *
 * \details This class is essentially useful to construct data for mesh generators.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Domain
{

struct Vertex : public Point<real_t> {
   real_t h;
   size_t label;
   int code;
   Vertex() : h(0.1), label(1), code(0) { }
};

struct Ln {
   Ln() : nb(0), n1(0), n2(0), Dcode(0), Ncode(0) { }
   size_t nb, n1, n2;
   int Dcode, Ncode;
   vector<real_t> s;
   vector<Point<real_t> > node;
};

struct Cont {
   size_t nb, first_line;
   vector<int> orientation;
   vector<size_t> line;
   Cont() : nb(0), first_line(0) { }
};

struct Sd { int code, contour, line, type, orient; }d;
struct LP { size_t i, j; int dc, nc; };
struct El { size_t n1, n2, n3, n4; int code; };

 public:

/// \brief Constructor of a null domain.
/// \details This constructor assigns maximal values of parameters.
    Domain();

/// \brief Constructor with an input file.
///  @param [in] file Input file in the XML format defining the domain
    Domain(const string& file);

/// \brief Destructor
    ~Domain();

/// \brief Set file containing Domain data
    void setFile(string file) { _domain_file = file; }

/// \brief Set space dimension
    void setDim(size_t d) { _dim = d; }
   
/// \brief Return space dimension
    size_t getDim() const { return _dim; }
   
/// \brief Set number of degrees of freedom
    void setNbDOF(size_t n) { _nb_dof = n; }
   
/// \brief Return number of degrees of freedom
    size_t getNbDOF() const { return _nb_dof; }

/// \brief Return number of vertices
    size_t getNbVertices() const { return _nb_vertices; }

/// \brief Return number of lines
    size_t getNbLines() const { return _nb_lines; }

/// \brief Return number of contours
    size_t getNbContours() const { return _nb_contours; }

/// \brief Return number of holes
    size_t getNbHoles() const { return _nb_holes; }

/// \brief Return number of sub-domains
    size_t getNbSubDomains() const { return _nb_sub_domains; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Get() { return get(); }
    void Get(const string& file) { get(file); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Read domain data interactively
    int get();

/// \brief Read domain data from a data file
/// @param [in] file Input file in Domain XML format
    void get(const string& file);

/// \brief Return reference to generated Mesh instance
    Mesh &getMesh() const { return *_theMesh; }
    
/// \brief Generate geometry file
    void genGeo(string file);

/// \brief Generate 2-D mesh
    void genMesh();
   
/// \brief Generate 2-D mesh and save in file (OFELI format)
/// @param [in] file File where the generated mesh is saved
    void genMesh(const string& file);

/** \brief Generate 2-D mesh and save geo, bamg and mesh file (OFELI format)
 * @param [in] geo_file Geo file
 * @param [in] bamg_file Bamg file
 * @param [in] mesh_file File where the generated mesh is saved
 */
    void genMesh(string geo_file,
                 string bamg_file,
                 string mesh_file);

/// \brief Generate 2-D mesh using the BAMG mesh generator
    void generateMesh();

/// \brief Operator <tt>*=</tt>
/// \details Rescale domain coordinates by myltiplying by a factor
/// @param [in] a Value to multiply by
    Domain &operator*=(real_t a);

/** \brief Insert a vertex
 *  @param [in] x x-coordinate of vertex
 *  @param [in] y y-coordinate of vertex
 *  @param [in] h mesh size around vertex
 *  @param [in] code code of coordinate
 */
    void insertVertex(real_t x,
                      real_t y,
                      real_t h,
                      int code);

/** \brief Insert a vertex (3-D case)
 *  @param [in] x x-coordinate of vertex
 *  @param [in] y y-coordinate of vertex
 *  @param [in] z z-coordinate of vertex
 *  @param [in] h mesh size around vertex
 *  @param [in] code code of coordinate
 */
    void insertVertex(real_t x,
                      real_t y,
                      real_t z,
                      real_t h,
                      int code);

/** \brief Insert a straight line
 *  @param [in] n1 Label of the first vertex of line
 *  @param [in] n2 Label of the second vertex of line
 *  @param [in] c Code to associate to created nodes (Dirichlet) or sides (Neumann) if < 0
 */
    void insertLine(size_t n1,
                    size_t n2,
                    int    c);

/** \brief Insert a straight line
 *  @param [in] n1 Label of the first vertex of line
 *  @param [in] n2 Label of the second vertex of line
 *  @param [in] dc Code to associate to created nodes (Dirichlet)
 *  @param [in] nc Code to associate to created sides (Neumann)
 */
    void insertLine(size_t n1,
                    size_t n2,
                    int    dc,
                    int    nc);

/** \brief Insert a circluar arc
 *  @param [in] n1 Label of vertex defining the first end of the arc
 *  @param [in] n2 Label of vertex defining the second end of the arc
 *  @param [in] n3 Label of vertex defining the center of the arc
 *  @param [in] c Code to associate to created nodes (Dirichlet) or sides (Neumann) if < 0
 */
    void insertCircle(size_t n1,
                      size_t n2,
                      size_t n3,
                      int    c);

/** \brief Insert a circluar arc
 *  @param [in] n1 Label of vertex defining the first end of the arc
 *  @param [in] n2 Label of vertex defining the second end of the arc
 *  @param [in] n3 Label of vertex defining the center of the arc
 *  @param [in] dc Code to associate to created nodes (Dirichlet)
 *  @param [in] nc Code to associate to created sides (Neumann)
 */
    void insertCircle(size_t n1,
                      size_t n2,
                      size_t n3,
                      int    dc,
                      int    nc);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Insert contour
/// @param [in] c Vector containing list of lines that describe the contour
    void insertContour(const vector<size_t>& c);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Insert hole
/// @param [in] c Vector containing list of lines that describe the hole
    void insertHole(const vector<size_t>& c);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Insert a required (imposed) vertex
/// @param [in] v Label of vertex
    void insertRequiredVertex(size_t v);

/// \brief Insert a required (imposed) edge (or line)
/// @param [in] e Label of line
    void insertRequiredEdge(size_t e);

/** \brief Insert subdomain
 *  @param [in] n
 *  @param [in] code
 */
    void insertSubDomain(size_t n,
                         int    code);

/** \brief Insert subdomain
 *  @param [in] ln Line label
 *  @param [in] orient Orientation (1 or -1)
 *  @param [in] code Subdomain code or reference
 */
    void insertSubDomain(size_t ln,
                         int    orient,
                         int    code);

/// \brief Set Number of degrees of freedom per node
    void setNbDOF(int nb_dof);

/// \brief Return minimum coordinates of vertices
    Point<real_t> getMinCoord() const;
   
/// \brief Return maximum coordinates of vertices
    Point<real_t> getMaxCoord() const;

/// \brief Return minimal value of mesh size
    real_t getMinh() const;

/// \brief Define output mesh file
/// @param [in] file String defining output mesh file
    void setOutputFile(string file) { _output_file = file; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int Rectangle();
    int Rectangle(real_t* x,
                  size_t  n1,
                  size_t  n2,
                  real_t  r,
                  int     c1,
                  int     c2,
                  int     c3,
                  int     c4);
    int Rectangle(real_t* x,
                  size_t  n1, 
                  size_t  n2,
                  real_t  r,
                  int*    c,
                  int*    cv);
    Mesh *getMesh() { return _theMesh; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    friend class XMLParser;

private:

   FFI *_ff;
   Mesh *_theMesh;
   string _domain_file, _geo_file, _bamg_file, _ofeli_file, _output_file;
   Fct _theFct;
   ofstream *_jf;
   const vector<string> _kw {"q$uit","he$lp","v$ertex","line$","circ$le","sub$domain","rectangle$",
                             "disk$","dv$","dl$","lis$t","save$","mesh$","dof$","dim$","curve$",
                             "con$tour","h$ole","EOF$"};
   vector<Vertex> _v, _nd, _bnd;
   vector<El> _el;
   vector<Ln> _l;
   vector<Cont> _h;
   vector<Sd> _sd;
   vector<Cont> _c;
   vector<LP> _ln;
   vector<size_t> _v_label, _l_label, _h_label;
   vector<size_t> _required_vertex, _required_edge;
   vector<int> _mark;
   size_t _dim, _nb_vertices, _nb_lines, _nb_contours, _nb_holes, _sub_domain;
   size_t _nb_dof, _nb_required_vertices, _nb_required_edges, _nb_sub_domains;
   size_t _ret_cont, _ret_line, _ret_save, _ret_sd;

   void getVertex();
   int getLine();
   int getCurve();
   void getCircle();
   int getContour();
   int getHole();
   int getSubDomain();
   int Disk();
   void deleteVertex();
   void deleteLine();
   void list();
   int saveAsEasymesh();
   int saveAsBamg();
   int saveAsTriangle();
   int Position(real_t  s,
                real_t& x,
                real_t& y,
                real_t& z);
   void gm();
   void gm1();
   void gm2();
   void gm3();
   void removeUnusedNodes();
   int remove(size_t          item,
              size_t&         length,
              vector<size_t>& set);
   void dof_code(int          mark,
                 vector<int>& code);
   void dof_code(int  mark,
                 int* code);
   int zero_code(const vector<int>& code);
   int setCode(size_t ne1,
               size_t ne2,
               size_t i,
               size_t j,
               int    c1,
               int    c2,
               int    c3,
               int    c4);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
