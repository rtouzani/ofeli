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

                    Definition of 'MeshUtil' for Mesh utilities

  ==============================================================================*/


#ifndef __MESH_UTIL_H
#define __MESH_UTIL_H

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#ifdef USE_PETSC
#include "linear_algebra/petsc/PETScVect.h"
#endif

using std::pair;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file MeshUtil.h
 *  \brief Definitions of utility functions for meshes.
 */

#ifdef USE_PETSC
template<class T_> class PETScVect;
#endif

class Mesh;
class Node;
class Element;
class Side;
class Edge;
class Grid;


#ifndef DOXYGEN_SHOULD_SKIP_THIS
enum { MIN_X, MAX_X, MIN_Y, MAX_Y, MIN_Z, MAX_Z };
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn size_t Label(const Node& nd)
 *  \ingroup Mesh
 *  \brief Return label of a given node
 *  @param [in] nd Reference to Node instance
 *  @return Label of node
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
size_t Label(const Node& nd);

/** \fn size_t Label(const Element& el)
 *  \ingroup Mesh
 * \brief Return label of a given element
 *  @param [in] el Reference to Element instance
 *  @return Label of element
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
size_t Label(const Element& el);

/** \fn size_t Label(const Side& sd)
 *  \ingroup Mesh
 *  \brief Return label of a given side
 *  @param [in] sd Reference to Side instance
 *  @return Label of side
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
size_t Label(const Side& sd);

/** \fn size_t Label(const Edge& ed)
 *  \ingroup Mesh
 *  \brief Return label of a given edge
 *  @param [in] ed Reference to Edge instance
 *  @return Label of edge
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
size_t Label(const Edge& ed);

/** \fn size_t NodeLabel(const Element& el, size_t n)
 *  \ingroup Mesh
 *  \brief Return global label of node local label in element
 *  @param [in] el Reference to Element instance
 *  @param [in] n Local label of node in element
 *  @return Global label of node
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
size_t NodeLabel(const Element& el, size_t n);

/** \fn size_t NodeLabel(const Side& sd, size_t n)
 *  \ingroup Mesh
 *  \brief Return global label of node local label in side
 *  @param [in] sd Reference to Side instance
 *  @param [in] n Local label of node in side
 *  @return Global label of node
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
size_t NodeLabel(const Side& sd, size_t n);

/** \fn Point<real_t> Coord(const Node& nd)
 *  \ingroup Mesh
 * \brief Return coordinates of a given node
 *  @param [in] nd Reference to Node instance
 *  @return Coordinates of node
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
Point<real_t> Coord(const Node& nd);

/** \fn int Code(const Node& nd, size_t i=1)
 *  \ingroup Mesh
 *  \brief Return code of a given (degree of freedom of) node
 *  @param [in] nd Reference to Node instance
 *  @param [in] i Label of dof [Default: <tt>1</tt>]
 *  @return Code of dof of node
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
int Code(const Node& nd, size_t i=1);

/** \fn int Code(const Element& el)
 *  \ingroup Mesh
 *  \brief Return code of a given element
 *  @param [in] el Reference to Element instance
 *  @return Code of element
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
int Code(const Element& el);

/** \fn int Code(const Side& sd, size_t i=1)
 *  \ingroup Mesh
 *  \brief Return code of a given (degree of freedom of) side
 *  @param [in] sd Reference to Side instance
 *  @param [in] i Label of dof [Default: <tt>1</tt>]
 *  @return Code of dof of side
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
int Code(const Side& sd, size_t i=1);

/** \fn operator==(const Element& el1, const Element& el2)
 *  \ingroup Mesh
 *  \brief Check equality between 2 elements
 *  @param [in] el1 Reference to first Side instance
 *  @param [in] el2 Reference to second Side instance
 *  @return <tt>true</tt> is elements are equal, <i>i.e.</i> if they have the same nodes,
 *  <tt>false</tt> if not.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
bool operator==(const Element& el1,
                const Element& el2);

/** \fn bool operator==(const Side& sd1, const Side& sd2)
 *  \ingroup Mesh
 *  \brief Check equality between 2 sides
 *  @param [in] sd1 Reference to first Side instance
 *  @param [in] sd2 Reference to second Side instance
 *  @return <tt>true</tt> is sides are equal, <i>i.e.</i> if they have the same nodes,
 *  <tt>false</tt> if not.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
bool operator==(const Side& sd1,
                const Side& sd2);

/** \fn void DeformMesh(Mesh& mesh, const Vect<real_t>& u, real_t a=1)
 *  \brief Calculate deformed mesh using a displacement field
 *  \ingroup Mesh
 *  @param [in,out] mesh Mesh instance. On output, node coordinates are modified to take into
 *  account the displacement
 *  @param [in] u Displacement field at nodes
 *  @param [in] a Maximal deformation rate. [Default: <tt>1</tt>]. A typical value is
 *  0.2 (<i>i.e.</i> 20%).
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void DeformMesh(Mesh&               mesh,
                const Vect<real_t>& u,
                real_t              rate=0.2);

#ifdef USE_PETSC
/** \fn void DeformMesh(Mesh& mesh, const PETScVect<real_t>& u, real_t a=1)
 *  \brief Calculate deformed mesh using a displacement field as instance of PETScVect
 *  \ingroup Mesh
 *  @param [in,out] mesh Mesh instance. On output, node coordinates are modified to take into
 *  account the displacement
 *  @param [in] u Displacement field at nodes
 *  @param [in] a Amplification factor [Default: <tt>0.2</tt>].
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void DeformMesh(Mesh&                    mesh,
                const PETScVect<real_t>& u,
                real_t                   a=1);
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct ND { 
   ND() : e1(0), e2(0) { nd[0]=nd[1]=nd[2]=nd[3]=0; }
   ND(size_t n1, size_t n2, size_t n3=0, size_t n4=0) : e1(0), e2(0)
   { 
      nd[0] = n1; nd[1] = n2; nd[2] = n3; nd[3] = n4; 
   }
   size_t n, e1, e2, nd[4];
};

ostream& operator<<(ostream& s,
                    ND&      d);
bool operator==(const ND& n1,
                const ND& n2);
void DOFCode(int    mark,
             size_t nb_dof,
             int*   code);
int isReferencedSide(const Side& sd);
void HexahedraToTetrahedra(Mesh &m1, Mesh &m2);
void QuadrilateralsToTriangles(Mesh &m1, Mesh &m2);
void Refine(Mesh &in_mesh, Mesh &out_mesh);
int BoundaryConditionCode(vector<string> &str, int *code, string s);
size_t init_side_node_numbering(int shape, vector<vector<size_t> >& nsd, int &sh);
int equal_sides(const class Side *sd1, const class Side *sd2);
int equal_sides(const class Side *sd, vector<size_t> &s);
void order_side_nodes(size_t ns, ND &s);
void order_edge_nodes(ND &s);
bool compare_sides(const ND &s1, const ND &s2);
bool compare_edges(const ND &s1, const ND &s2);
void order_edge_nodes(pair<size_t,size_t> &s);
void complete_sides(vector<ND> &p);
size_t clean_edges(vector<pair<size_t,size_t> > &p, vector<pair<size_t,size_t> > &q, const size_t &n);
size_t remove_internal_sides(vector<ND> &p);
void DofCode(int mark, size_t nb_dof, int *code);
void FindRoot(size_t &root, vector<long> &xadj, vector<size_t> &adjncy,
              vector<size_t> &mask, long &nlvl, vector<size_t> &xls, size_t *ls);
void RootLs(size_t &root, vector<long> &xadj, vector<size_t> &adjncy,
            vector<size_t> &mask, long &nlvl, vector<size_t> &xls, size_t *ls);
void RCM(size_t &root, vector<long> &xadj, vector<size_t> &adjncy,
         vector<size_t> &mask, size_t *perm, size_t &ccsize, vector<size_t> &deg);
void Degree(size_t &root, vector<long> &xadj, vector<size_t> &adjncy,
            vector<size_t> &mask, vector<size_t> &deg, size_t &ccsize,
            size_t *ls);
void MeshToGrid(Mesh &m, Grid &g, const Vect<real_t> &u, Vect<real_t> &ug, size_t dof=1);
void GridToMesh(Grid &g, Mesh &m, const Vect<real_t> &ug, Vect<real_t> &u, size_t dof=1);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \fn void MeshToMesh(Mesh &m1, Mesh &m2, const Vect<real_t> &u1, Vect<real_t> &u2, 
 *                      size_t nx, size_t ny=1, size_t nz=1, size_t dof=1)
 * \ingroup Mesh
 * \brief Function to redefine a vector defined on a mesh to a new mesh
 *
 * \details The program interpolates (piecewise linear) first the vector on a finer structured
 * grid. Then the values on the new mesh nodes are computed.
 *
 * \remark For efficiency the number of grid cells must be large enough so that interpolation provides
 * efficient accuracy
 *
 * @param [in] m1 Reference to the first mesh instance
 * @param [out] m2 Reference to the second mesh instance
 * @param [in] u1 Input vector of nodal values defined on first mesh
 * @param [out] u2 Output vector of nodal values defined on second mesh
 * @param [in] nx Number of cells in the <tt>x</tt>-direction in the fine structured grid
 * @param [in] ny Number of cells in the <tt>y</tt>-direction in the fine structured grid
 *                The default value of <tt>ny</tt> is <tt>0</tt>, i.e. a 1-D grid
 * @param [in] nz Number of cells in the <tt>z</tt>-direction in the fine structured grid
 *                The default value of <tt>nz</tt> is <tt>0</tt>, i.e. a 1-D or 2-D grid
 * @param [in] dof Label of degree of freedom of vector <tt>u</tt>. Only this dof is considered.
 *                 [Default: <tt>1</tt>]
 *
 * @note The input vector <tt>u1</tt> is a one degree of freedom per node vector, i.e. its
 * size must be equal (or greater than) the total number of nodes of mesh <tt>m1</tt>.
 * The size of vector <tt>u2</tt> is deduced from the mesh <tt>m2</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void MeshToMesh(Mesh&               m1,
                Mesh&               m2,
                const Vect<real_t>& u1,
                Vect<real_t>&       u2,
                size_t              nx,
                size_t              ny=0,
                size_t              nz=0,
                size_t              dof=1);


/** \fn void MeshToMesh(const Vect<real_t> &u1, Vect<real_t> &u2, 
 *                      size_t nx, size_t ny=1, size_t nz=1, size_t dof=1)
 * \ingroup Mesh
 * \brief Function to redefine a vector defined on a mesh to a new mesh
 *
 * \details The program interpolates (piecewise linear) first the vector on a finer structured
 * grid. Then the values on the new mesh nodes are computed.
 *
 * \remark For efficiency the number of grid cells must be large enough so that interpolation provides
 * efficient accuracy
 *
 * @param [in] u1 Input vector of nodal values defined on first mesh. This vector instance must
 *             contain Mesh instance
 * @param [out] u2 Output vector of nodal values defined on second mesh. This vector instance must
 *             contain Mesh instance
 * @param [in] nx Number of cells in the <tt>x</tt>-direction in the fine structured grid
 * @param [in] ny Number of cells in the <tt>y</tt>-direction in the fine structured grid
 *                The default value of <tt>ny</tt> is <tt>0</tt>, i.e. a 1-D grid
 * @param [in] nz Number of cells in the <tt>z</tt>-direction in the fine structured grid
 *                The default value of <tt>nz</tt> is <tt>0</tt>, i.e. a 1-D or 2-D grid
 * @param [in] dof Label of degree of freedom of vector <tt>u</tt>. Only this dof is considered.
 *                 [Default: <tt>1</tt>]
 *
 * @note The input vector <tt>u1</tt> is a one degree of freedom per node vector, i.e. its
 * size must be equal (or greater than) the total number of nodes of mesh <tt>m1</tt>.
 * The size of vector <tt>u2</tt> is deduced from the mesh <tt>m2</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void MeshToMesh(const Vect<real_t>& u1,
                Vect<real_t>&       u2,
                size_t              nx,
                size_t              ny=0,
                size_t              nz=0,
                size_t              dof=1);

/** \fn void MeshToMesh(Mesh &m1, Mesh &m2, const Vect<real_t> &u1, Vect<real_t> &u2, 
 *                      const Point<real_t> &xmin, const Point<real_t> &xmax, size_t nx, 
 *                      size_t ny, size_t nz, size_t dof=1)
 * \ingroup Mesh
 * \brief Function to redefine a vector defined on a mesh to a new mesh
 *
 * \details The program interpolates (piecewise linear) first the vector on a finer structured
 * grid. Then the values on the new mesh nodes are computed. In this function the grid rectangle is defined
 * so that this one can cover only a submesh of <tt>m1</tt>.
 *
 * \remark For efficiency the number of grid cells must be large enough so that interpolation provides
 * efficient accuracy
 *
 * @param [in] m1 Reference to the first mesh instance
 * @param [out] m2 Reference to the second mesh instance
 * @param [in] u1 Input vector of nodal values defined on first mesh
 * @param [out] u2 Output vector of nodal values defined on second mesh
 * @param [in] xmin Point instance containing minimal coordinates of the rectangle that defines the grid
 * @param [in] xmax Point instance containing maximal coordinates of the rectangle that defines the grid
 * @param [in] nx Number of cells in the <tt>x</tt>-direction in the fine structured grid
 * @param [in] ny Number of cells in the <tt>y</tt>-direction in the fine structured grid
 *                The default value of <tt>ny</tt> is <tt>0</tt>, i.e. a 1-D grid
 * @param [in] nz Number of cells in the <tt>z</tt>-direction in the fine structured grid
 *                The default value of <tt>nz</tt> is <tt>0</tt>, i.e. a 1-D or 2-D grid
 * @param [in] dof Label of degree of freedom of vector <tt>u</tt>. Only this dof is considered.
 *                 [Default: <tt>1</tt>]
 *
 * @note The input vector <tt>u1</tt> is a one degree of freedom per node vector, i.e. its
 * size must be equal (or greater than) the total number of nodes of mesh <tt>m1</tt>.
 * The size of vector <tt>u2</tt> is deduced from the mesh <tt>m2</tt>
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void MeshToMesh(Mesh&                m1,
                Mesh&                m2,
                const Vect<real_t>&  u1,
                Vect<real_t>&        u2,
                const Point<real_t>& xmin,
                const Point<real_t>& xmax,
                size_t               nx,
                size_t               ny,
                size_t               nz,
                size_t               dof=1);

/** \fn real_t getMaxSize(const Mesh &m)
 * \ingroup Mesh
 * \brief Return maximal size of element edges for given mesh
 * @param [in] m Reference to mesh instance
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
real_t getMaxSize(const Mesh& m);

/** \fn real_t getMinSize(const Mesh& m)
 * \ingroup Mesh
 * \brief Return minimal size of element edges for given mesh
 * @param [in] m Reference to mesh instance
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
real_t getMinSize(const Mesh& m);

/** \fn real_t getMinElementMeasure(const Mesh& m)
 *  \ingroup Mesh
 *  \brief Return minimal measure (length, area or volume) of elements of given mesh
 *  @param [in] m Reference to mesh instance
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
real_t getMinElementMeasure(const Mesh& m);

/** \fn real_t getMaxElementMeasure(const Mesh& m)
 *  \ingroup Mesh
 *  \brief Return maximal measure (length, area or volume) of elements of given mesh
 *  @param [in] m Reference to mesh instance
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
real_t getMaxElementMeasure(const Mesh &m);

/** \fn real_t getMinSideMeasure(const Mesh& m)
 * \ingroup Mesh
 * \brief Return minimal measure (length or area) of sides of given mesh
 * @param [in] m Reference to mesh instance
 * @note Use this function only if sides are present in the mesh and for 2-D meshes
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
real_t getMinSideMeasure(const Mesh &m);

/** \fn real_t getMaxSideMeasure(const Mesh& m)
 * \ingroup Mesh
 * \brief Return maximal measure (length or area) of sides of given mesh
 * @param [in] m Reference to mesh instance
 * @note Use this function only if sides are present in the mesh and for 2-D meshes
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
real_t getMaxSideMeasure(const Mesh &m);

/** \fn real_t getMeanElementMeasure(const Mesh& m)
 * \ingroup Mesh
 * \brief Return average measure (length, area or volume) of elements of given mesh
 * @param [in] m Reference to mesh instance
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
real_t getMeanElementMeasure(const Mesh &m);

/** \fn real_t getMeanSideMeasure(const Mesh& m)
 * \ingroup Mesh
 * \brief Return average measure (length or area) of sides of given mesh
 * @param [in] m Reference to mesh instance
 * @note Use this function only if sides are present in the mesh and for 2-D meshes
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
real_t getMeanSideMeasure(const Mesh &m);

/** \fn void setNodeCodes(Mesh& m, const string& exp, int code, size_t dof=1)
 * \ingroup Mesh
 * \brief Assign a given code to all nodes satisfying a boolean expression using
 * node coordinates
 * @param [in] m Reference to mesh instance
 * @param [in] exp Regular expression using <tt>x</tt>, <tt>y</tt>, and <tt>z</tt>
 * coordinates of nodes, according to <tt>exprtk</tt> parser
 * @param [in] code Code to assign
 * @param [in] dof Degree of freedom for which code is assigned [Default: <tt>1</tt>]
 */
void setNodeCodes(Mesh&         m,
                  const string& exp,
                  int           code,
                  size_t        dof=1);

/** \fn void setBoundaryNodeCodes(Mesh& m, const string& exp, int code, size_t dof=1)
 * \ingroup Mesh
 * \brief Assign a given code to all nodes on boundary that satisfy a boolean expression using
 * node coordinates
 * @param [in] m Reference to mesh instance
 * @param [in] exp Regular expression using <tt>x</tt>, <tt>y</tt>, and <tt>z</tt> 
 * coordinates of nodes, according to <tt>exprtk</tt> parser
 * @param [in] code Code to assign
 * @param [in] dof Degree of freedom for which code is assigned [Default: <tt>1</tt>]
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
void setBoundaryNodeCodes(Mesh&         m,
                          const string& exp,
                          int           code,
                          size_t        dof=1);

/** \fn int NodeInElement(const Node *nd, const Element *el)
 * \ingroup Mesh
 * \brief Say if a given node belongs to a given element
 * @param [in] nd Pointer to Node
 * @param [in] el Pointer to Element
 * @return Local label of the node if this one is found, <tt>0</tt> if not.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
int NodeInElement(const Node*    nd,
                  const Element* el);

/** \fn int NodeInSide(const Node *nd, const Side *sd)
 * \ingroup Mesh
 * \brief Say if a given node belongs to a given side
 * @param [in] nd Pointer to Node
 * @param [in] sd Pointer to Side
 * @return Local label of the node if this one is found, <tt>0</tt> if not.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
int NodeInSide(const Node* nd,
               const Side* sd);

/** \fn int SideInElement(const Side *sd, const Element *el)
 * \ingroup Mesh
 * \brief Say if a given side belongs to a given element
 * @param [in] sd Pointer to Side
 * @param [in] el Pointer to Element
 * @return Local label of the side if this one is found, <tt>0</tt> if not.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
int SideInElement(const Side*    sd,
                  const Element* el);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
