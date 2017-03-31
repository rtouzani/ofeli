/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

               Definition of class Mesh for finite element meshes

  ==============================================================================*/

#ifndef __MESH_H
#define __MESH_H

#include <vector>
using std::vector;

#include "OFELI_Config.h"
#include "util/macros.h"
#include "linear_algebra/Point.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "mesh/Edge.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Mesh.h
 *  \brief Definition file for class Mesh.
 */

/** \defgroup Global Global Variables
 *  \brief All global variables in the library
 */

/** \defgroup Mesh Finite Element Mesh
 *  \brief %Mesh management classes
 */

/// \ingroup Global
/// \brief A pointer to Node.
/// \details Useful for loops on nodes
    extern Node *theNode;

/// \ingroup Global
/// \brief A pointer to Element.
/// \details Useful for loops on elements
    extern Element *theElement;

/// \ingroup Global
/// \brief A pointer to Side.
/// \details Useful for loops on sides
    extern Side *theSide;
   
/// \ingroup Global
/// \brief A pointer to Edge.
/// \details Useful for loops on edges
    extern Edge *theEdge;

/*! \class Mesh
 *  \ingroup Mesh
 * \brief To store and manipulate finite element meshes.
 *
 * \details Class Mesh enables defining as an object a finite element mesh. A finite
 * element mesh is characterized by its nodes, elements and sides. Each of
 * these types of data constitutes a class in the OFELI library.
 *
 * The standard procedure to introduce the finite element mesh is to provide
 * an input file containing its data. For this, we have defined our own mesh
 * data file (following the XML syntax). Of course, a developer can write his own
 * function to read his finite element mesh file using the methods in Mesh.
 *
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern Element *the_element;
extern Node    *the_node;
extern Side    *the_side;
extern Edge    *the_edge;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

class Grid;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
typedef std::vector<Node *>    NodeSet;
typedef std::vector<Element *> ElementSet;
typedef std::vector<Side *>    SideSet;
typedef std::vector<Edge *>    EdgeSet;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

class Mesh
{

 public:

//----------------------------   BASIC OPERATIONS   -----------------------------

/// \brief Default constructor (Empty mesh)
    Mesh();

/** \brief Constructor using a mesh file
 *  @param [in] file File containing mesh data. The extension of the file yields the file format:
 *              The extension .m implies OFELI file format and .msh implies GMSH msh file.
 *  @param [in] bc Flag to remove (true) or not (false) imposed Degrees of Freedom [default: false]
 *  @param [in] opt Type of DOF support: To choose among enumerated values <tt>NODE_DOF</tt>, <tt>SIDE_DOF</tt> 
 *              or <tt>ELEMENT_DOF</tt>.\n
 *  Say if degrees of freedom (unknowns) are supported by nodes, sides or elements.
 *  @param [in] nb_dof Number of degrees of freedom per node [Default: <tt>1</tt>]. This value is meaningful only
 *  if other format than OFELI's one is used. Otherwise, the information is contained in the OFELI file format.
 */
    Mesh(const string& file,
         bool          bc=false,
         int           opt=NODE_DOF,
         int           nb_dof=1);

/** \brief Constructor for a 1-D mesh.
 *  The domain is the interval [0,L]
 *  @param [in] L Length of the interval
 *  @param [in] nb_el Number of elements to generate
 *  @param [in] p Degree of finite element polynomial (Default = 1)
 *  @param [in] nb_dof Number of degrees of freedom for each node (Default = 1)
 */
    Mesh(real_t L,
         size_t nb_el,
         size_t p=1,
         size_t nb_dof=1);

/** \brief Constructor for a uniform finite difference grid given by and instance of class Grid.
 *  @param [in] g Grid instance
 *  @param [in] opt Optional value to say which type of elements to generate
 *     - TRIANGLE: %Mesh elements are triangles
 *     - QUADRILATERAL: %Mesh elements are quadrilaterals [default]
 */
    Mesh(const Grid& g,
         int         opt=QUADRILATERAL);

/** \brief Constructor of dual mesh for a uniform finite difference grid given by and instance of 
 *  class Grid.
 *  @param [in] g Grid instance
 *  @param [in] shape Value to say which type of elements to generate
 *     - TRIANGLE: %Mesh elements are triangles
 *     - QUADRILATERAL: %Mesh elements are quadrilaterals [default]
 *  @param [in] opt This argument can take any value. It is here only to distinguish from the
 *  other constructor using Grid instance.
 *  @remark This constructor is to be used to obtain a dual mesh from a structured grid. It is
 *  mainly useful if a cell centered finite volume method is used.
 */
    Mesh(const Grid& g,
         int         shape,
         int         opt);

/** \brief Constructor for a uniform 1-D finite element mesh.
 *  \details The domain is the line (xmin,xmax)
 *  @param [in] xmin Minimal coordinate
 *  @param [in] xmax Maximal coordinate
 *  @param [in] ne Number of elements
 *  @param [in] c1 Code for the first node (x=xmin)
 *  @param [in] c2 Code for the last node (x=xmax)
 *  @param [in] opt Flag to generate elements as well (if not zero) [Default: 0].
 *  @remark The option opt can be set to 0 if the user intends to use finite differences.
 */
    Mesh(real_t xmin,
         real_t xmax,
         size_t ne,
         int    c1,
         int    c2,
         int    opt=0);

/** \brief Constructor for a uniform 2-D structured finite element mesh.
 *  \details The domain is the rectangle (xmin,xmax)x(ymin,ymax)
 *  @param [in] xmin Minimal x-coordinate
 *  @param [in] xmax Maximal x-coordinate
 *  @param [in] ymin Minimal y-coordinate
 *  @param [in] ymax Maximal y-coordinate
 *  @param [in] nx Number of subintervals on the x-axis
 *  @param [in] ny Number of subintervals on the y-axis
 *  @param [in] cx0 Code for nodes generated on the line x=x0 if >0, for sides on this line if <0
 *  @param [in] cxN Code for nodes generated on the line x=xN if >0, for sides on this line if <0
 *  @param [in] cy0 Code for nodes generated on the line y=y0 if >0, for sides on this line if <0
 *  @param [in] cyN Code for nodes generated on the line y=yN if >0, for sides on this line if <0
 *  @param [in] opt Flag to generate elements as well (if not zero) [Default: 0]. 
 *  If the flag is not 0, it can take one of the enumerated values: TRIANGLE or QUADRILATERAL, 
 *  with obvious meaning.
 *  @remark The option opt can be set to 0 if the user intends to use finite differences.
 */
    Mesh(real_t xmin,
         real_t xmax,
         real_t ymin,
         real_t ymax,
         size_t nx,
         size_t ny,
         int    cx0,
         int    cxN,
         int    cy0,
         int    cyN,
         int    opt=0);

/** \brief Constructor for a uniform 3-D structured finite element mesh.
 *  \details The domain is the parallepiped (xmin,xmax)x(ymin,ymax)x(zmin,zmax)
 *  @param [in] xmin Minimal x-coordinate
 *  @param [in] xmax Maximal x-coordinate
 *  @param [in] ymin Minimal y-coordinate
 *  @param [in] ymax Maximal y-coordinate
 *  @param [in] zmin Minimal z-coordinate
 *  @param [in] zmax Maximal z-coordinate
 *  @param [in] nx Number of subintervals on the x-axis
 *  @param [in] ny Number of subintervals on the y-axis
 *  @param [in] nz Number of subintervals on the z-axis
 *  @param [in] cx0 Code for nodes generated on the line x=xmin if >0, for sides on this line if <0
 *  @param [in] cxN Code for nodes generated on the line x=xmax if >0, for sides on this line if <0
 *  @param [in] cy0 Code for nodes generated on the line y=ymin if >0, for sides on this line if <0
 *  @param [in] cyN Code for nodes generated on the line y=ymax if >0, for sides on this line if <0
 *  @param [in] cz0 Code for nodes generated on the line z=zmin if >0, for sides on this line if <0
 *  @param [in] czN Code for nodes generated on the line z=zmax if >0, for sides on this line if <0
 *  @param [in] opt Flag to generate elements as well (if not zero) [Default: 0]. 
 *  If the flag is not 0, it can take one of the enumerated values: HEXAHEDRON or TETRAHEDRON, 
 *  with obvious meaning.
 *  @remark The option opt can be set to 0 if the user intends to use finite differences.
 */
    Mesh(real_t xmin,
         real_t xmax,
         real_t ymin,
         real_t ymax,
         real_t zmin,
         real_t zmax,
         size_t nx,
         size_t ny,
         size_t nz,
         int    cx0,
         int    cxN,
         int    cy0,
         int    cyN,
         int    cz0,
         int    czN,
         int    opt);

/** \brief Constructor that extracts the mesh of a rectangular region from an initial mesh.
 *  \details This constructor is useful for zooming purposes for instance.
 *  @param [in] m Initial mesh from which the submesh is extracted
 *  @param [in] x_bl Coordinate of bottom left vertex of the rectangle
 *  @param [in] x_tr Coordinate of top right vertex of the rectangle
 */
    Mesh(const Mesh&          m,
         const Point<real_t>& x_bl,
         const Point<real_t>& x_tr);

/** \brief Constructor that copies the input mesh and selects given degrees of freedom.
 *  \details This constructor is to be used for coupled problems where each
 *  subproblem uses a choice of degrees of freedom.
 *  @param [in] mesh Initial mesh from which the submesh is extracted
    @param [in] opt Type of DOF support: To choose among enumerated values <tt>NODE_DOF</tt>, <tt>SIDE_DOF</tt> or <tt>ELEMENT_DOF</tt>.
 *  @param [in] dof1 Label of first degree of freedom to select to the output mesh
 *  @param [in] dof2 Label of last degree of freedom to select to the output mesh
    @param [in] bc Flag to remove (<tt>true</tt>) or not (<tt>false</tt>) imposed Degrees of Freedom [Default: <tt>false</tt>]
 */
    Mesh(const Mesh& mesh,
         int         opt,
         size_t      dof1,
         size_t      dof2,
         bool        bc=false);

/// \brief Copy Constructor
/// @param [in] ms Mesh instance to copy
    Mesh(const Mesh& ms);

/// Destructor
    ~Mesh();

//----------------------------   MODIFIERS   --------------------------------

/** \brief Define space dimension.
 *  Normally, between 1 and 3.
 *  @param [in] dim Space dimension to set (must be between 1 and 3)
 */
    void setDim(size_t dim) { _dim = dim; }

/** \brief Define Verbose Parameter.
 *  Controls output details
 *  @param [in] verb verbosity parameter (Must be between 0 and 10)
 */
    void setVerbose(int verb) { _verb = verb; }

/// \brief Add a node to mesh.
/// @param [in] nd Pointer to Node to add
    void Add(Node* nd);

/// \brief Add an element to mesh.
/// @param [in] el Pointer to Element to add
    void Add(Element* el);

/// \brief Add a side to mesh.
/// @param [in] sd Pointer to Side to add
    void Add(Side* sd);

/// \brief Add an edge to mesh.
/// @param [in] ed Pointer to Edge to add
    void Add(Edge* ed);
    
/// \brief Operator <tt>*=</tt>
/// \details Rescale mesh coordinates by myltiplying by a factor
/// @param [in] a Value to multiply by
    Mesh &operator*=(real_t a);

/** \brief Read mesh data in file
 *  \details Mesh file must be in <tt>OFELI</tt> format. See "File Formats" page
 *  @param [in] mesh_file Mesh file name
 */
    void get(const string& mesh_file);

/** \brief Read mesh data in file with giving its format
 *  \details File format can be chosen among a variety of choices. See "File Formats" page
 *  @param [in] mesh_file Mesh file name
 *  @param [in] ff File format: Integer to chose among enumerated values:
 *  <tt>OFELI_FF</tt>, <tt>GMSH</tt>, <tt>MATLAB</tt>, <tt>EASYMESH</tt>, <tt>GAMBIT</tt>, <tt>BAMG</tt>, 
 *  <tt>NETGEN</tt>, <tt>TRIANGLE_FF</tt>
 *  @param [in] nb_dof Number of degrees of freedom per node (Default value: 1)
 */
    void get(const string& mesh_file,
             int           ff,
             int           nb_dof=1);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void Get(const string& mesh_file) { get(mesh_file); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Define supports of degrees of freedom
 *  @param [in] opt DOF type:
 *  <ul>
 *    <li><tt>NODE_DOF</tt>: Degrees of freedom are supported by nodes
 *    <li><tt>SIDE_DOF</tt>: Degrees of freedom are supported by sides
 *    <li><tt>EDGE_DOF</tt>: Degrees of freedom are supported by edges
 *    <li><tt>ELEMENT_DOF</tt>: Degrees of freedom are supported by elements
 *  </ul>
 *  @param [in] nb_nodes Number of nodes on sides or elements (default=1).
 *  This parameter is useful only if dofs are supported by sides or elements
 *  @note This member function creates all mesh sides if the option <tt>ELEMENT_DOF</tt> or
 *  <tt>SIDE_DOF</tt> is selected. So it not necessary to call getAllSides() after
 */
    void setDOFSupport(int opt,
                       int nb_nodes=1);

/** \brief Define number of degrees of freedom for each node
 *  @param [in] nb_dof Number of degrees of freedom (unknowns) for each mesh node
 *  (Default value is <tt>1</tt>)
 *  @note This function first declares nodes as unknown supports, sets the
 *  number of degrees of freedom and renumbers equations
 */
    void setNbDOFPerNode(size_t nb_dof=1);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Set nodes as supports for degrees of freedom
    void setNodesForDOF() { setDOFSupport(NODE_DOF); }

/** \brief Set sides as supports for degrees of freedom
 *  @note This member function creates all mesh sides. So it not necessary
 *  to call getAllSides() after
 */
    void setSidesForDOF() { setDOFSupport(SIDE_DOF); }

/// \brief Set edges as supports for degrees of freedom
    void setEdgesForDOF() { setDOFSupport(EDGE_DOF); }

/** \brief Set elements as supports for degrees of freedom
 *  @note This member function creates all mesh sides. So it not necessary
 *  to call getAllSides() after
 */
    void setElementsForDOF() { setDOFSupport(ELEMENT_DOF); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Define a point in the domain.
 *  This function makes sense only if boundary mesh is given without
 *  internal mesh (Case of Boundary Elements)
 *  @param [in] x Coordinates of point to define
 */
    void setPointInDomain(Point<real_t> x);

/// \brief Eliminate equations corresponding to imposed DOF
    void removeImposedDOF() { _no_imposed_dof = true; NumberEquations(); }

/** \brief Renumber Equations
 *  @param [in] dof Label of degree of freedom for which numbering is performed.
 *  Default value (<tt>0</tt>) means that all degrees of freedom are taken into account
 */
    size_t NumberEquations(size_t dof=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Select degrees of freedom and renumber corresponding equations
 *  \details This function is to be used when degrees of freedom are grouped
 *  for different systems of equations
 *  @param [in] dof_type Type of support of dof. To be given among the enumerated
 *  values: <tt>NODE_DOF</tt>, <tt>ELEMENT_DOF</tt>, <tt>SIDE_DOF</tt> 
 *  @param [in] dof1 Label of first degree of freedom for which numbering is performed.
 *  @param [in] dof2 Label of second degree of freedom for which numbering is performed.
 *  @param [in] bc Flag to remove (<tt>true</tt>) or not (<tt>false</tt>) imposed Degrees of Freedom
 *  [default: <tt>false</tt>]
 */
    void selectDOF(int    dof_type,
                   size_t dof1,
                   size_t dof2,
                   bool   bc=false);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Renumber Equations
 *  @param [in] dof Label of degree of freedom for which numbering is performed.
 *  @param [in] c code for which degrees of freedom are enforced.
 */
    size_t NumberEquations(size_t dof,
                           int    c);

/// \brief Determine all mesh sides.
/// @return Number of all sides.
    int getAllSides(int opt=0);

/// \brief Return the number of nodes on each side.
    size_t getNbSideNodes() const { return _nb_side_nodes; }

/// \brief Return the number of nodes in each element.
    size_t getNbElementNodes() const { return _nb_element_nodes; }

/// \brief Determine all boundary sides
/// @return Number of boundary sides.
    int getBoundarySides();

/** \brief Create list of boundary sides.
 *  \details This function is useful to loop over boundary sides without testing
 *  Once this one is called, the function getNbBoundarySides() is available.
 *  Moreover, looping over boundary sides is available via the member
 *  functions topBoundarySide() and getBoundarySide()
 *  @return Number of boundary sides.
 */
    int createBoundarySideList();

/// \brief Determine all boundary nodes
///  @return n Number of boundary nodes.
    int getBoundaryNodes();

/** \brief Create list of internal sides (not on the boundary).
 *  \details This function is useful to loop over internal sides without testing
 *  Once this one is called, the function getNbInternalSides() is available.
 *  Moreover, looping over internal sides is available via the member
 *  functions topInternalSide() and getInternalSide()
 *  @return n Number of internal sides.
 */
    int createInternalSideList();

/// \brief Determine all edges
/// @return Number of all edges.
    int getAllEdges();

/** \brief Create node neighboring elements
 *  \details This function is generally useful when, for a numerical method, one looks
 *  for a given node to the list of elements that share this node.
 *  Once this function is invoked, one can retrieve the list of neighboring elements of
 *  any node (Node::getNeigEl)
 */
    void getNodeNeighborElements();

/** \brief Create element neighboring elements
 *  \details This function creates for each element the list of elements that share
 *  a side with it.
 *  Once this function is invoked, one can retrieve the list of neighboring elements of
 *  any element (Element::getNeigborElement)
 */
    void getElementNeighborElements();

/** \brief Associate material to code of element
 *  @param [in] code Element code for which material is assigned
 *  @param [in] mname Name of material
 */
    void setMaterial(int           code,
                     const string& mname);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Imbed mesh into a coarse mesh ms.
/// @param [in] ms Instance to coarse mesh to imbed in
/// @param [in] test_el
    void inCoarse(Mesh& ms,
                  bool  test_el=false);

/// \brief Imbed mesh into fine mesh ms.
/// @param [in] ms Instance to fine mesh to imbed in
/// @param [in] test_el
    void inFine(Mesh& ms,
                bool  test_el=false);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Renumber mesh nodes according to reverse Cuthill Mc Kee algorithm
/// @param [in] m Memory size needed for matrix graph (default value is GRAPH_MEMORY, see OFELI_Config.h)
    void Reorder(size_t m=GRAPH_MEMORY) { RenumberNodes(m); }

/** \brief Add a node by giving its label and an array containing its coordinates
 *  @param [in] num Label of node to add
 *  @param [in] x C-array of node coordinates
 */
    void Add(size_t  num,
             real_t* x);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Reset mesh by removing all inactive nodes, elements and sides
/// @remark: All information on adaptivity is lost.
    void Reset();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Set discontinuous approximations
    void setDiscontinuous(size_t p);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Remove a node given by its label.
 *  \details This function does not release the space previously occupied
 *  @param [in] label Label of node to delete
 */
    void DeleteNode(size_t label) { Delete(getPtrNode(label)); }

/** \brief Remove an element given by its label.
 *  \details This function does not release the space previously occupied
 *  @param [in] label Label of element to delete
 */
    void DeleteElement(size_t label) { Delete(getPtrElement(label)); }

/** \brief Remove a side given by its label.
 *  \details This function does not release the space previously occupied
 *  @param [in] label Label of side to delete
 */
    void DeleteSide(size_t label) { Delete(getPtrSide(label)); }

/** \brief Remove a node given by its pointer.
 *  \details This function does not release the space previously occupied
 *  @param [in] nd Pointer to node to delete
 */
    void Delete(Node* nd);

/** \brief Remove a node given by its pointer.
 *  \details This function does not release the space previously occupied
 *  @param [in] el Pointer to element to delete
 */
    void Delete(Element* el);

/** \brief Remove a side given by its pointer.
 *  \details This function does not release the space previously occupied
 *  @param [in] sd Pointer to side to delete
 */
    void Delete(Side* sd);

/** \brief Remove an edge given by its pointer.
 *  \details This function does not release the space previously occupied
 *  @param [in] ed Pointer to edge to delete
 */
    void Delete(Edge* ed);

/** \brief Renumber a node
 *  @param [in] n1 Old label
 *  @param [in] n2 New label
 */
    void RenumberNode(size_t n1,
                      size_t n2);

/** \brief Renumber an element
 *  @param [in] n1 Old label
 *  @param [in] n2 New label
 */
    void RenumberElement(size_t n1,
                         size_t n2);

/** \brief Renumber a side
 *  @param [in] n1 Old label
 *  @param [in] n2 New label
 */
    void RenumberSide(size_t n1,
                      size_t n2);

/** \brief Renumber an edge
 *  @param [in] n1 Old label
 *  @param [in] n2 New label
 */
    void RenumberEdge(size_t n1,
                      size_t n2);

/** \brief Set viewing window for nodes
 *  @param [in] n1 First node to view
 *  @param [in] n2 last node to view
 */
    void setNodeView(size_t n1,
                     size_t n2);

/** \brief Set viewing window for elements
 *  @param [in] n1 First element to view
 *  @param [in] n2 last element to view
 */
    void setElementView(size_t n1,
                        size_t n2);

/** \brief Set viewing window for sides
 *  @param [in] n1 First side to view
 *  @param [in] n2 last side to view
 */
    void setSideView(size_t n1,
                     size_t n2);
   
/** \brief Set viewing window for edges
 *  @param [in] n1 First edge to view
 *  @param [in] n2 last edge to view
 */
    void setEdgeView(size_t n1,
                     size_t n2);

/// \brief Initialize list of mesh nodes using the input vector
/// @param [in] nl vector instance that contains the list of pointers to nodes
    void setList(const std::vector<Node *>& nl);

/// \brief Initialize list of mesh elements using the input vector
/// @param [in] el vector instance that contains the list of pointers to elements
    void setList(const std::vector<Element *>& el);

/// \brief Initialize list of mesh sides using the input vector
/// @param [in] sl vector instance that contains the list of pointers to sides
    void setList(const std::vector<Side *>& sl);

/** \brief Rescale mesh by multiplying node coordinates by constants
 *  \details This function can be used e.g. for changing coordinate units
 *  @param [in] sx Factor to multiply by <tt>x</tt> coordinates
 *  @param [in] sy Factor to multiply by <tt>y</tt> coordinates [Default: <tt>sx</tt>]
 *  @param [in] sz Factor to multiply by <tt>z</tt> coordinates [Default: <tt>sx</tt>]
 */
    void Rescale(real_t sx, real_t sy=0., real_t sz=0.);

//-----------------------------   INSPECTORS  ----------------------------------

/// \brief Return Verbose Parameter
    int getVerbose() const { return _verb; }

/// \brief Return space dimension
    size_t getDim() const { return _dim; }

/// \brief Return number of nodes
    size_t getNbNodes() const { return _nb_nodes; }

/// \brief Return number of marked nodes
    size_t getNbMarkedNodes() const { return _nb_marked_nodes; }

/// \brief Return number of vertices
    size_t getNbVertices() const { return _nb_vertices; }

/// \brief Return total number of degrees of freedom (DOF)
    size_t getNbDOF() const { return _nb_dof; }

/// \brief Return number of equations
    size_t getNbEq() const { return _nb_eq; }

/// \brief Return number of equations for the i-th set of degrees of freedom
    size_t getNbEq(int i) const { return _dof_nbeq[i-1]; }

/// \brief Return number of elements
    size_t getNbElements() const { return _nb_elements; }

/// \brief Return number of sides
    size_t getNbSides() const { return _nb_sides; }

/// \brief Return number of sides
    size_t getNbEdges() const { return _nb_edges; }

/// \brief Return number of boundary sides.
/// \details This function is valid if member function \b getAllSides or \b getBoundarySides
/// has been invoked before
    size_t getNbBoundarySides() const { return _nb_boundary_sides; }

/// \brief Return number of internal sides.
/// \details This function is valid if member functions \b getAllSides and \b createInternalSideList
/// have been invoked before
    size_t getNbInternalSides() const { return _nb_internal_sides; }

/// \brief Return number of materials
    size_t getNbMat() const { return _nb_mat; }

/// \brief Add mid-side nodes
/// \details This is function is valid for triangles only
/// @param [in] g Option to say of barycentre node is to be added (>0) or not (=0)
    void AddMidNodes(int g=0);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void AddNodes(int p=2);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return maximum coordinates of nodes
    Point<real_t> getMaxCoord() const;

/// \brief Return minimum coordinates of nodes
    Point<real_t> getMinCoord() const;

/** \brief Replace node in the mesh
 *  \details If the node label exists already, the existing node pointer will
 *  be replaced by the current one. If not, an error message is displayed.
 *  @param [in] nd Pointer to node
 */
    void set(Node* nd);

/** \brief Replace element in the mesh
 *  \details If the element label exists already, the existing element pointer will
 *  be replaced by the current one. If not, an error message is displayed.
 *  @param [in] el Pointer to element
 */
    void set(Element* el);

/** \brief Choose side in the mesh
 *  \details If the side label exists already, the existing side pointer will
 *  be replaced by the current one. If not, an error message is displayed.
 *  @param [in] sd Pointer to side
 */
    void set(Side* sd);

/// \brief Return information about DOF type.
/// @return true if DOF are supported by nodes, <tt>false</tt> otherwise
    bool NodesAreDOF() const { return _set_nodes; }

/// \brief Return information about DOF type.
/// @return true if DOF are supported by sides, <tt>false</tt> otherwise
    bool SidesAreDOF() const { return _set_sides; }

/// \brief Return information about DOF type.
/// @return true if DOF are supported by edges, <tt>false</tt> otherwise
    bool EdgesAreDOF() const { return _set_edges; }

/// \brief Return information about DOF type.
/// @return true if DOF are supported by elements, <tt>false</tt> otherwise
    bool ElementsAreDOF() const { return _set_elements; }

/** \brief Return information on dof support
 *  Return an integer according to enumerated values: NODE_DOF, ELEMENT_DOF
 *  SIDE_DOF
 */
    int getDOFSupport() const;
 
#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Return pointer to an element (in an imbedded mesh) to which node <tt>*nd</tt> belongs.
/// This information is valid if utility function has been used.
    Element* InCoarse(const Node* nd) const { return _node_in_coarse_element[nd->n()-1]; }

/// \brief Return pointer to an element (in an imbedded mesh) to which node <tt>*nd</tt> belongs.
/// This information is valid if utility function has been used.
    Element* InFine(const Node* nd) const { return _node_in_fine_element[nd->n()-1]; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Write mesh data on file
/// @param [in] mesh_file %Mesh file name
    void put(const string& mesh_file) const;

/** \brief Write mesh data on file in various formats
 *  \details File format depends on the extension in file name
 * @param [in] mesh_file %Mesh file name
 * If the extension is '.m', the output file is an OFELI file
 * If the extension is '.gpl', the output file is a Gnuplot file
 * If the extension is '.msh' or '.geo', the output file is a Gmsh file
 * If the extension is '.vtk', the output file is a VTK file
 */
    void save(const string& mesh_file) const;

/// \brief Return true if imposed DOF count in equations, false if not
    bool withImposedDOF() const { return !_no_imposed_dof; }

/// \brief Return true is mesh is structured, false if not
    bool isStructured() const { return _is_structured; }

/// \brief Return new label of node of a renumbered node
    size_t getNodeNewLabel(size_t n) const { return _node_new_label[n-1]; }

/// \brief Fill vector <tt>nl</tt> with list of pointers to nodes
/// @param [out] nl Instance of class vector that contain on output the list
    void getList(vector<Node *>& nl) const;

/// \brief Fill vector <tt>el</tt> with list of pointers to elements
/// @param [out] el Instance of class vector that contain on output the list
    void getList(vector<Element *>& el) const;

/// \brief Fill vector <tt>sl</tt> with list of pointers to sides
/// @param [out] sl Instance of class vector that contain on output the list
    void getList(vector<Side *>& sl) const;

//---------------------------------  ITERATORS  --------------------------------

/// \brief Return pointer to node with label <tt>i</tt>.
    Node *getPtrNode(size_t i) const { return _nodes[i-1]; }

/// \brief Return refenrece to node with label <tt>i</tt>
    Node &getNode(size_t i) const { return *(_nodes[i-1]); }

/// \brief Return pointer to element with label <tt>i</tt>
    Element *getPtrElement(size_t i) const { return _elements[i-1]; }

/// \brief Return reference to element with label <tt>i</tt>
    Element &getElement(size_t i) const { return *(_elements[i-1]); }

/// \brief Return pointer to side with label <tt>i</tt>
    Side *getPtrSide(size_t i) const { return _sides[i-1]; }

/// \brief Return reference to side with label <tt>i</tt>
    Side &getSide(size_t i) const { return *(_sides[i-1]); }

/// \brief Return pointer to edge with label <tt>i</tt>
    Edge *getPtrEdge(size_t i) const { return _edges[i-1]; }

/// \brief Return reference to edge with label <tt>i</tt>
    Edge &getEdge(size_t i) const { return *(_edges[i-1]); }

/// \brief Return label of <tt>i</tt>-th node
/// @param [in] i Node index
    size_t getNodeLabel(size_t i) const { return _nodes[i-1]->n(); }

/// \brief Return label of <tt>i</tt>-th element
/// @param [in] i Element index
    size_t getElementLabel(size_t i) const { return _elements[i-1]->n(); }

/// \brief Return label of <tt>i</tt>-th side
/// @param [in] i Side index
    size_t getSideLabel(size_t i) const { return _sides[i-1]->n(); }

/// \brief Return label of <tt>i</tt>-th edge
/// @param [in] i Edge index
    size_t getEdgeLabel(size_t i) const { return _edges[i-1]->n(); }

/// \brief Reset list of nodes at its top position (Non constant version)
    void topNode() const { _node_it = 0; }

/// \brief Reset list of boundary nodes at its top position (Non constant version)
    void topBoundaryNode() const { _node_it = 0; }

/// \brief Reset list of marked nodes at its top position (Non constant version)
    void topMarkedNode() const { _node_it = 0; }

/// \brief Reset list of elements at its top position (Non constant version)
    void topElement() const { _element_it = 0; }

/// \brief Reset list of sides at its top position (Non constant version)
    void topSide() const { _side_it = 0; }

/// \brief Reset list of boundary sides at its top position (Non constant version)
    void topBoundarySide() const { _side_it = 0; }

/// \brief Reset list of intrenal sides at its top position (Non constant version)
    void topInternalSide() const { _side_it = 0; }

/// \brief Reset list of edges at its top position (Non constant version)
    void topEdge() const { _edge_it = 0; }

/// \brief Reset list of boundary edges at its top position (Non constant version)
    void topBoundaryEdge() const { _edge_it = 0; }

/// \brief Return pointer to current node and move to next one (Non constant version)
    Node *getNode() const
    {
       if (_node_it==_nb_nodes)
          return NULL;
       else
          return _nodes[_node_it++];
    }

/// \brief Return pointer to current boundary node and move to next one (Non constant version)
    Node* getBoundaryNode() const
    {
       if (_node_it==_nb_boundary_nodes)
          return NULL;
       else
          return _boundary_nodes[_node_it++];
    }

/// \brief Return pointer to current marked node and move to next one (Non constant version)
    Node* getMarkedNode() const
    {
       if (_node_it==_nb_marked_nodes)
          return NULL;
       else
          return _marked_nodes[_node_it++];
    }

/// \brief Return pointer to current element and move to next one (Non constant version)
    Element* getElement() const
    {
       if (_element_it==_nb_elements)
          return NULL;
       else
          return _elements[_element_it++];
    }

/** \brief Return pointer to current element and move to next one (Non constant version)
 * \details This function returns pointer to the current element only is this one is
 * active. Otherwise it goes to the next active element (To be used when adaptive meshing is involved)
 */
    Element* getActiveElement() const
    {
       if (_element_it==_nb_elements)
          return NULL;
       else {
          while (_elements[_element_it]->isActive()==false)
             _element_it++;
          return _elements[_element_it++];
       }
    }

/// \brief Return pointer to current side and move to next one (Non constant version)
    Side* getSide() const
    { 
       if (_side_it==_nb_sides)
          return NULL;
       else
          return _sides[_side_it++];
    }

/// \brief Return pointer to current boundary side and move to next one (Non constant version)
    Side* getBoundarySide() const
    {
       if (_side_it==_nb_boundary_sides)
          return NULL;
       else
          return _boundary_sides[_side_it++];
    }

/// \brief Return pointer to current internal side and move to next one (Non constant version)
    Side* getInternalSide() const
    {
       if (_side_it==_nb_internal_sides)
          return NULL;
       else
          return _internal_sides[_side_it++];
    }

/// \brief Return pointer to current edge and move to next one (Non constant version)
    Edge* getEdge() const
    {
       if (_edge_it==_nb_edges)
          return NULL;
       else
          return _edges[_edge_it++];
    }

/// \brief Return pointer to current boundary edge and move to next one (Non constant version)
    Edge* getBoundaryEdge() const
    {
       if (_edge_it==_nb_boundary_edges)
          return NULL;
       else
          return _boundary_edges[_edge_it++];
    }

/// \brief Determine shape of elements
/// Return Shape index (see enum ElementShape) if all elements have the same shape, 0 if not.
    int getShape() const;
 
//---------------------------------  OPERATORS  --------------------------------

/// \brief Operator () : Return pointer to i-th element
    Element* operator()(size_t i) const { return getPtrElement(i); }

/// \brief Operator [] : Return pointer to i-th node
    Node* operator[](size_t i) const { return getPtrNode(i); }

/// \brief Operator () : Return pointer to i-th node of n-th element
    size_t operator()(size_t i, size_t n) const { return getPtrElement(n)->getNodeLabel(i); }

/// \brief Operator = : Assign a Mesh instance
    Mesh& operator=(Mesh& ms);

/** \brief Refine mesh.
 *  Subdivide each triangle into 4 subtriangles. This member function is valid for
 *  2-D triangular meshes only.
 *  @param [in] in_mesh Input mesh
 *  @param [out] out_mesh Output mesh
 */
    friend void Refine(Mesh& in_mesh,
                       Mesh& out_mesh);

    friend class XMLParser;
    friend ostream& operator<<(ostream& s, const Mesh& ms);

 private:
   size_t            _nb_nodes, _nb_elements, _nb_sides, _nb_edges, _nb_side_nodes, _nb_element_nodes;
   size_t            _dim, _nb_dof, _nb_vertices, _first_dof, _nb_mat;
   size_t            _nb_boundary_nodes, _nb_boundary_sides, _nb_internal_sides, _nb_boundary_edges;
   size_t            _nb_marked_nodes, _nb_eq;
   size_t            _max_nb_nodes, _max_nb_elements, _max_nb_sides, _max_nb_edges;
   mutable size_t    _node_it, _side_it, _element_it, _edge_it;
   vector<size_t>    _node_old_label, _node_new_label, _element_old_label;
   mutable size_t    _const_node_it, _const_side_it, _const_element_it, _const_edge_it;
   int               _verb;
   ElementSet        _elements;
   SideSet           _sides, _boundary_sides, _internal_sides;
   EdgeSet           _edges, _boundary_edges;
   NodeSet           _nodes, _boundary_nodes, _marked_nodes;
   int               _code[MAX_NBDOF_NODE], _dof_nbeq[MAX_NBDOF_NODE];
   vector<Element *> _node_in_coarse_element, _node_in_fine_element;
   bool              _set_nodes, _set_sides, _set_elements, _set_edges;
   bool              _no_imposed_dof, _is_structured;
   bool              _all_sides_created, _boundary_sides_created, _all_edges_created, _boundary_edges_created;
   bool              _boundary_nodes_created;
   bool              _node_neighbor_elements_created, _element_neighbor_elements_created;
   int               _code_mat[MAX_NB_MATERIALS];
   string            _mat[MAX_NB_MATERIALS];
   unsigned long     _available_memory;
   size_t            _n_view1, _n_view2, _e_view1, _e_view2, _s_view1, _s_view2, _ed_view1, _ed_view2;

   void FindSideNeighbors();
   void RenumberNodes(size_t m=GRAPH_MEMORY);
   void RenumberNodes(vector<size_t> &perm, size_t m=GRAPH_MEMORY);
   unsigned long FindGraph(vector<long> &xadj, vector<size_t> &adjncy);
   void GenRCM(vector<long> &xadj, vector<size_t> &adjncy, size_t *perm,
               vector<size_t> &mask, vector<size_t> &xls);
   void checkNodeLabels();
   void checkElementLabels();
   void checkSideLabels();
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// The MVC 6 linker doesn't want here a template function: Don't know why !!!
inline bool _compare_n(Node const* a, Node const* b) { return (a->n() < b->n()); }
inline bool _compare_e(Element const* a, Element const* b) { return (a->n() < b->n()); }
inline bool _compare_s(Side const* a, Side const* b) { return (a->n() < b->n()); }

inline bool _node_compare(Node* const &a, Node* const &b) { return (a->n() < b->n()); }
inline bool _element_compare(Element* const &a, Element* const &b) { return (a->n() < b->n()); }
inline bool _side_compare(Side* const &a, Side* const &b) { return (a->n() < b->n()); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \fn ostream & operator<<(ostream& s, const Mesh &ms)
 *  \brief Output mesh data.
 *  \ingroup Mesh
 */
    ostream& operator<<(ostream&    s,
                        const Mesh& ms);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
