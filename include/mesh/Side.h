/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: yo!u can redistribute it and/or modify
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

                            Definition of class 'Side'

  ==============================================================================*/


#ifndef __SIDE_H
#define __SIDE_H

#include <vector>
using std::vector;

#include <string>
using std::string;

#include "OFELI_Config.h"
#include "io/output.h"
#include "mesh/Node.h"
#include "mesh/Edge.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Side.h
 *  \brief Definition file for class Side.
 */


/*! \enum BCType
 * To select special boundary conditions.
 */
enum BCType {
  PERIODIC_A     =  9999,     /*!< Periodic Boundary condition  (first side)      */
  PERIODIC_B     = -9999,     /*!< Periodic Boundary condition  (second side)     */
  CONTACT_BC     =  9998,     /*!< Contact Boundary conditions                    */
  CONTACT_M      =  9997,     /*!< Contact Boundary condition, set as master side */
  CONTACT_S      = -9997,     /*!< Contact Boundary condition, set as slave side  */
  SLIP           =  9996      /*!< Slip Boundary condition                        */
};

class Element;

/*! \class Side
 *  @ingroup Mesh
 * \brief To store and treat finite element sides (edges in 2-D or faces in 3-D)
 *
 * \details Defines a side of a finite element mesh. The sides are given in particular
 * by their shapes and a list of nodes. Each node can be accessed by the member
 * function \b getPtrNode(). The string defining the element shape must be chosen
 * according to the following list:
 *
 *  Shape            Shape name     Dimension    Min. number of nodes
 *  Line             line             3                 2
 *  Triangle         tria             3                 3
 *  Quadrilateral    quad             3                 4
 * 
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class Side
{

 public:

/*! \enum SideType
 * To select side type (boundary side or not).
 */
enum SideType {
   INTERNAL_SIDE     =  0,  /*!< Internal side */
   EXTERNAL_BOUNDARY =  1,  /*!< Side on external boundary */
   INTERNAL_BOUNDARY =  2   /*!< Side on internal boundary */
};

//----------------------------   BASIC OPERATIONS   -----------------------------

/// \brief Default Constructor
    Side();

/** \brief Constructor initializing side label and shape.
 *  @param [in] label Label to assign to side.
 *  @param [in] shape Shape of side (See class description).
 */
    Side(size_t        label,
         const string& shape);

/** \brief Constructor initializing side label and shape.
 *  @param [in] label to assign to side.
 *  @param [in] shape of side (See enum ElementShape in Mesh).
 */
    Side(size_t label,
         int    shape);

/// \brief Copy constructor
    Side(const Side& sd);

/// \brief Destructor
    ~Side() { }

//-------------------------------   MODIFIERS  ---------------------------------

/// \brief Insert a node at end of list of nodes of side
    void Add(Node* node);

/// \brief Insert an edge at end of list of edges of side
    void Add(Edge* edge);

/// \brief Define label of side.
    void setLabel(size_t i) { _label = i; }

/// \brief Define First DOF
    void setFirstDOF(size_t n) { _first_dof = n; }

/// \brief Set number of degrees of freedom (DOF).
    void setNbDOF(size_t nb_dof) { _nb_dof = nb_dof; }

/** \brief Define label of DOF.
 *  @param [in] i DOF index
 *  @param [in] dof Its label
 */
    void DOF(size_t i,
             size_t dof) { _dof[i-1] = dof; }

/** \brief Define number of DOF.
 *  @param [in,out] first_dof Label of the first DOF in input that is actualized
 *  @param [in] nb_dof Number of DOF
 */
    void setDOF(size_t& first_dof,
                size_t  nb_dof);

/** \brief Assign code to a DOF
 *  @param [in] dof DOF to which code is assigned
 *  @param [in] code Code to assign
 */
    void setCode(size_t dof,
                 int    code) { _code[dof-1] = code; }

/// \brief Replace a node at a given local label
    void Replace(size_t label,
                 Node*  node);

/** \brief Set pointer to neighbor element.
 *  @param [in] el Pointer to element to add as a neigbor element
 *  @remark This function adds the pointer <tt>el</tt> only if this one is not a null pointer
 */
    void Add(Element* el) { if (el) _el[_neig_el++] = el; }

/** \brief Set pointer to neighbor element.
 *  @param [in] el Pointer to element to set as a neighbor element
 *  @param [in] i Local number of neighbor element
 *  @remark This function differs from the Add by the fact that the local
 *  label of neighbor element is given
 */
    void set(Element* el,
             size_t   i) { _el[i-1] = el; }

/// \brief Assign a node given by its pointer as the <tt>i</tt>-th node of side
    void setNode(size_t i,
                 Node*  node);

/// \brief Say that the side is on the boundary
    void setOnBoundary();

//-----------------------------   INSPECTORS  ----------------------------------

/// \brief Return side's shape
    int getShape() const { return _shape; }

/// \brief Return label of side
    size_t getLabel() const { return _label; }
   
/// \brief Return label of side
    size_t n() const { return _label; }

/// \brief Return number of side nodes
    size_t getNbNodes() const { return _nb_nodes; }

/// \brief Return number of side vertices
    size_t getNbVertices() const;

/// \brief Return number of side equations
    size_t getNbEq() const { return _nb_eq; }

/// \brief Return number of DOF
    size_t getNbDOF() const { return _nb_dof; }

/// \brief Return code for a given DOF of node.
/// @param [in] dof Local label of degree of freedom. [Default: <tt>1</tt>]
    int getCode(size_t dof=1) const { return _code[dof-1]; }

/// \brief Return label of <tt>i</tt>-th dof
    size_t getDOF(size_t i) const { return _dof[i-1]; }

/// \brief Return label of first dof of node.
    size_t getFirstDOF() const { return _first_dof; }

/// \brief Return pointer to node of local label <tt>i</tt>.
    Node* getPtrNode(size_t i) const { return _node[i-1]; }

/// \brief Operator ().
/// \details Return pointer to node of local label <tt>i</tt>.
    Node* operator()(size_t i) const { return _node[i-1]; }
   
/// \brief Return global label of node with given local label
    size_t getNodeLabel(size_t i) const { return _node[i-1]->n(); }

/// \brief Return pointer to i-th side neighboring element
///  @param [in] i Local label of neighbor element (must be equal to <tt>1</tt> or <tt>2</tt>).
    Element* getNeighborElement(size_t i) const { return _el[i-1]; }

/** \brief Return pointer to other neighboring element than given one
 *  @param [in] el Pointer to a given neighbor element
 *  @remark If the side is on the boundary this function returns null pointer
 */
    Element* getOtherNeighborElement(Element* el) const;

/** \brief Return normal vector to side
 *  \details The normal vector is oriented from the first neighbor element to
 *  the second one.
 *  \warning The norm of this vector is equal to the measure of the side
 *  (length of the edge in 2-D and area of the face in 3-D), and
 *  To get the unit normal, use rather the member function getUnitNormal.
 */
    Point<real_t> getNormal() const;

/** \brief Return unit normal vector to side
 *  \details The unit normal vector is oriented from the first neighbor
 *  element to the second one.
 *  @remark The norm of this vector is equal to one.
 */
    Point<real_t> getUnitNormal() const;

/** \brief Boundary side or not
 *  \details Returns <tt>1</tt> or <tt>-1</tt> if side is on boundary
 *  Depending on whether the first or the second neighbor element is defined
 *  Returns <tt>0</tt> if side is an inner one
 *  @remark This member function is valid only if member function
 *  \b Mesh::getAllSides() or \b Mesh::getBoundarySides() has been called before.
 */
    int isOnBoundary() const;

/// \brief Say if side has a nonzero code or not
    int isReferenced();

/** \brief Return measure of side.
 *  \details This member function returns length or area of side.
 *  In case of quadrilaterals it returns determinant of Jacobian of
 *  mapping between reference and actual side
 */
    real_t getMeasure() const;

/// \brief Return coordinates of center of side.
    Point<real_t> getCenter() const;

/** \brief Say if a given node belongs to current side
 *  @param [in] nd Pointer to searched node
 *  @return index (local label) of node if found, <tt>0</tt> if not
 */
    size_t Contains(const Node *nd) const;

/// \brief Set side is active (default) or not if argument is <tt>false</tt>
    void setActive(bool opt=true) { _active = opt; }

/// \brief Return <tt>true</tt> or <tt>false</tt> whether side is active or not 
    bool isActive() const { return _active; }

/** \brief Return side level
 *  %Side level increases when side is refined (starting from <tt>0</tt>).
 *  If the level is <tt>0</tt>, then the element has no father
 */
    int getLevel() const { return _level; }

/** \brief Assign side as child of current one and assign current side as father
 *  \details This function is principally used when refining is invoked 
 *  (<i>e.g.</i> for mesh adaption)
 *  @param [in] sd Pointer to side to assign
 */
    void setChild(Side* sd);

/// \brief Return pointer to parent side
/// Return null if no parent
    Side* getParent() const { return _parent; }

/// \brief Return pointer to <tt>i</tt>-th child side
/// Returns null pointer is no childs
    Side* getChild(size_t i) const;

/// \brief Return number of children of side
    size_t getNbChilds() const { return _nb_childs; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    int getGlobalCode() const;
    friend class Mesh;
    friend void Refine(Mesh &in_mesh, Mesh &out_mesh);
    friend ostream& operator<<(ostream&    s,
                               const Mesh& ms);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

private:

    bool           _active;
    vector<int>    _code;
    vector<size_t> _dof;
    mutable int    _on_boundary;
    size_t         _nb_dof, _label, _nb_nodes, _nb_eq, _first_dof, _neig_el;
    size_t         _nb_edges, _nb_childs;
    int            _shape, _level;
    Node           *_node[MAX_NB_SIDE_NODES], *_side_dof[MAX_NB_SIDE_NODES];
    Edge           *_ed[4];
    Element        *_el[2];
    Side           *_child[4], *_parent;
    void shape_index(const string &shape);
};


//-----------------------------------------------------------------------------
// Associated functions
//-----------------------------------------------------------------------------

/** \fn ostream & operator<<(ostream& s, const Side &sd)
 *  \brief Output side data
 *  \ingroup Mesh
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    ostream& operator<<(ostream&    s,
                        const Side& sd);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
