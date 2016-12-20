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

                           Definition of class 'Element'

  ==============================================================================*/

#ifndef __ELEMENT_H
#define __ELEMENT_H

#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
using std::string;

#include <map>
using std::map;

#include "OFELI_Config.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

  
/*! \file Element.h
 *  \brief Definition file for class Element.
 */

class Node;
class Side;
class Mesh;

/** \class Element
 *  \ingroup Mesh
 *  \brief To store and treat finite element geometric information
 *
 *  Class Element enables defining an element of a finite element mesh. The element is given in
 *  particular by its shape and a list of nodes. Each node can be accessed by the
 *  member function getPtrNode. Moreover, class Mesh can generate for each element its list of sides.
 *  The string that defines the element shape must be chosen according to the following list :
 *
 *  \htmlinclude "ElementDescription.html"
 *
 *  @remark Once a Mesh instance is constructed, one has access for each Element of the mesh
 *  to pointers to element sides provided the member function getAllSides of Mesh has been
 *  invoked. With this, an element can be tested to see if it is on the boundary, i.e. if it
 *  has at least one side on the boundary
 */

class Element
{

 public:

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor
    Element();

/** \brief Constructor initializing label, shape of element.
 *  @param [in] label Label to assign to element.
 *  @param [in] shape Shape of element (See class description).
 */
    Element(size_t        label,
            const string& shape);

/** \brief Constructor initializing label, shape of element.
 *  @param [in] label Label to assign to element.
 *  @param [in] shape Shape of element (See enum ElementShape in Mesh)
 */
    Element(size_t label,
            int    shape);

/** \brief Constructor initializing label, shape and code of element.
 *  @param [in] label Label to assign to element.
 *  @param [in] shape Shape of element (See class description).
 *  @param [in] c Code to assign to element (useful for media properties).
 */
    Element(size_t        label,
            const string& shape,
            int           c);

/** \brief Constructor initializing label, shape and code of element.
 *  @param [in] label Label to assign to element.
 *  @param [in] shape Shape of element (See enum ElementShape in Mesh).
 *  @param [in] c Code to assign to element (useful for media properties).
 */
    Element(size_t label,
            int    shape,
            int    c);

/// \brief Copy constructor
    Element(const Element& el);

/// \brief Destructor
    ~Element() { ; }

//-------------------------------   MODIFIERS  ---------------------------------

/// \brief Define label of element.
/// @param [in] i Label to assign to element
    void setLabel(size_t i) { _label = i; }

/// \brief Define code of element.
/// @param [in] c Code to assign to element.
    void setCode(int c) { _code = c; }

/** \brief Define code by a boolean algebraic expression invoking coordinates of element nodes
 *  @param [in] exp Boolean algebraic expression as required by <tt>fparser</tt>
 *  @param [in] code Code to assign to node if the algebraic expression is true
 */
    void setCode(const string& exp,
                 int           code);

/// \brief Insert a node at end of list of nodes of element.
/// @param [in] node Pointer to Node instance.
    void Add(Node* node);

/// \brief Insert a node and set its local node number.
/// @param node [in] Pointer to Node instance
/// @param [in] n Element node number to assign
    void Add(Node* node,
             int   n);

/// \brief Replace a node at a given local label.
/// @param [in] label Node to replace.
/// @param [in] node Pointer to Node instance to copy to current instance.
    void Replace(size_t label,
                 Node*  node);

/// \brief Replace a side at a given local label.
/// @param [in] label Side to replace.
/// @param [in] side Pointer to Side instance to copy to current instance.
    void Replace(size_t label,
                 Side*  side);

/// \brief Assign Side to Element.
/// @param [in] sd Pointer to Side instance.
    void Add(Side* sd);

/// \brief Assign Side to Element with assigned local label.
/// @param [in] sd Pointer to Side instance.
/// @param [in] k Local label.
    void Add(Side* sd,
             int   k);

/// \brief Add a neighbor element.
/// @param [in] el Pointer to Element instance
    void Add(Element* el);

/// \brief Add a neighbor element and set its label.
/// @param [in] el Pointer to Element instance
/// @param [in] n Neighbor element number to assign
    void set(Element* el,
             int      n);

/// \brief Define label of DOF.
/// @param [in] i Index of DOF.
/// @param [in] dof Label of DOF to assign.
    void setDOF(size_t i,
                size_t dof)
    { _dof[i-1] = dof; }

/// \brief Assign code to a DOF.
/// @param [in] dof Index of dof for assignment.
/// @param [in] code Code to assign.
    void setCode(size_t dof,
                 int    code)
    { _dof_code[dof-1] = code; }

/// \brief Assign a node given by its pointer as the i-th node of element
    void setNode(size_t i,
                 Node*  node);

/// \brief Set number of degrees of freedom of element
    void setNbDOF(size_t i) { _nb_dof = i; }

/// \brief Set label of first DOF in element
    void setFirstDOF(size_t i) { _first_dof = i; }

//-----------------------------   INSPECTORS  ----------------------------------

/// \brief Return element shape
    int getShape() const { return _shape; }

/// \brief Return label of element
    size_t getLabel() const { return _label; }

/// \brief Return label of element
    size_t n() const { return _label; }

/// \brief Return code of element
    int getCode() const { return _code; }

/// \brief Return number of element nodes
    size_t getNbNodes() const { return _nb_nodes; }

/// \brief Return number of element vertices
    size_t getNbVertices() const;

/// \brief Return number of element sides (Constant version)
    size_t getNbSides() const { return _nbs; }

/// \brief Return number of element equations
    size_t getNbEq() const { return _nb_eq; }

/// \brief return element nb of DOF
    size_t getNbDOF() const;

/// \brief Return element DOF label
    size_t getDOF(size_t i=1) const { return _first_dof+i-1; }

/// \brief Return element first DOF label
    size_t getFirstDOF() const { return _first_dof; }

/// \brief Return global label of node of local label <tt>i</tt>.
    size_t getNodeLabel(size_t n) const;

/// \brief Return global label of side of local label <tt>i</tt>.
    size_t getSideLabel(size_t n) const;

/// \brief Return pointer to node of label <tt>i</tt> (Local labelling).
    Node *getPtrNode(size_t i) const;
   
/// \brief Operator ().
/// \details Return pointer to node of local label <tt>i</tt>.
    Node *operator()(size_t i) const { return getPtrNode(i); }
   
/// \brief Return pointer to side of label <tt>i</tt> (Local labelling).
    Side *getPtrSide(size_t i) const;

/** \brief Say if element contains given node
 *  \details This function tests if the element contains a node with the same pointer at the sought one
 *  @param [in] nd Pointer to Node instance
 *  @return Local node label in element. If <tt>0</tt>, the element does not contain this node
 */
    int Contains(const Node* nd) const;
   
/** \brief Say if element contains given node
 *  \details This function tests if the element contains a node with the same label at the sought one
 *  @param [in] nd Reference to Node instance
 *  @return Local node label in element. If <tt>0</tt>, the element does not contain this node
 */
    int Contains(const Node& nd) const;

/** \brief Say if element contains given side
 *  \details This function tests if the element contains a side with the same pointer at the sought one
 *  @param [in] sd Pointer to Side instance
 *  @return Local side label in element. If <tt>0</tt>, the element does not contain this side
 */
    int Contains(const Side* sd) const;
   
/** \brief Say if element contains given side
 *  \details This function tests if the element contains a side with the same label at the sought one
 *  @param [in] sd Reference to Side instance
 *  @return Local side label in element. If <tt>0</tt>, the element does not contain this side
 */
    int Contains(const Side& sd) const;

/** \brief Return pointer to element Neighboring element.
 *  @param [in] i Index of element to look for.
 *  @note This method returns valid information only if the Mesh member function
 *  Mesh::getElementNeighborElements() has been called before.
 */
    Element *getNeighborElement(size_t i) const { return _neig_el[i-1]; }

/** \brief Return number of neigboring elements.
 *  \details @note This method returns valid information only if the Mesh member function
 *  Mesh::getElementNeighborElements() has been called before.
 */
    size_t getNbNeigElements() const { return _nb_neig_el; }

/** \brief Return measure of element.
 *  \details This member function returns length, area or volume of element.
 *  In case of quadrilaterals and hexahedrals it returns determinant of Jacobian of
 *  mapping between reference and actual element
 */
    real_t getMeasure() const;

/// \brief Return outward unit normal to i-th side of element.
/// \details Sides are ordered [node_1,node_2], [node_2,node_3], ...
    Point<real_t> getUnitNormal(size_t i) const;

/** \brief Say if current element is a boundary element or not.
 *  @note this information is available only if boundary elements were determined
 *  i.e. if member function Mesh::getBoundarySides or Mesh::getAllSides has
 *  been invoked before.
 */
    bool isOnBoundary() const;

/// \brief Operator ().
/// \details Return pointer to node of local label <tt>i</tt>.
    Node *operator()(size_t i) { return getPtrNode(i); }

/** \brief Initialize information on element sides.
 *  \details This function is to be used to initialize loops over sides.
 *  @param [in] n Label of side.
 *  @param [in] nd Array of pointers to nodes of the side (<tt>nd[0], nd[1], ...</tt> point
 *  to first, second nodes, ...
 */
    int setSide(size_t  n,
                size_t* nd);

/// \brief Return <tt>true</tt> or <tt>false</tt> whether element is active or not 
    bool isActive() const;

/** \brief Return element level
 *  Element level decreases when element is refined (starting from 0).
 *  If the level is 0, then the element has no father
 */
    int getLevel() const { return _level; }

/** \brief Assign element as child of current one and assign current element as father
 *  This function is principally used when refining is invoked (e.g. for mesh adaption)
 *  @param [in] el Pointer to element to assign
 */
    void setChild(Element* el);

/// \brief Return pointer to i-th child element
/// Return null pointer is no childs
    Element *getChild(size_t i) const;

/// \brief Return number of children of element
    size_t getNbChilds() const { return _nb_childs; }

/// \brief Return pointer to parent element
/// Return null if no parent
    Element *getParent() const { return _parent; }

/** \brief Check if a given node belongs to current element
 *  @param [in] nd Pointer to node to locate
 *  @return local label of node if this one is found, 0 otherwise
 */
    size_t IsIn(const Node* nd);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    void setGlobalToLocal();
    size_t getLocal(size_t i) { return _g2l[i]; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

    friend class Mesh;
    friend void Refine(Mesh &in_mesh, Mesh &out_mesh);
    friend ostream& operator<<(ostream& s, const Mesh& ms);

 private :

   bool _active;
   size_t _nb_nodes, _nb_eq, _nb_sides, _label, _nbs, _nb_neig_el, _nb_childs;
   size_t _nb_dof, _first_dof;
   int _level, _code, _shape, _dof_code[MAX_NBDOF_SIDE];
   size_t _dof[MAX_NBDOF_SIDE];
   Node *_node[MAX_NB_ELEMENT_NODES];
   Side *_side[MAX_NB_ELEMENT_SIDES];
   Element *_neig_el[MAX_NB_ELEMENT_SIDES];
   Element *_child[4], *_parent;
   map<size_t,size_t> _g2l;

   void calculate_nb_sides();
   void shape_index(const string &shape);
};

///////////////////////////////////////////////////////////////////////////////
//                             Associated Functions                          //
///////////////////////////////////////////////////////////////////////////////

/** \fn ostream & operator<<(ostream& s, const Element &el)
 *  \ingroup Mesh
 *  \brief Output element data
 */
    ostream& operator<<(ostream&       s,
                        const Element& el);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
