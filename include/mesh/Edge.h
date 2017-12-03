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

                            Definition of class 'Edge'

  ==============================================================================*/


#ifndef __EDGE_H
#define __EDGE_H

#include <stdlib.h>
#include <string>

#include "OFELI_Config.h"
#include "mesh/Node.h"
#include "io/output.h"

using std::string;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Edge.h
 *  \brief Definition file for class Edge.
 */

class Side;

/*! \class Edge
 *  \ingroup Mesh
 * \brief To describe an edge.
 *
 * \details Defines an edge of a 3-D finite element mesh. The edges are given
 * in particular by a list of nodes. Each node can be accessed by the
 * member function getPtrNode.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class Edge
{

 public:

//----------------------------   BASIC OPERATIONS   -----------------------------

/// \brief Default Constructor.
/// \details Initializes data to zero
    Edge();

/// \brief Constructor with label.
/// \details Define an edge by giving its <tt>label</tt>
    Edge(size_t label);

/// \brief Copy constructor.
    Edge(const Edge& ed);

/// \brief Destructor.
    ~Edge();

//-------------------------------   MODIFIERS  ---------------------------------

/// \brief Insert a node at end of list of nodes of edge.
    void Add(Node* node);

/// \brief Assign label of edge.
    void setLabel(size_t i) { _label = i; }

/// \brief Define First DOF
    void setFirstDOF(size_t n) { _first_dof = n; }

/// \brief Define number of DOF of edge
    void setNbDOF(size_t nb_dof) { _nb_dof = nb_dof; }
   
/** \brief Define label of DOF.
 *  @param [in] i DOF index
 *  @param [in] dof Its label
 */
    void DOF(size_t i,
             size_t dof)
    { _dof[i-1] = dof; }
   
/** \brief Define number of DOF.
 *  @param [in,out] first_dof Label of the first DOF in input that is actualized
 *  @param [in] nb_dof Number of DOF
 */
    void setDOF(size_t& first_dof,
                size_t  nb_dof);
   
/** \brief Assign code <tt>code</tt> to DOF number <tt>dof</tt>.
    @param [in] dof index of dof for assignment.
    @param [in] code Value of code to assign.
*/
    void setCode(size_t dof,
                 int    code)
    { _code[dof-1] = code; }

/// \brief Add side pointed by <tt>sd</tt> to list of edge sides.
    void AddNeighbor(Side* sd)
    { _sd[_neig_sd++] = sd; }

//-----------------------------   INSPECTORS  ----------------------------------

/// \brief Return label of edge
    size_t getLabel() const { return _label; }
   
/// \brief Return label of edge
    size_t n() const { return _label; }
   
/// \brief Return number of edge equations
    size_t getNbEq() const { return _nb_eq; }

/// \brief Return number of DOF.
    size_t getNbDOF() const { return _nb_dof; }

/// \brief Return code for a given DOF of node.
/// \details Default value is <tt>1</tt>
    int getCode(size_t dof=1) const { return _code[dof-1]; }

/// \brief Return label of <tt>i</tt>-th DOF.
    size_t getDOF(size_t i) const { return _dof[i-1]; }

/// \brief Return number of first dof of node
    size_t getFirstDOF() const { return _first_dof; }

/// \brief List of element nodes
    Node* getPtrNode(size_t i) const { return _node[i-1]; }
   
/// \brief Operator ().
/// \details Return pointer to node of local label <tt>i</tt>.
    Node* operator()(size_t i) const { return _node[i-1]; }
   
/// \brief Return node label.
/// @param [in] i Local label of node for which global label is returned
    size_t getNodeLabel(size_t i) const { return _node[i-1]->n(); }

/// \brief Return pointer to neighbor i-th side.
    Side* getNeighborSide(size_t i) const { return _sd[i-1]; }

/** \brief Say if current edge is a boundary edge or not.
 *  \details Note this information is available only if boundary edges were determined. 
 *  See class Mesh
 */
    int isOnBoundary() const;

/// \brief Say that the edge is on the boundary
    void setOnBoundary();

/// \brief Operator ().
/// \details Returns pointer to node of local label <tt>i</tt>
    Node *operator()(size_t i) { return getPtrNode(i); }

 private:

   std::vector<int>    _code;
   std::vector<size_t> _dof;
   mutable int         _on_boundary;
   size_t              _nb_nodes, _nb_dof, _label, _nb_eq, _first_dof, _neig_sd;
   Node                *_node[2];
   Side                *_sd[2];
};


/// \fn ostream & operator<<(ostream& s, const Edge &ed)
/// \brief Output edge data.
/// \ingroup Mesh
    ostream& operator<<(ostream&    s,
                        const Edge& ed);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */


#endif
