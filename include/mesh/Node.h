/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

                           Definition of for class 'Node'

  ==============================================================================*/

#ifndef __NODE_H
#define __NODE_H

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <iostream>
using std::ostream;
using std::endl;
using std::cerr;
using std::vector;

#include <iomanip>
using std::setw;

#include <string>
using std::string;

#include "OFELI_Config.h"
#include "linear_algebra/Point.h"
#include "io/Fct.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Node.h
 *  \brief Definition file for class Node.
 */

class Element;
class Mesh;

/*! \class Node
 * \ingroup Mesh
 * \brief To describe a node.
 *
 * A node is characterized by its label, its coordinates, its number
 * of degrees of freedom (DOF) and codes that are associated to each DOF.
 *
 * @remark Once the mesh is constructed, information on neighboring elements of node
 * can be retrieved (see appropriate member functions). However, the member function
 * getNodeNeighborElements of Mesh must have been called before. If this is not the
 * case, the program crashes down since no preliminary checking is done for efficiency
 * reasons.
 */

 class Node
{

 public:

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor.
/// \details Initialize data to zero
    Node();

/** \brief Constructor with label and coordinates.
 *  @param [in] label Label of node
 *  @param [in] x %Node coordinates
 */
    Node(size_t               label,
         const Point<real_t>& x);

/// \brief Copy Constructor
    Node(const Node& node);

/// \brief Destructor
    ~Node() { }

//-------------------------------   MODIFIERS  ---------------------------------

/// \brief Define label of node
    void setLabel(size_t label) { _label = label; }

/// \brief Define number of DOF
    void setNbDOF(size_t n);

/// \brief Define First DOF
    void setFirstDOF(size_t n) { _first_dof = n; }

/** \brief Define code for a given DOF of node.
 *  @param [in] dof DOF index
 *  @param [in] code Code to assign to DOF
 */
    void setCode(size_t dof,
                 int    code) { _code[dof-1] = code; }

/// \brief Define codes for all node DOFs.
/// @param [in] code vector instance that contains code for each DOF of current node
    void setCode(const vector<int>& code)
    {
       for (size_t i=0; i<_nb_dof; i++)
          _code[i] = code[i];
    }

/// \brief Define codes for all node DOFs.
/// @param [in] code C-array that contains code for each DOF of current node
    void setCode(int* code)
    {
       for (size_t i=0; i<_nb_dof; i++)
          _code[i] = code[i];
    }

/** \brief Define code by a boolean algebraic expression invoking node coordinates
 *  @param [in] exp Boolean algebraic expression as required by <tt>fparser</tt>
 *  @param [in] code Code to assign to node if the algebraic expression is true
 *  @param [in] dof Degree of Freedom for which code is assigned [Default: 1]
 */
    void setCode(const string& exp,
                 int           code,
                 size_t        dof=1);

/** \brief Set i-th coordinate.
 *  @param [in] i Coordinate index (1..3)
 *  @param [in] x Coordinate value
 */
    void setCoord(size_t i,
                  real_t x) { _x[i-1] = x; }

/** \brief Define label of DOF.
 *  @param [in] i DOF index
 *  @param [in] dof Label of DOF
 */
    void DOF(size_t i,
             size_t dof) { _dof[i-1] = dof; }

/** \brief Define number of DOF.
 *  @param [in,out] first_dof Label of the first DOF in input that is actualized
 *  @param [in] nb_dof Number of DOF
 */
    void setDOF(size_t& first_dof,
                size_t  nb_dof);

/// \brief Set node as boundary node.
/// \details This function is mostly internally used (Especially in class Mesh)
    void setOnBoundary() { _on_boundary = true; }

//-----------------------------   INSPECTORS  ----------------------------------

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Return label of node
    size_t getLabel() const { return _label; }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return label of node
    size_t n() const { return _label; }

/// \brief Return number of degrees of freedom (DOF)
    size_t getNbDOF() const { return _nb_dof; }

/// \brief Return code for a given DOF of node.
/// @param [in] dof label of degree of freedom for which code is to be returned. Default value is 1.
    int getCode(size_t dof) const { return _code[dof-1]; }

/// \brief Return global code of node.
    int getCode() const;

/// \brief Return i-th coordinate of node.
/// i = 1..3
    real_t getCoord(size_t i) const { return _x[i-1]; }

/// \brief Return coordinates of node.
/// \details Return value is an instance of class Point
    Point<real_t> getCoord() const { return Point<real_t>(_x[0],_x[1],_x[2]); }

/// \brief Return x-coordinate of node.
    real_t getX() const { return _x[0]; }

/// \brief Return y-coordinate of node.
    real_t getY() const { return _x[1]; }

/// \brief Return z-coordinate of node.
    real_t getZ() const { return _x[2]; }

/// \brief Return coordinates of node.
/// \details Return value is an instance of class Point
    Point<real_t> getXYZ() const { return Point<real_t>(_x[0],_x[1],_x[2]); }

/// \brief Return label of i-th dof.
    size_t getDOF(size_t i) const { return _dof[i-1]; }

/** \brief Return number of neighbor elements
 *  \details Neighbor elements are those that share node.
 *  Note that the returned information is valid only if the Mesh member function
 *  \b getNodeNeighborElements() has been invoked before
 */
    size_t getNbNeigEl() const { return _nb_neig_el; }

/// \brief Return i-th neighbor element
/// \details Note that the returned information is valid only if the Mesh member function
/// \b getNodeNeighborElements() has been invoked before
    Element* getNeigEl(size_t i) const { return _el[i-1]; }

/// \brief Return label of first DOF of node.
    size_t getFirstDOF() const { return _first_dof; }

/** \brief Say if node is a boundary node.
 *  \details Note this information is available only if boundary sides (and nodes)
 *  were determined (See class Mesh).
 */
    bool isOnBoundary() const { return _on_boundary; }

/// \brief Add element pointed by <tt>el</tt> as neighbor element to node
    void Add(Element* el);

/** \brief Assign a level to current node
 *  \details This member function is useful for mesh adaption.\n
 *  Default node's level is zero
 */
    void setLevel(int level) { _level = level; }

/** \brief Return node level
 *  \details Node level decreases when element is refined (starting from 0).
 *  If the level is 0, then the element has no parents
 */
    int getLevel() const { return _level; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    friend class Mesh;
    friend void Refine(Mesh &in_mesh, Mesh &out_mesh);
    friend ostream& operator<<(ostream&    s,
                               const Mesh& ms);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 private:

   size_t            _nb_dof, _first_dof, _label, _nb_neig_el, _neig_i, _dof[MAX_NBDOF_NODE];
   int               _level, _code[MAX_NBDOF_NODE];
   Point<real_t>     _x;
   vector<Element *> _el;
   Fct               _theFct;
   bool              _on_boundary;
   void Neig() { _nb_neig_el++; }
   void Add() { }
      //   void Add() { _el = new Element * [_nb_neig_el]; }
};


//-----------------------------------------------------------------------------
// Associated functions
//-----------------------------------------------------------------------------

/** \fn ostream & operator<<(ostream &s, const Node &nd)
 *  \brief Output node data.
 *  \ingroup Mesh
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    ostream & operator<<(ostream&    s,
                         const Node& nd);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
