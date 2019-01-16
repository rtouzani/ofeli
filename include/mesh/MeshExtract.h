/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

               Definition of classes for extracting submeshes

  ==============================================================================*/

#ifndef __MESH_EXTRACT_H
#define __MESH_EXTRACT_H

#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file MeshExtract.h
 *  \brief Definition file for classes for extracting submeshes.
 */


/*! \class NodeList
 *  \ingroup Mesh
 *  \brief Class to construct a list of nodes having some common properties.
 *
 * This class enables choosing multiple selection criteria by using
 * function <tt>select...</tt> However, the intersection of these properties 
 * must be empty.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class NodeList {

 public:

/// \brief Constructor using a Mesh instance
    NodeList(Mesh& ms);

/// \brief Destructor
    ~NodeList() { }

/** \brief Select nodes having a given code for a given degree of freedom
 *  @param [in] code Code that nodes share
 *  @param [in] dof Degree of Freedom label [Default: <tt>1</tt>]
 */
    void selectCode(int code,
                    int dof=1);

/** \brief Unselect nodes having a given code for a given degree of freedom
 *  @param [in] code Code of nodes to exclude
 *  @param [in] dof Degree of Freedom label [Default: <tt>1</tt>]
 */
    void unselectCode(int code,
                      int dof=1);

/** \brief Select nodes having given coordinates
 *  @param [in] x <tt>x</tt>-coordinate that share the selected nodes
 *  @param [in] y <tt>y</tt>-coordinate that share the selected nodes [Default: <tt>ANY</tt>]
 *  @param [in] z <tt>z</tt>-coordinate that share the selected nodes [Default: <tt>ANY</tt>]
 *
 *  Coordinates can be assigned the value <tt>ANY</tt>. This means that any coordinate
 *  value is accepted. For instance, to select all nodes with <tt>x=0</tt>, use
 *  \b selectCoordinate(0.,ANY,ANY);
 */
    void selectCoordinate(real_t x,
                          real_t y=ANY,
                          real_t z=ANY);

/// \brief Return number of selected nodes
    size_t getNbNodes() const { return _nb_nodes; }

/// \brief Reset list of nodes at its top position (Non constant version)
    void top() { _it = 0; }

/// \brief Reset list of nodes at its top position (Constant version)
    void top() const { _const_it = 0; }

/// \brief Return pointer to current node and move to next one (Non constant version)
    Node *get()
    {
       if (_it==_nb_nodes)
          return NULL;
       else
          return _theList[_it++];
    }

/// \brief Return pointer to current node and move to next one (Constant version)
    Node *get() const
    {
       if (_const_it==_nb_nodes)
          return NULL;
       else
          return _theList[_const_it++];
    }

 private:

   Mesh *_ms;
   NodeSet _theList;
   size_t _it, _nb_nodes, _nbn;
   mutable size_t _const_it;
};


/*! \class ElementList
 *  \ingroup Mesh
 * \brief Class to construct a list of elements having some common properties.
 *
 * \details This class enables choosing multiple selection criteria by using
 * function <tt>select...</tt> However, the intersection of these properties 
 * must be empty.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class ElementList {

 public:

/// \brief Constructor using a Mesh instance
    ElementList(Mesh& ms);

/// \brief Destructor
    ~ElementList() { }

/// \brief Select elements having a given code
    void selectCode(int code);

/// \brief Unselect elements having a given code
/// @param [in] code Code of elements to exclude
    void unselectCode(int code);

/** \brief Select elements having a given level
 *  @param [in] level Level of elements to select
 *
 *  Elements having a given level (for mesh adaption) are selected in a list
 */
    void selectLevel(int level);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Select active elements
/// \details Useful for mesh adaption: Only active elements are selected
    void selectActive();

/** \brief Select active elements
 *  @param [in] e Vector that contains for each element the estimator value
 *  @param [in] threshold If the estimator in an element is above this value
 *  then this element is selected
 */
    void selectMarked(const Vect<real_t>& e,
                            real_t        threshold);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Return number of selected elements
    size_t getNbElements() const { return _nb_elements; }

/// \brief Reset list of elements at its top position (Non constant version)
    void top() { _it = 0; }

/// \brief Reset list of elements at its top position (Constant version)
    void top() const { _const_it = 0; }

/// \brief Return pointer to current element and move to next one (Non constant version)
    Element *get()
    {
       if (_it==_nb_elements)
          return NULL;
       else
          return _theList[_it++];
    }

/// \brief Return pointer to current element and move to next one (Constant version)
    Element *get() const
    {
       if (_const_it==_nb_elements)
          return NULL;
       else
          return _theList[_const_it++];
    }

 private:

   Mesh *_ms;
   ElementSet _theList;
   size_t _it, _nb_elements, _nbn;
   mutable size_t _const_it;
};


/*! \class SideList
 *  \ingroup Mesh
 * \brief Class to construct a list of sides having some common properties.
 *
 * \details This class enables choosing multiple selection criteria by using
 * function <tt>select...</tt> However, the intersection of these properties 
 * must be empty.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class SideList {

 public:

/// \brief Constructor using a Mesh instance
    SideList(Mesh& ms);

/// \brief Destructor
    ~SideList() { }

/** \brief Select sides having a given code for a given degree of freedom
 *  @param [in] code Code that sides share
 *  @param [in] dof Degree of Freedom label [Default: <tt>1</tt>]
 */
    void selectCode(int code,
                    int dof=1);

/** \brief Unselect sides having a given code for a given degree of freedom
 *  @param [in] code Code of sides to exclude
 *  @param [in] dof Degree of Freedom label [Default: <tt>1</tt>]
 */
    void unselectCode(int code,
                      int dof=1);

/// \brief Return number of selected sides
    size_t getNbSides() const { return _nb_sides; }

/// \brief Reset list of sides at its top position (Non constant version)
    void top() { _it = 0; }

/// \brief Reset list of sides at its top position (Constant version)
    void top() const { _const_it = 0; }

/// \brief Return pointer to current side and move to next one (Non constant version)
    Side *get()
    {
       if (_it==_nb_sides)
          return NULL;
       else
          return _theList[_it++];
    }

/// \brief Return pointer to current side and move to next one (Constant version)
    Side *get() const
    {
       if (_const_it==_nb_sides)
          return NULL;
       else
          return _theList[_const_it++];
    }

 private:

   Mesh *_ms;
   SideSet _theList;
   size_t _it, _nb_sides, _nbn;
   mutable size_t _const_it;
};


/*! \class EdgeList
 *  \ingroup Mesh
 * \brief Class to construct a list of edges having some common properties.
 *
 * \details This class enables choosing multiple selection criteria by using
 * function <tt>select...</tt> However, the intersection of these properties 
 * must be empty.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class EdgeList {

 public:

/// \brief Constructor using a Mesh instance
    EdgeList(Mesh& ms);

/// \brief Destructor
    ~EdgeList() { }

/** \brief Select edges having a given code for a given degree of freedom
 *  @param [in] code Code that edges share
 *  @param [in] dof Degree of Freedom label [Default: <tt>1</tt>]
 */
    void selectCode(int code,
                    int dof=1);

/** \brief Unselect edges having a given code for a given degree of freedom
 *  @param [in] code Code of edges to exclude
 *  @param [in] dof Degree of Freedom label [Default: <tt>1</tt>]
 */
    void unselectCode(int code,
                      int dof=1);

/// \brief Return number of selected edges
    size_t getNbEdges() const { return _nb_edges; }

/// \brief Reset list of edges at its top position (Non constant version)
    void top() { _it = 0; }

/// \brief Reset list of edges at its top position (Constant version)
    void top() const { _const_it = 0; }

/// \brief Return pointer to current edge and move to next one (Non constant version)
    Edge *get()
    {
       if (_it==_nb_edges)
          return NULL;
       else
          return _theList[_it++];
    }

/// \brief Return pointer to current edge and move to next one (Constant version)
    Edge *get() const
    {
       if (_const_it==_nb_edges)
          return NULL;
       else
          return _theList[_const_it++];
    }

 private:

    Mesh *_ms;
    EdgeSet _theList;
    size_t _it, _nb_edges, _nbn;
    mutable size_t _const_it;
};


/** \fn ostream & operator<<(ostream& s, const NodeList& nl)
 *  \brief Output NodeList instance
 *  \ingroup Mesh
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    ostream& operator<<(ostream&        s,
                        const NodeList& nl);

/** \fn ostream & operator<<(ostream& s, const ElementList &el)
 *  \brief Output ElementList instance
 *  \ingroup Mesh
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    ostream& operator<<(ostream&           s,
                        const ElementList& el);

/** \fn ostream & operator<<(ostream& s, const SideList& sl)
 *  \brief Output SideList instance
 *  \ingroup Mesh
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    ostream& operator<<(ostream&        s,
                        const SideList& sl);
   
/** \fn ostream & operator<<(ostream& s, const EdgeList& el)
 * \brief Output EdgeList instance
 * \ingroup Mesh
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
    ostream& operator<<(ostream&        s,
                        const EdgeList& el);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
