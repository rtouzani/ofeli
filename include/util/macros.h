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

                                    OFELI's macros

  ==============================================================================*/

#ifndef __MACROS_H
#define __MACROS_H

#include "OFELI_Config.h"
using std::exception;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define MESH_NODES          for (topNode(); (theNode=getNode());)
#define MESH_ELEMENTS       for (topElement(); (theElement=getElement());)
#define MESH_SIDES          for (topSide(); (theSide=getSide());)
#define MESH_BOUNDARY_SIDES for (topBoundarySide(); (theSide=getBoundarySide());)
#define MESH_EDGES          for (topEdge(); (theEdge=getEdge());)
#define MAT(M,A)            (M &)(*A)
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \def CATCH_EXCEPTION
 *  \ingroup Util
 *  This macro can be inserted after a <tt>try</tt> loop to catch a thrown 
 *  exception.
 */
#define CATCH_EXCEPTION                                         \
   catch(OFELIException &e) {                                   \
      std::cout << "OFELI error: " << e.what() << endl;         \
      return 1;                                                 \
   }                                                            \
   catch(runtime_error &e) {                                    \
      std::cout << "OFELI Runtime error: " << e.what() << endl; \
      return 1;                                                 \
   }                                                            \
   catch( ... ) {                                               \
      std::cout << "OFELI Unexpected error: " << endl;          \
      return 1;                                                 \
   }

/*! \def TheNode
 *  \ingroup Mesh
 *  A macro that gives the instance pointed by \a theNode
 */
#define TheNode (*theNode)

/*! \def TheElement
 *  \ingroup Mesh
 * A macro that gives the instance pointed by \a theElement
 */
#define TheElement (*theElement)

/*! \def TheSide
 *  \ingroup Mesh
 * A macro that gives the instance pointed by \a theSide
 */
#define TheSide (*theSide)

/*! \def TheEdge
 *  \ingroup Mesh
 * A macro that gives the instance pointed by \a theEdge
 */
#define TheEdge (*theEdge)

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#define node_loop(m) for (auto const& the_node: (m)->theNodes)
#define boundary_node_loop(m) for (auto const& the_node: (m)->theBoundaryNodes)
#define element_loop(m) for (auto const& the_element: (m)->theElements)
#define side_loop(m) for (auto const& the_side: (m)->theSides)
#define boundary_side_loop(m) for (auto const& the_side: (m)->theBoundarySides)
#define edge_loop(m) for (auto const& the_edge: (m)->theEdges)

#define MESH_EL element_loop(_theMesh)
#define MESH_ND node_loop(_theMesh)
#define MESH_SD side_loop(_theMesh)
#define MESH_BD_SD boundary_side_loop(_theMesh)
#define MESH_BD_ND boundary_node_loop(_theMesh)
#define MESH_ED edge_loop(_theMesh)

#define The_element (*the_element)
#define The_node (*the_node)
#define The_side (*the_side)
#define The_edge (*the_edge)

#define element_label (the_element->n())
#define node_label (the_node->n())
#define side_label (the_side->n())
#define edge_label (the_edge->n())

#define SQR(x) ((x)*(x))

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \def ElementLoop(m)
 *  \ingroup Mesh
 * A macro to loop on mesh elements
 * \a m : Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theElement</tt> to current Element
 */
#define ElementLoop(m) for (auto const& theElement: (m).theElements)

/*! \def SideLoop(m)
 *  \ingroup Mesh
 * A macro to loop on mesh sides
 * \a m : Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theSide</tt> to current Element
 */
#define SideLoop(m) for (auto const& theSide: (m).theSides)

/*! \def EdgeLoop(m)
 *  \ingroup Mesh
 * A macro to loop on mesh edges
 * \a m : Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theEdge</tt> to current Edge
 */
#define EdgeLoop(m) for (auto const& theEdge: (m).theEdges)

/*! \def NodeLoop(m)
 *  \ingroup Mesh
 * A macro to loop on mesh nodes
 * \a m : Instance of Mesh
 * @note: Each iteration updates the pointer \a theNode to current Node
 */
#define NodeLoop(m) for (auto const& theNode: (m).theNodes)

/*! \def BoundaryNodeLoop(m)
 *  \ingroup Mesh
 * A macro to loop on mesh nodes
 * <tt>m</tt>: Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theNode</tt> to current Node
 */
#define BoundaryNodeLoop(m) for (auto const& theNode: (m).theBoundaryNodes)

/*! \def BoundarySideLoop(m)
 *  \ingroup Mesh
 * A macro to loop on mesh boundary sides
 * <tt>m</tt>: Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theSide</tt> to current Node
 */
#define BoundarySideLoop(m) for (auto const& theSide: (m).theBoundarySides)

/*! \def ElementNodeLoop
 *  \ingroup Mesh
 * A macro to loop on element nodes
 * @param el Instance of Element
 * @param nd Pointer to pointed node
 */
#define ElementNodeLoop(el,nd) for (auto const& nd: (el).theNodes)

/*! \def ElementSideLoop
 *  \ingroup Mesh
 * A macro to loop on element sides
 * @param el Instance of Element
 * @param sd Pointer to pointed side
 */
#define ElementSideLoop(el,sd) for (auto const& sd: (el).theSides)

/*! \def theNodeLabel
 *  \ingroup Mesh
 * A macro that returns node label in a loop using macro \a MeshNodes
 */
#define theNodeLabel theNode->n()

/*! \def theSideLabel
 *  \ingroup Mesh
 * \brief A macro that returns side label in a loop using macro <tt>MeshSides</tt>
 */
#define theSideLabel theSide->n()

/*! \def theSideNodeLabel(i)
 *  \ingroup Mesh
 * \brief A macro that returns label of <tt>i</tt>-th node of side using macro <tt>MeshSides</tt>
 */
#define theSideNodeLabel(i) theSide->getNodeLabel(i)

/*! \def theElementLabel
 *  \ingroup Mesh
 * \brief A macro that returns element label in a loop using macro <tt>MeshElements</tt>
 */
#define theElementLabel theElement->n()

/*! \def theElementNodeLabel(i)
 *  \ingroup Mesh
 * \brief A macro that returns label of <tt>i</tt>-th node of element using macro <tt>MeshElements</tt>
 */
#define theElementNodeLabel(i) theElement->getNodeLabel(i)


/*! \def TIME_LOOP
 *  \ingroup Solver
 * \brief A macro to loop on time steps to integrate on time
 * \a ts : Time step
 * \a t  : Initial time value updated at each time step
 * \a ft : Final time value
 * \a n  : Time step index
 */
#define TIME_LOOP(ts,t,ft,n)   \
          n = 1;               \
          for (real_t t=0; t<ft+0.01*ts; t+=ts, ++n)


/*! \def TimeLoop
 *  \ingroup Solver
 * \brief A macro to loop on time steps to integrate on time.
 * \details It uses the following global variables defined in \b OFELI:
 * <tt>theStep, theTime, theTimeStep, theFinalTime</tt>
 */
#define TimeLoop \
          OFELI::NbTimeSteps = int(OFELI::theFinalTime/OFELI::theTimeStep);   \
          for (OFELI::theTime=OFELI::theTimeStep, theStep=1;                  \
               theTime<OFELI::theFinalTime+0.001*OFELI::theTimeStep;          \
               OFELI::theTime+=OFELI::theTimeStep, ++OFELI::theStep)

/*! \def IterationLoop
 *  \ingroup Solver
 * \brief A macro to loop on iterations for an iterative procedure.
 * \details It uses the following global variables defined in \b OFELI:
 * <tt>theIteration, MaxNbIterations, Converged</tt>
 * @warning The variable <tt>theIteration</tt> must be zeroed before using this macro 
 */
#define IterationLoop \
          while (++theIteration<MaxNbIterations && Converged==false)

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
