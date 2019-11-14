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


/*! \def PARSE(exp,var)
 *  \ingroup Util
 * A macro that parses a regular expression \a exp using
 * the variables in the string \a var.
 * For instance, to parse the function \a sin(x+y) one must
 * declare \a PARSE("sin(x+y)","x,y")
 */
#define PARSE(exp,var) theParser.Parse(exp,var)

/*! \def EVAL(d)
 *  \ingroup Util
 * A macro that evaluates a parsed regular expression
 * For instance, with a declaration \a PARSE("sin(x+y)","x,y")
 * the data \a x=1 and \a y=2 using this function must
 * be evaluated as follows: \a EVAL(d) with \a d[0]=1,
 * \a d[1]=2
 */
#define EVAL(d)    theParser.Eval(d)


#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define EVAL_ERR theParser.EvalError()
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \def CATCH_EXCEPTION
 *  \ingroup Util
 *  This macro can be inserted after a <tt>try</tt> loop to catch a thrown 
 *  exception.
 */
#define CATCH_EXCEPTION                                      \
   catch(OFELIException &e) {                                \
      std::cout << "OFELI exception: " << e.what() << endl;  \
      return 1;                                              \
   }                                                         \
   catch(runtime_error &e) {                                 \
      std::cout << "Runtime exception: " << e.what() << endl;\
      return 1;                                              \
   }                                                         \
   catch( ... ) {                                            \
      std::cout << "Unexpected Exception: " << endl;         \
      return 1;                                              \
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
/*! \def MeshElementLoop(mesh,element)
 *  \ingroup Mesh
 * A macro to loop on mesh elements
 * \a mesh : Instance of Mesh
 * \a element : Pointer to pointed element
 */
#define MeshElementLoop(mesh,element) for ((mesh).topElement(); ((element)=(mesh).getElement());)

#define mesh_elements(mesh) for ((mesh).topElement(); ((the_element)=(mesh).getElement());)
#define mesh_nodes(mesh) for ((mesh).topNode(); ((the_node)=(mesh).getNode());)
#define mesh_sides(mesh) for ((mesh).topSide(); ((the_side)=(mesh).getSide());)
#define mesh_edges(mesh) for ((mesh).topEdge(); ((the_edge)=(mesh).getEdge());)
#define mesh_boundary_sides(mesh) for ((mesh).topBoundarySide(); ((the_side)=(mesh).getBoundarySide());)
#define The_element (*the_element)
#define The_node (*the_node)
#define The_side (*the_side)
#define The_edge (*the_edge)
#define element_label (the_element->n())
#define node_label (the_node->n())
#define side_label (the_side->n())
#define edge_label (the_edge->n())
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \def MeshElements(mesh)
 *  \ingroup Mesh
 * A macro to loop on mesh elements
 * <tt>mesh</tt>: Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theElement</tt> to current Element
 */
#define MeshElements(mesh) for ((mesh).topElement(); (theElement=(mesh).getElement());)

/*! \def MeshActiveElements(mesh)
 *  \ingroup Mesh
 * A macro to loop on mesh active elements
 * \a mesh : Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theElement</tt> to current Element
 * @note: This macro is necessary only if adaptive meshing is used
 */
#define MeshActiveElements(mesh) for ((mesh).topElement(); (theElement=(mesh).getActiveElement());)

/*! \def MeshNodeLoop(mesh,node)
 *  \ingroup Mesh
 * A macro to loop on mesh nodes
 * <tt>mesh</tt>: Instance of Mesh
 * <tt>node</tt>: Pointer to pointed node
 */
#define MeshNodeLoop(mesh,node) for ((mesh).topNode(); ((node)=(mesh).getNode());)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*! \def MeshBoundaryNodeLoop(mesh,node)
 *  \ingroup Mesh
 * A macro to loop on mesh nodes
 * <tt>mesh</tt>: Instance of Mesh
 * \a node</tt>: Pointer to pointed node
 *
 * @note: Boundary node list must has been created by Mesh::getBoundaryNodes
 */
#define MeshBoundaryNodeLoop(mesh,node) for ((mesh).topBoundaryNode(); ((node)=(mesh).getBoundaryNode());)
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \def MeshNodes(mesh)
 *  \ingroup Mesh
 * A macro to loop on mesh nodes
 * \a mesh : Instance of Mesh
 * @note: Each iteration updates the pointer \a theNode to current Node
 */
#define MeshNodes(mesh) for ((mesh).topNode(); (theNode=(mesh).getNode());)

/*! \def MeshBoundaryNodes(mesh)
 *  \ingroup Mesh
 * A macro to loop on mesh nodes
 * <tt>mesh</tt>: Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theNode</tt> to current Node
 */
#define MeshBoundaryNodes(mesh) for ((mesh).topBoundaryNode(); (theNode=(mesh).getBoundaryNode());)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*! \def MeshSideLoop(mesh,side)
 *  \ingroup Mesh
 * A macro to loop on mesh sides
 * <tt>mesh</tt>: Instance of Mesh
 * <tt>side</tt>: Pointer to pointed side
 */
#define MeshSideLoop(mesh,side) for ((mesh).topSide(); ((side)=(mesh).getSide());)
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \def MeshSides(mesh)
 *  \ingroup Mesh
 * A macro to loop on mesh sides
 * <tt>mesh</tt>: Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theSide</tt> to current Side
 */
#define MeshSides(mesh) for ((mesh).topSide(); (theSide=(mesh).getSide());)

/*! \def MeshSideSet(mesh)
 *  \ingroup Mesh
 * A macro to loop on a subset of mesh sides
 * <tt>sl</tt>: Instance of SideList class
 * @note: Each iteration updates the pointer <tt>theSide</tt> to current Side
 */
#define MeshSideSet(sl) for ((sl).top(); (theSide=(sl).get());)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*! \def MeshBoundarySideLoop(mesh,side)
 *  \ingroup Mesh
 * A macro to loop on mesh boundary sides
 * \a mesh : Instance of Mesh
 * \a side : Pointer to pointed boundary side
 * Important: List of boundary sides must have been previously created by using class SideList
 */
#define MeshBoundarySideLoop(mesh,side) for ((mesh).topBoundarySide(); ((side)=(mesh).getBoundarySide());)
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \def MeshBoundarySides(mesh)
 *  \ingroup Mesh
 * A macro to loop on mesh boundary sides
 * <tt>mesh</tt>: Instance of Mesh
 *
 * Notes:
 *   - List of boundary sides must have been previously created by using class SideList
 *   - Each iteration updates the pointer <tt>theSide</tt> to current Side
 */
#define MeshBoundarySides(mesh) for ((mesh).topBoundarySide(); (theSide=(mesh).getBoundarySide());)

/*! \def MeshEdges(mesh)
 *  \ingroup Mesh
 * A macro to loop on mesh edges
 * <tt>mesh</tt>: Instance of Mesh
 * @note: Each iteration updates the pointer <tt>theEdge</tt> to current Edge
 */
#define MeshEdges(mesh) for ((mesh).topEdge(); (theEdge=(mesh).getEdge());)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*! \def ElementNodeLoop
 *  \ingroup Mesh
 * A macro to loop on element nodes
 * @param el Instance of Element
 * @param nd Pointer to pointed node
 */
#define ElementNodeLoop(el,nd) for (el.topNode(); ((nd)=el.getNode()); )
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define MESH_EL MeshElements(*_theMesh)
#define MESH_ND MeshNodes(*_theMesh)
#define MESH_SD MeshSides(*_theMesh)
#define MESH_BD_SD MeshBoundarySides(*_theMesh)
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


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
          NbTimeSteps = int(theFinalTime/theTimeStep); \
          for (theTime=theTimeStep, theStep=1; theTime<theFinalTime+0.001*theTimeStep; theTime+=theTimeStep, ++theStep)

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
