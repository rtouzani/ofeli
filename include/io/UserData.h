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

              Template Class 'UserData' for user data prescription

  ==============================================================================*/


#ifndef __USERDATA_H
#define __USERDATA_H

#include <float.h>
#include <stdlib.h>
#include <math.h>

#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Point.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file UserData.h
 *  \brief Definition file for abstract class UserData.
 */

/*! \class UserData
 *  \ingroup Util
 *  \brief Abstract class to define by user various problem data.
 *
 * The user has to implement a class that inherits from the present one where
 * the virtual functions are implemented.
 *
 * \tparam <T_> Data type (real_t, float, complex<real_t>, ...)
 */

template <class T_> class UserData
{

 public:

/// \brief Default Constructor
    UserData()
    {
       _time = 0.;
    }

/// \brief Constructor using mesh instance
/// @param mesh Reference to Mesh instance
    UserData(const class Mesh& mesh)
    {
       _theMesh = &mesh;
       _time = 0.;
    }

/// \brief Destructor
    virtual ~UserData() { }

/// \brief Set time value
    void setTime(real_t time) { _time = time; }

/** \brief Set Dirichlet Boundary Conditions
 *  \details This function loops over all nodes and calls for each node the member
 *  function BoundaryCondition to assign the value defined by it
 *  @param [out] b Vector that contains boundary conditions at nodes
 *  This vector must be sized before invoking this function
 */
    void setDBC(Vect<T_>& b)
    {
       size_t k = 1;
       MeshNodes(*_theMesh) {
          for (size_t i=1; i<=theNode->getNbDOF(); i++)
             b(k++) = BoundaryCondition(theNode->getCoord(),theNode->getCode(i),_time,i);
       }
    }

/** \brief Set initial data
 *  \details This function loops over all nodes and calls for each node the member
 *  function InitialData to assign the value defined by it
 *  @param [out] b Vector that contains initial data at nodes
 *  This vector must be sized before invoking this function
 */
    void setInitialData(Vect<T_>& b)
    {
       MeshNodes(*_theMesh) {
          for (size_t j=1; j<=theNode->getNbDOF(); j++)
             b(theNode->getDOF(j)) = InitialData(theNode->getCoord(),j);
       }
    }

/** \brief Set Nodewise Body Force using a Vect instance
 *  @param [in] b Vector containing body forces at nodes to impose
 *  \details This function loops over all nodes and calls for each node the member
 *  function BodyForce to assign the value defined by it
 *  @param [out] b Vector that contains body forces at nodes
 *  This vector must be sized before invoking this function
 */
    void setBodyForce(Vect<T_>& b)
    {
       MeshNodes(*_theMesh) {
          for (size_t j=1; j<=theNode->getNbDOF(); j++)
             b(TheNode.getDOF(j)) = BodyForce(TheNode.getCoord(),_time,j);
       }
    }

/** \brief Set Surface Force
 *  @param [in] b Vector containing surface forces at nodes to impose
 *  \details This function loops over all nodes and calls for each node the member
 *  function SurfaceForce to assign the value defined by it
 *  @param [out] b Vector that contains body forces at nodes
 *  This vector must be sized before invoking this function
 */
    void setSurfaceForce(Vect<T_>& b)
    {
       MeshNodes(*_theMesh) {
          for (size_t j=1; j<=theNode->getNbDOF(); j++)
             b(TheNode.getDOF(j)) = SurfaceForce(TheNode.getCoord(),TheNode.getCode(j),_time,j);
       }
    }

/** \brief Define boundary condition to impose at point of coordinates <tt>x</tt>, with code
 *  <tt>code</tt> at time <tt>time</tt>, for DOF <tt>dof</tt>
 *  \details Function to implement by user
 *  @param [in] x Coordinates of point at which the value is to be prescribed
 *  @param [in] code Code of node for which the value is to be prescribed
 *  @param [in] time Value of time [Default: <tt>0.</tt>]
 *  @param [in] dof Corresponding degree of freedom [Default: <tt>1</tt>]
 *  @return Value of boundary condition to prescribe corresponding to these parameters
 */
    virtual T_ BoundaryCondition(Point<real_t> x,
                                 int           code,
                                 real_t        time=0.,
                                 size_t        dof=1)
    {
       x = 0;
       time = 0;
       dof = 1;
       code = 0;
       return T_(0);
    }

/** \brief Define body force to impose at point of coordinates <tt>x</tt>, with code 
 *  <tt>code</tt> at time <tt>time</tt>, for DOF <tt>dof</tt>
 *  \details Function to implement by user
 *  @param [in] x Coordinates of point at which the body force is given
 *  @param [in] time Value of time [Default: <tt>0.</tt>]
 *  @param [in] dof Corresponding degree of freedom [Default: <tt>1</tt>]
 *  @return Value of body force corresponding to these parameters
 */
    virtual T_ BodyForce(Point<real_t> x,
                         real_t        time=0.,
                         size_t        dof=1)
    {
       x = 0;
       time = 0;
       dof = 1;
       return T_(0);
    }

/** \brief Define surface force to impose at point of coordinates <tt>x</tt>, with 
 *  code <tt>code</tt> at time <tt>time</tt>, for DOF <tt>dof</tt>
 *  \details Function to implement by user
 *  @param [in] x Coordinates of point at which the surface force is given
 *  @param [in] code Code of node for which the surface force is given
 *  @param [in] time Value of time [Default: <tt>0.</tt>]
 *  @param [in] dof Corresponding degree of freedom [Default: <tt>1</tt>]
 *  @return Value of surface force corresponding to these parameters
 */
    virtual T_ SurfaceForce(Point<real_t> x,
                            int           code,
                            real_t        time=0.,
                            size_t        dof=1)
    {
       x = 0;
       time = 0;
       dof = 1;
       code = 0;
       return T_(0);
    }

/** \brief Define initial data to impose at point of coordinates <tt>x</tt>, for DOF <tt>dof</tt>
 *  \details Function to implement by user
 *  @param [in] x Coordinates of point at which the surface force is given
 *  @param [in] dof Corresponding degree of freedom [Default: <tt>1</tt>]
 *  @return Value of initial data corresponding to these parameters
 */
    virtual T_ InitialData(Point<real_t> x,
                           size_t        dof=1)
    {
       x = 0;
       dof = 1;
       return T_(0.);
    }

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   real_t      _time;
   const Mesh *_theMesh;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
