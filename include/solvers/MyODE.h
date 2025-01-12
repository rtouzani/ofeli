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

                   Template Class 'MyODE' for user specified ode

  ==============================================================================*/


#ifndef __MY_ODE_H
#define __MY_ODE_H

#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"


namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file MyODE.h
 *  \brief Definition file for abstract class MyODE.
 */

/*! \class MyODE
 *  \ingroup Solver
 *  \brief Abstract class to define by user specified (system of) ODE
 *
 *  \details The user has to implement a class that inherits from the present one where
 *  the virtual functions are implemented.
 *
 *  \author Rachid Touzani
 *  \copyright GNU Lesser Public License
 */

class MyODE
{

 public:

/// \brief Default Constructor
    MyODE() { }

/// \brief Constructor using mesh instance
/// @param mesh Reference to Mesh instance
    MyODE(const Mesh& mesh)
    {
       _theMesh = &mesh;
    }

/// \brief Destructor
    virtual ~MyODE() { }

/** \brief Virtual member function to define function defining ODE in case of a single ODE
 *  @param [in] t Time at which the ode function is evaluated
 *  @param [in] y Unknown
 *  @return Value of function
 */
    virtual real_t Function(real_t        t,
                            const real_t& y)
    { return 0.; }

/** \brief Virtual member function to define function defining ODE
 *  @param [in] t Time at which the ode function is evaluated
 *  @param [in] y Vector of unknowns
 *  @param [in] i component of function to define [Default: <tt>1</tt>].
 *  @return Value of function
 *  @warning The component must not be larger than vector size
 */
    virtual real_t Function(real_t              t,
                            const Vect<real_t>& y,
                            int                 i=1)
    { return 0.; }

/** \brief Virtual member function to define derivative of function in case of a single ODE
 *  @param [in] t Time at which the derivative function is evaluated
 *  @param [in] y Unknown
 *  @return Value of derivative
 */
    virtual real_t Jacobian(real_t        t,
                            const real_t& y)
    { return 0.; }

/** \brief Virtual member function to define partial derivatives of function
 *  @param [in] t Time at which the jacobian function is evaluated
 *  @param [in] y Vector of unknowns
 *  @param [in] i Function component [Default: <tt>1</tt>]
 *  @param [in] j Index of partial derivative [Default: <tt>1</tt>]
 *  @return Value of partial derivative
 */
    virtual real_t Jacobian(real_t              t,
                            const Vect<real_t>& y,
                            int                 i=1,
                            int                 j=1)
    { return 0.; }

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   const Mesh *_theMesh;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
