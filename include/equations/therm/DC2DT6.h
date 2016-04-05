/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

                            Definition of class DC2DT6
    for 2-D diffusion-convection scalar equation using 6-node triangular element

  ==============================================================================*/


#ifndef __DC2DT6_H
#define __DC2DT6_H


#include "equations/therm/Equa_Therm.h"
#include "io/UserData.h"
#include "shape_functions/Triang6S.h"
#include "shape_functions/Line3.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file DC2DT6.h
 *  \brief Definition file for class DC2DT6.
 */

/*! \class DC2DT6
 *  \ingroup Therm
 *  \brief Builds finite element arrays for thermal diffusion and convection in 2-D domains
 *  using 6-Node triangles.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that will be multiplied by the contribution of the current element.
 *  This makes possible testing different algorithms.
 *
 */

class DC2DT6 : virtual public Equa_Therm<real_t,6,6,3,3>
{

 public :

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    DC2DT6();

/// \brief Constructor for an element.
/// @param [in] el Pointer to element.
    DC2DT6(const Element* el);

/// \brief Constructor for a boundary side.
/// @param [in] sd Pointer to side.
    DC2DT6(const Side* sd);

/** \brief Constructor for an element (Transient case).
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 */
    DC2DT6(const Element*      el,
           const Vect<real_t>& u,
                 real_t        time=0.);

/** \brief Constructor for an element (transient case) with specification of time 
 *  integration scheme.
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>1</tt>]
 *  @param [in] deltat Value of time step
 *  @param [in] scheme Time Integration Scheme:
 *  <ul>
 *     <li><tt>FORWARD_EULER</tt>: Forward Euler scheme
 *     <li><tt>BACKWARD_EULER</tt>: Backward Euler scheme,
 *     <li><tt>CRANK_NICOLSON</tt>: Crank-Nicolson Euler scheme.
 *  </ul>
 */
    DC2DT6(const Element*      el,
           const Vect<real_t>& u,
                 real_t        time,
                 real_t        deltat,
                 int           scheme);

/** \brief Constructor for a boundary side (transient case).
 *  @param [in] sd Pointer to side.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 */
    DC2DT6(const Side*         sd,
           const Vect<real_t>& u,
                 real_t        time=0.);

/** \brief Constructor for a side (transient case) with specification of time integration scheme.
 *  @param [in] sd Pointer to side.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value
 *  @param [in] deltat Value of time step
 *  @param [in] scheme Time Integration Scheme: To be chosen among the enumerated values:
 *  <ul>
 *     <li><tt>FORWARD_EULER</tt>: Forward Euler scheme
 *     <li><tt>BACKWARD_EULER</tt>: Backward Euler scheme,
 *     <li><tt>CRANK_NICOLSON</tt>: Crank-Nicolson Euler scheme.
 *  </ul>
 */
    DC2DT6(const Side*         sd,
           const Vect<real_t>& u,
                 real_t        time,
                 real_t        deltat,
                 int           scheme);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Overall constructor
    DC2DT6(Mesh&              ms,
           SkSMatrix<real_t>& a,
           Vect<real_t>&      b);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Destructor
    ~DC2DT6();

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Build equation
    void build();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add lumped capacity matrix to left-hand side after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
//    void LCapacityToLHS(real_t coef=1);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add lumped capacity contribution to right-hand side after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
//    void LCapacityToRHS(real_t coef=1);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add lumped capacity contribution to left and right-hand sides after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
//    void LCapacity(real_t coef) { LCapacityToLHS(coef); LCapacityToRHS(coef); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add Consistent capacity matrix to left-hand side after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
//    void CapacityToLHS(real_t coef=1);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add Consistent capacity contribution to right-hand side after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
//    void CapacityToRHS(real_t coef=1);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add Consistent capacity contribution to left and right-hand sides after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
//    void Capacity(real_t coef=1) { CapacityToLHS(coef); CapacityToRHS(coef); }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Add diffusion matrix to left hand side after multiplying it by coefficient <tt>coef</tt>
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void Diffusion(real_t coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(real_t coef=1) { coef = 1; }

/** \brief Add convection matrix to left hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] v Constant velocity vector.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(Point<real_t>& v,
                    real_t         coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  \details Case where velocity field is given by a vector <tt>v</tt>
 *  @param [in] v Velocity vector.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(const Vect<real_t>& v,
                          real_t        coef=1);

/** \brief Add convection contribution to right-hand side after multiplying it by coefficient
 *  <tt>coef</tt>.
 *  \details To be used for explicit convection term.
 *  @param [in] v Velocity vector.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void _RHSConvection(const Point<real_t>& v,
                              real_t         coef=1);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Add body right-hand side term to right hand side after multiplying it by coefficient \c coef.
    void BodyRHS(UserData<real_t>& ud);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Add body right-hand side term to right hand side
 *  @param [in] b Local vector (of size <tt>6</tt>) containing source at element nodes
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size 6 or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BodyRHS(const Vect<real_t>& b,
                       int           opt=GLOBAL_ARRAY);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] b Local vector (of size \c 3) containing source at side nodes.
 */
    void BoundaryRHS(UserData<real_t>& ud);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] sf Vector containing source at side nodes
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size 3 or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BoundaryRHS(const Vect<real_t>& sf,
                           int           opt=GLOBAL_ARRAY);

 protected:
    void set(const Element *el);
    void set(const Side *sd);

 private:
   Point<real_t>  _x[6], _s[3];
   real_t         _a3;
   Triang6S       *_tr;
   Line3          *_ln;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
