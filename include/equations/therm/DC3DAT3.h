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

                           Definition of class DC3DAT3
       for 3-D diffusion-convection scalar equation in axisymmetric geometry
                         using 3-node triangular element

  ==============================================================================*/


#ifndef __DC3DAT3_H
#define __DC3DAT3_H


#include "equations/therm/Equa_Therm.h"
#include "linear_algebra/LocalMatrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file DC3DAT3.h
 *  \brief Definition file for class DC3DAT3.
 */

/*! \class DC3DAT3
 *  \ingroup Therm
 *  \brief Builds finite element arrays for thermal diffusion and convection in 3-D domains 
 *  with axisymmetry using 3-Node triangles.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  \c coef that will be multiplied by the contribution of the current element. 
 *  This makes possible testing different algorithms.
 *
 */


class DC3DAT3 : virtual public Equa_Therm<real_t,3,3,2,2>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    DC3DAT3();

/// \brief Constructor for an element.
/// @param [in] el Pointer to element.
    DC3DAT3(const Element* el);

/// \brief Constructor for a boundary side.
/// @param [in] sd Pointer to side.
    DC3DAT3(const Side* sd);

/** \brief Constructor for an element (transient case).
 *  @param [in] el Pointer to element
 *  @param [in] u Vect instance that contains solution at previous time step
 *  @param [in] time Current time value [Default: <tt>0</tt>]
 */
    DC3DAT3(const Element*      el,
            const Vect<real_t>& u,
                  real_t        time=0.);

/** \brief Constructor for an element (transient case) with specification of time 
 *  integration scheme.
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value.
 *  @param [in] deltat Value of time step
 *  @param [in] scheme Time Integration Scheme ():
 *  <ul>
 *     <li><tt>FORWARD_EULER</tt> for Forward Euler scheme, <tt>BACKWARD_EULER</tt> for
 *         Backward Euler scheme,
 *     <li><tt>CRANK_NICOLSON</tt> for Crank-Nicolson Euler scheme.
 *  </ul>
 */
    DC3DAT3(const Element*      el,
            const Vect<real_t>& u,
                  real_t        time,
                  real_t        deltat,
                  int           scheme);

/** \brief Constructor for a boundary side (transient case).
 *  @param [in] sd Pointer to side
 *  @param [in] u Vect instance that contains solution at previous time step
 *  @param [in] time Current time value [Default: <tt>0</tt>]
 */
    DC3DAT3(const Side*         sd,
            const Vect<real_t>& u,
                  real_t        time=0.);

/** \brief Constructor for a side (transient case) with specification of time integration scheme.
 *  @param [in] sd Pointer to side
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value.
 *  @param [in] deltat Value of time step
 *  @param [in] scheme Time Integration Scheme (enumerated values) :
 *  <ul>
 *     <li><tt>FORWARD_EULER</tt>: Forward Euler scheme
 *     <li><tt>BACKWARD_EULER</tt>: Backward Euler scheme
 *     <li><tt>CRANK_NICOLSON</tt>: Crank-Nicolson Euler scheme
 *  </ul>
 */
    DC3DAT3(const Side*         sd,
            const Vect<real_t>& u,
                  real_t        time,
                  real_t        deltat,
                  int           scheme);

/// \brief Destructor
    ~DC3DAT3();

/// \brief Add lumped capacity matrix to left-hand side after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LCapacityToLHS(real_t coef=1);

/// \brief Add lumped capacity contribution to right-hand side after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
    void LCapacityToRHS(real_t coef=1);

/** \brief Add lumped capacity contribution to left and right-hand sides after multiplying
 *  it by coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void LCapacity(real_t coef=1) { LCapacityToLHS(coef); LCapacityToRHS(coef); }

/** \brief Add Consistent capacity matrix to left-hand side after multiplying it by
 *  coefficient<tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void CapacityToLHS(real_t coef=1);

/** \brief Add Consistent capacity contribution to right-hand side after multiplying it
 *  by coefficient <tt>coef</tt>.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void CapacityToRHS(real_t coef=1);

/** \brief Add Consistent capacity contribution to left and right-hand sides after 
 *  multiplying it by coefficient <tt>coef</tt>.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Capacity(real_t coef=1) { CapacityToLHS(coef); CapacityToRHS(coef); }

/** \brief Add diffusion matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Diffusion(real_t coef=1);

/** \brief Add diffusion matrix to left-hand side after multiplying it by coefficient <tt>coef</tt>
 *  \details Case where the diffusivity matrix is given as an argument
 *  @param [in] diff Instance of class DMatrix containing diffusivity matrix
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
    void Diffusion(const LocalMatrix<real_t,2,2>& diff,
                         real_t                   coef=1);

/** \brief Add diffusion contribution to right-hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  \details To be used for explicit diffusion term
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>]
 */
   void DiffusionToRHS(real_t coef=1);

/** \brief Add body right-hand side term to right-hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] ud Instance of UserData or of an inherited class. Contains
 *  a member function that provides body source.
 */
    void BodyRHS(UserData<real_t>& ud);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] b Local vector (of size <tt>3</tt>) containing source at element nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BodyRHS(const Vect<real_t>& b,
                       int           opt=GLOBAL_ARRAY);

/// \brief Add boundary right-hand side term to right hand side.
/// @param [in] flux Value of flux to impose on the side
    void BoundaryRHS(real_t flux);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] sf Vector containing source at side nodes
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size 2 or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BoundaryRHS(const Vect<real_t>& sf,
                           int           opt=GLOBAL_ARRAY);

/// \brief Return gradient of a vector in element.
/// @param [in] u Vector for which gradient is computed.
    Point<real_t> & Grad(const Vect<real_t>& u);

/// \brief Build the linear system without solving
    void build();

 protected:
    void set(const Element *el);
    void set(const Side *sd);

 private:
    real_t        _r[3], _h;
    Point<real_t> _grad;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
