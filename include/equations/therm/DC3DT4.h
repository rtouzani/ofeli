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

                         Definition of class DC3DT4
   for 3-D diffusion-convection scalar equation using 4-node tetrahedral element

  ==============================================================================*/


#ifndef __DC3DT4_H
#define __DC3DT4_H


#include "equations/therm/Equa_Therm.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file DC3DT4.h
 *  \brief Definition file for class DC3DT4.
 */

/*! \class DC3DT4
 *  \ingroup Therm
 *  \brief Builds finite element arrays for thermal diffusion and convection in 3-D domains
 *  using 4-Node tetrahedra.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that will be multiplied by the contribution of the current element.
 *  This makes possible testing different algorithms.
 */


class DC3DT4 : public Equa_Therm<real_t,4,4,3,3>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    DC3DT4() { }

/// \brief Constructor for an element.
/// @param [in] el Pointer to element.
    DC3DT4(const Element* el);

/// \brief Constructor for a boundary side.
/// @param sd [in] Pointer to side.
    DC3DT4(const Side* sd);

/** \brief Constructor for an element (transient case).
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 */
    DC3DT4(const Element*      el,
           const Vect<real_t>& u,
           real_t              time=0.);

/** \brief Constructor for a boundary side (transient case).
 *  @param [in] sd Pointer to side.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value [Default: <tt>0</tt>].
 */
    DC3DT4(const Side*         sd,
           const Vect<real_t>& u,
           real_t              time=0.);

/** \brief Constructor for an element (transient case) with specification of time integration scheme.
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value.
 *  @param [in] deltat Value of time step
 *  @param [in] scheme Time Integration Scheme:
 *  <ul>
 *     <li><tt>FORWARD_EULER</tt>:  Forward Euler scheme
 *     <li><tt>BACKWARD_EULER</tt>: Backward Euler scheme
 *     <li><tt>CRANK_NICOLSON</tt>: Crank-Nicolson Euler scheme
 *  </ul>
 */
    DC3DT4(const Element*      el,
           const Vect<real_t>& u,
           real_t              time,
           real_t              deltat,
           int                 scheme);

/** \brief Constructor for a side (transient case) with specification of time integration scheme.
 *  @param [in] sd Pointer to side.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value.
 *  @param [in] deltat Value of time step
 *  @param [in] scheme Time Integration Scheme ():
 *  <ul>
 *     <li>FORWARD_EULER:  for Forward Euler scheme
 *     <li>BACKWARD_EULER: for Backward Euler scheme
 *     <li>CRANK_NICOLSON: for Crank-Nicolson Euler scheme
 *  </ul>
 */
    DC3DT4(const Side*         sd,
           const Vect<real_t>& u,
           real_t              time,
           real_t              deltat,
           int                 scheme);

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    DC3DT4(Mesh& ms);

/// \brief Destructor
    ~DC3DT4();

/// \brief Build the linear system without solving
    void build();

/** \brief Add lumped capacity contribution to left and right-hand sides after multiplying
 *  it by coefficient <tt>coef</tt>.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void LCapacity(real_t coef=1.) { LCapacityToLHS(coef); LCapacityToRHS(coef); }

/** \brief Add lumped capacity matrix to left-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void LCapacityToLHS(real_t coef=1);

/** \brief Add lumped capacity contribution to right-hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void LCapacityToRHS(real_t coef=1);

/** \brief Add Consistent capacity contribution to left and right-hand sides after 
 *  multiplying it by coefficient <tt>coef</tt>.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Capacity(real_t coef=1) { CapacityToLHS(coef); CapacityToRHS(coef); }

/** \brief Add consistent capacity matrix to left-hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void CapacityToLHS(real_t coef=1);

/** \brief Add consistent capacity contribution to right-hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void CapacityToRHS(real_t coef=1);

/// \brief Add diffusion matrix to left hand side after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
    void Diffusion(real_t coef=1);

/** \brief Add diffusion matrix to left hand side after multiplying it by coefficient <tt>coef</tt>
 *  \details Case where the diffusivity matrix is given as an argument.
 *  @param [in] diff Diffusion matrix (class DMatrix).
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Diffusion(const DMatrix<real_t>& diff,
                   real_t                 coef=1);

/** \brief Add diffusion contribution to right hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  \details To be used for explicit diffusion term
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void DiffusionToRHS(real_t coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(real_t coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] v Constant velocity vector
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(const Point<real_t>& v,
                    real_t               coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  \details Case where velocity field is given by a vector <tt>v</tt>.
 *  @param [in] v Velocity vector.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(const Vect<Point<real_t> >& v,
                    real_t                      coef=1);

/** \brief Add convection contribution to right-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  \details To be used for explicit convection term.
 *  @param [in] v Velocity vector.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void RHS_Convection(const Point<real_t>& v,
                        real_t               coef=1.);

/** \brief Add body right-hand side term to right hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] ud Instance of UserData or of an inherited class. Contains
 *  a member function that provides body source.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void BodyRHS(UserData<real_t>& ud,
                 real_t            coef=1);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] b Local vector containing source at element nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>4</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BodyRHS(const Vect<real_t>& b,
                 int                 opt=GLOBAL_ARRAY);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] ud Instance of UserData or of an inherited class. Contains
 *  a member function that provides body source.
 *  @param [in] coef Value by which the added term is multiplied [Default: <tt>1</tt>].
 */
    void BoundaryRHS(UserData<real_t>& ud,
                     real_t            coef=1);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  \details Case where body source is given by a vector
 *  @param [in] b Vector containing source at side nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BoundaryRHS(const Vect<real_t>& b,
                     int                 opt=GLOBAL_ARRAY);

/// \brief Add boundary right-hand side flux to right hand side.
/// @param [in] flux Vector containing source at side nodes.
    void BoundaryRHS(real_t flux);

/// \brief Return (constant) heat flux in element
    Point<real_t> Flux() const;

/// \brief Return gradient of vector <tt>u</tt> in element.
/// <tt>u</tt> is a local vector.
    Point<real_t> Grad(const Vect<real_t>& u) const;

/** \brief Add contribution of periodic boundary condition (by a penalty technique).
 *  \details Boundary nodes where periodic boundary conditions are to be imposed must
 *  have codes equal to <tt>PERIODIC_A</tt> on one side and <tt>PERIODIC_B</tt> on 
 *  the opposite side.
 *  @param [in] coef Value of penalty parameter [Default: <tt>1.e20</tt>].
 */
    void Periodic(real_t coef=1.e20);

 protected:
   void set(const Element *el);
   void set(const Side *sd);

 private:
   Point<real_t> _grad;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
