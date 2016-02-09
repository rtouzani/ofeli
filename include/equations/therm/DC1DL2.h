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

                            Definition of class DC1DL2
       for 1-D diffusion-convection scalar equation using 2-node line element

  ==============================================================================*/


#ifndef __DC1DL2_H
#define __DC1DL2_H


#include "equations/therm/Equa_Therm.h"

namespace OFELI {

/*! \file DC1DL2.h
 *  \brief Definition file for class DC1DL2.
 */

/*! \class DC1DL2
 *  \ingroup Therm
 *  \brief Builds finite element arrays for thermal diffusion and convection in 1-D 
 *  using 2-Node elements.
 *
 *  \details Note that members calculating element arrays have as an argument a real
 *  <tt>coef</tt> that will be multiplied by the contribution of the current element. 
 *  This makes possible testing different algorithms.
 */

class DC1DL2 : public Equa_Therm<real_t,2,2,1,1>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    DC1DL2();

/// \brief Constructor for an element.
    DC1DL2(const Element* el);

/** \brief Constructor for an element (transient case).
    @param el [in] Pointer to element
    @param u [in] Vect instance that contains solution at previous time step
    @param time [in] Current time value (Default value is <tt>0</tt>)
*/
    DC1DL2(const Element*      el,
           const Vect<real_t>& u,
                 real_t        time=0.);

/** \brief Constructor for an element (transient case) with specification of time integration scheme.
 *  @param [in] el Pointer to element.
 *  @param [in] u Vect instance that contains solution at previous time step.
 *  @param [in] time Current time value (Default value is <tt>0</tt>).
 *  @param [in] deltat Value of time step
 *  @param [in] scheme Time Integration Scheme:
 *  <ul>
 *    <li> <tt>FORWARD_EULER</tt> for Forward Euler scheme
 *    <li> <tt>BACKWARD_EULER</tt> for Backward Euler scheme
 *    <li> <tt>CRANK_NICOLSON</tt> for Crank-Nicolson Euler scheme
 *  </ul>
 */
    DC1DL2(const Element*      el,
           const Vect<real_t>& u,
                 real_t        time,
                 real_t        deltat,
                 int           scheme);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// Constructors for a mesh
    DC1DL2(Mesh&         mesh,
           Vect<real_t>& b,
           real_t&       t,
           real_t&       ts)
          : Equation<real_t,2,2,1,1>(mesh,b,t,ts) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Destructor
    ~DC1DL2();

/// \brief Build the linear system without solving
    void build();

/** \brief Add lumped capacity matrix to left-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void LCapacityToLHS(real_t coef=1);

/** \brief Add lumped capacity contribution to right-hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void LCapacityToRHS(real_t coef=1);

/** \brief Add lumped capacity contribution to left and right-hand sides after multiplying 
 *  it by coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void LCapacity(real_t coef) { LCapacityToLHS(coef); LCapacityToRHS(coef); }

/** \brief Add Consistent capacity matrix to left-hand side after multiplying it by 
 *  coefficient \c coef.
 *  @param [in] coef Coefficient to multiply by added term [default: 1]
 */
    void CapacityToLHS(real_t coef=1);

/** \brief Add Consistent capacity contribution to right-hand side after multiplying it 
 *  by coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void CapacityToRHS(real_t coef=1);

/** \brief Add Consistent capacity contribution to left and right-hand sides after 
 *  multiplying it by coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Capacity(real_t coef=1);

/** \brief Add diffusion matrix to left hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Diffusion(real_t coef=1);

/** \brief Add diffusion contribution to right hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  \details To be used for explicit diffusion term
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void DiffusionToRHS(real_t coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] v Constant velocity vector
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Convection(const real_t& v,
                          real_t  coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  \details Case where velocity field is given by a vector <tt>v</tt>
 *  @param [in] v Velocity vector
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Convection(const Vect<real_t>& v,
                          real_t        coef=1);

/** \brief Add convection matrix to left-hand side after multiplying it by coefficient 
 *  <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void Convection(real_t coef=1);

/** \brief Add convection contribution to right-hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  \details To be used for explicit convection term.
 *  @param [in] v Velocity vector
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void ConvectionToRHS(const real_t& v,
                               real_t  coef=1);

/** \brief Add convection contribution to right-hand side after multiplying it by coefficient
 *  <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined 
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void ConvectionToRHS(real_t coef=1);

/** \brief Add body right-hand side term to right hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] ud Instance of UserData or of a derived class. Contains
 *  a member function that provides body source.
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void BodyRHS(UserData<real_t>& ud,
                 real_t            coef=1);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] b Vector containing source at element nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BodyRHS(const Vect<real_t>& b,
                       int           opt=GLOBAL_ARRAY);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  @param [in] ud Instance of \b UserData or of an inherited class. Contains
 *  a member function that provides body source.
 *  @param [in] coef Coefficient to multiply by added term [default: <tt>1</tt>]
 */
    void BoundaryRHS(UserData<real_t>& ud,
                     real_t            coef=1);

/// \brief Add boundary right-hand side flux to right hand side.
/// @param [in] flux Vector containing source at side nodes.
    void BoundaryRHS(real_t flux);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] b Vector containing source at side nodes.
 *  @param [in] opt Vector is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> or global
 *  (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes [Default: <tt>GLOBAL_ARRAY</tt>].
 */
    void BoundaryRHS(const Vect<real_t>& b,
                           int           opt=GLOBAL_ARRAY);

/// \brief Return (constant) heat flux in element.
    real_t Flux() const;

/** \brief Set equation input data
 *  @param [in] opt Parameter that selects data type for input. This parameter
 *  is to be chosen in the enumerated variable EqDataType
 *  <ul>
 *     <li><tt>INITIAL_FIELD</tt>: Initial temperature
 *     <li><tt>BOUNDARY_CONDITION_DATA</tt>: Boundary condition (Dirichlet)
 *     <li><tt>SOURCE_DATA</tt>: Heat source
 *     <li><tt>FLUX_DATA</tt>: Heat flux (Neumann boundary condition)
 *     <li><tt>VELOCITY</tt>: Velocity vector (for the convection term)
 *  </ul>
 *  @param [in] u Vector containing input data
 */
    void setInput(EqDataType    opt,
                  Vect<real_t>& u);

 private:
    TrMatrix<real_t> _A;
    void set(const Element* el);
};

} /* namespace OFELI */

#endif
