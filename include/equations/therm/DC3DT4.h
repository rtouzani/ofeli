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


class DC3DT4 : public Equa_Therm<4,4,3,3>
{

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    DC3DT4();

/// \brief Constructor using Mesh data
/// @param [in] ms Mesh instance
    DC3DT4(Mesh& ms);

/** \brief Constructor using Mesh and initial condition
 *  @param [in] ms Mesh instance
 *  @param [in] u Vect instance containing initial solution
 */
    DC3DT4(Mesh&         ms,
           Vect<real_t>& u);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Constructors for a mesh
    DC3DT4(Mesh&         mesh,
           Vect<real_t>& b,
           real_t        init_time,
           real_t        final_time,
           real_t        time_step)
   : Equation<4,4,3,3>(mesh,b,init_time,final_time,time_step) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Destructor
    ~DC3DT4();

/** \brief Add lumped capacity matrix to element matrix after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void LCapacity(real_t coef=1);

/** \brief Add consistent capacity matrix to element matrix after multiplying it by 
 *  coefficient <tt>coef</tt>
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Capacity(real_t coef=1);

/// \brief Add diffusion matrix to element matrix after multiplying it by coefficient \c coef.
/// @param [in] coef Coefficient to multiply by added term (default value = 1).
    void Diffusion(real_t coef=1);

/** \brief Add diffusion matrix to element matrix after multiplying it by coefficient <tt>coef</tt>
 *  \details Case where the diffusivity matrix is given as an argument.
 *  @param [in] diff Diffusion matrix (class DMatrix).
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Diffusion(const DMatrix<real_t>& diff,
                   real_t                 coef=1);

/** \brief Add convection matrix to element matrix after multiplying it by coefficient
 *  <tt>coef</tt>
 *  \details Case where velocity field has been previouly defined
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(real_t coef=1);

/** \brief Add convection matrix to element matrix after multiplying it by coefficient
 *  <tt>coef</tt>
 *  @param [in] v Constant velocity vector
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(const Point<real_t>& v,
                    real_t               coef=1);

/** \brief Add convection matrix to element matrix after multiplying it by coefficient
 *  <tt>coef</tt>
 *  \details Case where velocity field is given by a vector <tt>v</tt>.
 *  @param [in] v Velocity vector.
 *  @param [in] coef Coefficient to multiply by added term [Default: <tt>1</tt>].
 */
    void Convection(const Vect<Point<real_t> >& v,
                    real_t                      coef=1);

/** \brief Add body right-hand side term to right hand side.
 *  @param [in] f Vector containing source at nodes.
 */
    void BodyRHS(const Vect<real_t>& f);

/** \brief Add boundary right-hand side term to right hand side after multiplying it by
 *  coefficient <tt>coef</tt>
 *  \details Case where body source is given by a vector
 *  @param [in] f Vector containing source at nodes.
 */
    void BoundaryRHS(const Vect<real_t>& f);

/// \brief Add boundary right-hand side flux to right hand side.
/// @param [in] flux Vector containing source at side nodes.
    void BoundaryRHS(real_t flux);

/// \brief Return (constant) heat flux in element
    Point<real_t> Flux() const;

/// \brief Compute gradient of solution.
/// @param [in] g Elementwise vector containing gradient of solution.
    void Grad(Vect<Point<real_t> >& g);

/** \brief Add contribution of periodic boundary condition (by a penalty technique).
 *  \details Boundary nodes where periodic boundary conditions are to be imposed must
 *  have codes equal to <tt>PERIODIC_A</tt> on one side and <tt>PERIODIC_B</tt> on 
 *  the opposite side.
 *  @param [in] coef Value of penalty parameter [Default: <tt>1.e20</tt>].
 */
    void Periodic(real_t coef=1.e20);

 private:
   void set(const Element* el);
   void set(const Side* sd);
   Point<real_t> _grad;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
