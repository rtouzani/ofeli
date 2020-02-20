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

                         Definition of class Laplace3DT4
              for 3-D Laplace equation using 4-node tetrahedral element

  ==============================================================================*/


#ifndef __LAPLACE_3DT4_H
#define __LAPLACE_3DT4_H

#include "equations/laplace/Equa_Laplace.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Laplace3DT4.h
 *  \brief Definition file for class Laplace3DT4.
 */

/*! \class Laplace2DT3
 *  \ingroup Laplace
 *  \brief To build element equation for the Laplace equation
 *  using the 3-D tetrahedral element (<tt>P<sub>1</sub></tt>).
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Laplace3DT4 : virtual public Equa_Laplace<real_t,4,4,3,3> {

 public:

/// \brief Default constructor. Initializes an empty equation
    Laplace3DT4();

/** \brief Constructor with mesh.
 *  @param [in] ms Mesh instance
 */
    Laplace3DT4(Mesh& ms);

/** \brief Constructor using mesh and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in] u Reference to solution vector instance
 */
    Laplace3DT4(Mesh&         ms,
                Vect<real_t>& u);

/// \brief Destructor
    ~Laplace3DT4() { }

/// \brief Add finite element matrix to left-hand side
/// @param [in] coef Value to multiply by the added matrix
    void LHS();

/// \brief Add body source term to right-hand side
/// @param [in] f Vector containing the source given function at mesh nodes
    void BodyRHS(const Vect<real_t>& f);

/// \brief Add boundary source term to right-hand side
/// @param [in] h Vector containing the source given function at mesh nodes
    void BoundaryRHS(const Vect<real_t>& h);

/** \brief Build global stiffness and mass matrices for the eigen system
 *  @param [in] opt Flag to choose a lumed mass matrix (0) or consistent (1) [Default: <tt>0</tt>]
 */
    void buildEigen(int opt=0);

/** \brief Perform post calculations
 *  @param [in] u Solution at nodes
 *  @param [out] p Vector containing gradient at elements
 */
    void Post(const Vect<real_t>&   u,
              Vect<Point<real_t> >& p);

 private:
   void set(const Element* el);
   void set(const Side* sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
