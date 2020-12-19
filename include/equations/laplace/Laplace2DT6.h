/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                         Definition of class Laplace2DT6
              for 2-D Laplace equation using 6-node triangular element

  ==============================================================================*/


#ifndef __LAPLACE_2DT6_H
#define __LAPLACE_2DT6_H

#include "equations/laplace/Equa_Laplace.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Laplace2DT6.h
 *  \brief Definition file for class Laplace2DT6.
 */

/*! \class Laplace2DT6
 *  \ingroup Laplace
 *  \brief To build element equation for the Laplace equation
 *  using the 2-D triangle element (<tt>P<sub>2</sub></tt>).
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Triang6S;
class Line3;
  
class Laplace2DT6 : virtual public Equa_Laplace<real_t,6,6,3,3> {

 public:

/// \brief Default constructor.
    Laplace2DT6();

/** \brief Constructor with mesh.
 *  @param [in] ms Mesh instance
 */
    Laplace2DT6(Mesh& ms);

/** \brief Constructor using mesh and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in] u Problem right-hand side
 */
    Laplace2DT6(Mesh&         ms,
                Vect<real_t>& u);

/// \brief Destructor
    ~Laplace2DT6() { }

/// \brief Add finite element matrix to left-hand side
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

 private:

   void set(const Element* el);
   void set(const Side* sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
