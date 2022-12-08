/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                      Definition of class Laplace1DL3
              for 1-D Laplace equation using 3-node line element (P2)

  ==============================================================================*/


#ifndef __LAPLACE_1DL3_H
#define __LAPLACE_1DL3_H

#include <float.h>
#include <stdlib.h>
#include <math.h>

#include "OFELI_Config.h"
#include "equations/laplace/Equa_Laplace.h"
#include "mesh/Element.h"
#include "linear_algebra/BMatrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Laplace1DL3.h
 *  \brief Definition file for class Laplace1DL3.
 */

/*! \class Laplace1DL3
 *  \ingroup Laplace
 *  \brief To build element equation for the 1-D elliptic equation
 *  using the 3-Node line (<tt>P<sub>2</sub></tt>).
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Laplace1DL3 : virtual public Equa_Laplace<3,3,1,1> {

public :

    using Equa::setInput;

/// \brief Default constructor. Initializes an empty equation
    Laplace1DL3();

/// \brief Constructor using mesh instance
/// @param [in] ms Mesh instance
    Laplace1DL3(Mesh& ms);

/** Constructor using mesh instance and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance that contains, after execution of \b run() the solution
 */
    Laplace1DL3(Mesh&         ms,
                Vect<real_t>& u);

/// \brief Destructor
    ~Laplace1DL3();

/// \brief Compute element matrix
    void LHS();

/// \brief Add Right-hand side contribution
/// @param [in] f Vector of right-hand side of the Poisson equation at nodes
    void BodyRHS(const Vect<real_t>& f);

/** \brief Add Neumann contribution to Right-Hand %Side
 *  @param [in] h Vector with size the total number of nodes. The first entry stands for
 *  the force at the first node (Neumann condition) and the last entry is the force at
 *  the last node (Neumann condition)
 */
    void BoundaryRHS(const Vect<real_t>& h);

/** \brief Set Traction data
 *  @param [in] f Value of traction (Neumann boundary condition)
 *  @param [in] lr Option to choose location of the traction (<tt>-1</tt>: Left end, 
 *  <tt>1</tt>: Right end)
 */
    void setTraction(real_t f,
                     int    lr);

private:
   LocalVect<Point<real_t>,3> _dSh[3];
   real_t _lsf, _rsf;
   BMatrix<real_t> _A;
   void set(const Element* el);
   void set(const Side* sd) { }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
