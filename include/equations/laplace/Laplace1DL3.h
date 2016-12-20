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
 */

class Laplace1DL3 : virtual public Equa_Laplace<real_t,3,3,1,1> {

public :

    using AbsEqua<real_t>::setInput;

/** Constructor using mesh instance and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance that contains, after execution of \b run() the solution
 */
    Laplace1DL3(Mesh&         ms,
                Vect<real_t>& u);

/// \brief Constructor for an element
    Laplace1DL3(Element* el);

/// \brief Destructor
    ~Laplace1DL3();

/// \brief Add finite element matrix to left hand side
/// @param [in] coef Value to multiply by the added matrix
    void Matrix(real_t coef=1.);

/// \brief Add Right-hand side contribution
/// @param [in] f Vector of right-hand side of the Poisson equation at nodes
    void BodyRHS(const Vect<real_t>& f);

/** \brief Add Neumann contribution to Right-Hand %Side
 *  @param [in] n Parameter to select equal to <tt>0</tt> if the condition is at the left
 *  end of the domain and different if it is at the right of it
 *  @param [in] p Value of flux to add
 *  @note This member function is to be invoked only for the first or last element
 */
    void BoundaryRHS(int    n,
                     real_t p);

/** \brief Set Traction data
 *  @param [in] f Value of traction (Neumann boundary condition)
 *  @param [in] lr Option to choose location of the traction (<tt>-1</tt>: Left end, 
 *  <tt>1</tt>: Right end)
 */
    void setTraction(real_t f,
                     int    lr);

/** Run solution procedure
 *  This function is to be called when the constructor \b Laplace1DL2(mesh,u)
 *  is used.
 *  @return return code for the solution of the linear system
 */
    int run();

private:
   LocalMatrix<real_t,3,3> _dSh;
   real_t _lsf, _rsf;
   BMatrix<real_t> _A;
   void set(Element *el);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
