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

                      Definition of class Laplace1DL2
              for 1-D Laplace equation using 2-node line element (P1)

  ==============================================================================*/


#ifndef __LAPLACE_1DL2_H
#define __LAPLACE_1DL2_H

#include <float.h>
#include <stdlib.h>
#include <math.h>

#include "OFELI_Config.h"
#include "equations/laplace/Equa_Laplace.h"
#include "mesh/Element.h"

namespace OFELI {

/*! \file Laplace1DL2.h
 *  \brief Definition file for class Laplace1DL2.
 */

/*! \class Laplace1DL2
 *  \ingroup Laplace
 *  \brief To build element equation for a 1-D elliptic equation
 *  using the 2-Node line element (<tt>P<sub>1</sub></tt>).
 */

class Laplace1DL2 : virtual public Equa_Laplace<real_t,2,2,1,1> {

public :

    using AbsEqua<real_t>::setInput;

/// \brief Constructor for an element
    Laplace1DL2(Element* el);

/** Constructor using mesh instance and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance that contains, after execution of \b run() the solution
 */
    Laplace1DL2(Mesh&         ms,
                Vect<real_t>& u);

/// \brief Destructor
    ~Laplace1DL2() { }

/// \brief Add finite element matrix to left hand side
/// @param [in] coef Value to multiply by the added matrix
    void Matrix(real_t coef=1.);

/// \brief Add Right-Hand %Side Contribution
/// @param [in] f Vector containing the source given function at mesh nodes
    void BodyRHS(const Vect<real_t>& f);

/** \brief Add Neumann contribution to Right-Hand %Side
 *  @param [in] n Parameter to select equal to <tt>0</tt> if the condition is at the left
 *  end of the domain and different if it is at the right of it
 *  @param [in] p Value of flux to add
 *  @note This member function is to be called only for the first or last element
 */
    void BoundaryRHS(int    n,
                     real_t p);

/** \brief Set Dirichlet boundary data
 *  @param [in] f Value to assign
 *  @param [in] lr Option to choose location of the value (<tt>-1</tt>: Left end,
 *  <tt>1</tt>: Right end)
 */
    void setBoundaryCondition(real_t f,
                              int    lr);

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
    real_t _lsf, _rsf, _lbc, _rbc;
    bool   _is_lbc, _is_rbc;
    TrMatrix<real_t> _A;
    void set(const Element* el);
};

} /* namespace OFELI */

#endif
