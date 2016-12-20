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

                      Definition of class Laplace2DFVT
         for 2-D Laplace equation using the classical Finite Volume scheme

  ==============================================================================*/


#ifndef __LAPLACE_2DFVT_H
#define __LAPLACE_2DFVT_H

#include <float.h>
#include <stdlib.h>
#include <math.h>

#include "equations/laplace/Equa_Laplace.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Laplace2DFVT.h
 *  \brief Definition file for class Laplace.
 */

/*! \class Laplace2DFVT
 *  \ingroup Laplace
 *  \brief To build and solve the Laplace equation using a standard Finite Volume method.
 */

class Element;
class Triang3;
class Line2;

class Laplace2DFVT : virtual public Equa_Laplace<real_t,3,3,2,2>
{

 public:

   using Equa_Laplace<real_t,3,3,2,2>::run;
   using Equa_Laplace<real_t,3,3,2,2>::build;

/** \brief Standard constructor
 *  @param [in] ms Mesh instance
 *  @param [in] b Vect instance that contains Right-hand side
 *  @param [in] u Vect instance that contains solution
 */
    Laplace2DFVT(Mesh&         ms,
                 Vect<real_t>& b,
                 Vect<real_t>& u);

/** \brief Standard constructor
 *  @param [in] ms Mesh instance.
 *  The mesh must have been assigned the attribute <tt>ELEMENT_DOF</tt> to say that 
 *  unknowns are supported by elements.
 *  @param [in] A Problem matrix to be stored in sparse format (class SpMatrix)
 *  @param [in] b Vect instance that contains Right-hand side
 */
    Laplace2DFVT(Mesh&             ms,
                 SpMatrix<real_t>& A,
                 Vect<real_t>&     b);

/// \brief Destructor
    ~Laplace2DFVT();

/** \brief Check whether triangles are Delaunay ones
 *  @param [in] verb Output (<tt>>0</tt>) or not (<tt>0</tt>) list of failing elements
 *  @return ret Number of Non Delaunay triangles
 */
    int checkDelaunay(int verb=0);

/// \brief Build the linear system of equations
    void build(const Vect<real_t>& f);

/// \brief Build and solve the linear system of equations
    int run(const Vect<real_t>& f);

/// \brief Calculate left-hand side
    void LHS(const Element* e1,
             const Element* e2);

/// \brief Add right-hand side Contribution
    void RHS(const Vect<real_t>& f);

 private:

   bool                _in_place;
   Point<real_t>       _xc1, _xc2;
   real_t              _m1, _m2, _ll, _dd;
   SpMatrix<real_t>    *_A;
   Vect<real_t>        *_b, *_u;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
