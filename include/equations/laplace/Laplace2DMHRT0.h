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

                          Definition of class Laplace2DMHRT0
                for 2-D Laplace equation using 3-Node Mixed Hybrid Lowest 
                         degree Raviart-Thomas Finite element

  ==============================================================================*/


#ifndef __LAPLACE_2DMHRT0_H
#define __LAPLACE_2DMHRT0_H

#include "equations/laplace/Equa_Laplace.h"
#include "linear_algebra/Assembly.h"
#include "io/UserData.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {

/*! \file Laplace2DMHRT0.h
 *  \brief Definition file for class Laplace2DMHRT0.
 */

/*! \class Laplace2DMHRT0
 *  \ingroup Laplace
 *  \brief To build element equation for the 2-D elliptic equation
 *  using the Mixed Hybrid finite element at lowest degree (Raviart-Thomas <tt>RT<sub>0</sub></tt>).
 */

class Laplace2DMHRT0 : virtual public Equa_Laplace<real_t,3,3,2,2> {

 public:

/// \brief Default Constructor.
/// \details Constructs an empty equation.
    Laplace2DMHRT0() { }

/** \brief Constructor with problem data.
 *  @param [in] ms Mesh instance
 *  @param [in] A Problem matrix in Sparse format.
 *  This matrix must be zeroed before calling the constructor
 *  @param [in] b Problem right-hand side
 */
    Laplace2DMHRT0(Mesh&             ms,
                   SpMatrix<real_t>& A,
                   Vect<real_t>&     b);

/// \brief Destructor
    ~Laplace2DMHRT0() { }

/** \brief Define Diffusivity (or permeability) matrix.
 *  \details By default (if this function is not called) the identity matrix
 *  (Laplace equation) is used.
 *  @param [in] K Diffusivity matrix as \b LocalMatrix instance.
 *  Must be symmetric positive definite
 */
    void setDiffusivity(const LocalMatrix<real_t,2,2>& K);

/** \brief Build global matrix and right-hand side.
 *  \details The problem matrix and right-hand side are the ones used in the constructor.
 *  They are updated in this member function.
 */
    void build();

/** \brief Perform post calculations
 *  @param [in] lambda Solution (Lagrange multiplier) calculated at edges
 *  @param [in] f Vector containing the right-hand side of the Laplace equation
 *  @param [in] v Vector containing solution at mesh nodes
 *  @param [in] p Vector containing gradient at elements
 *  @param [in] u Vector containing solution at elements
 */
    void Post(const Vect<real_t>&         lambda,
              const Vect<real_t>&         f,
                    Vect<real_t>&         v, 
                    Vect<Point<real_t> >& p,
                    Vect<real_t>&         u);

/** \brief Solve the linear system of equations using the Conjugate Gradient iterative method.
 *  \details The matrix is preconditioned by an ILU method.
 *  @param [out] u Vector containing the solution at all sides (Sides where boundary conditions
 *  are prescribed are included).
 *  @return Number of performed iterations in the CG method. Note that the maximal number
 *  is <tt>1000</tt> and the tolerance is <tt>1.e-8</tt>
 */
    int solve(Vect<real_t>& u);

 private:

   SpMatrix<real_t>           *_A;
   LocalMatrix<real_t,2,2>    _K, _IK;
   Side                       *_sd1, *_sd2, *_sd3;
   LocalVect<Point<real_t>,3> _n, _ce;
   real_t                     _area;
   Point<real_t>              _c;
   Vect<size_t>               _vv;
   size_t NotOnSide(Element* el,
                    size_t   n1,
                    size_t   n2);
   void Set(const Element* el);
   void LM_LHS();
   void LM_RHS();
   void SolAtNodes(      Element*      el,
                   const Vect<real_t>& lambda,
                         Vect<real_t>& v);
   void getNormals(const Triang3 &tr);
   real_t getIP(const LocalMatrix<real_t,2,2>& A,
                const Point<real_t>&           a,
                const Point<real_t>&           b);
   Point<real_t> getP(const LocalMatrix<real_t,2,2>& A,
                      const Point<real_t>&           a);
};

} /* namespace OFELI */

#endif
