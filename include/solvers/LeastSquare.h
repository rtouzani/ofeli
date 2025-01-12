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

                         Definition of class 'LeastSquare' for 
                    least square approximation of a set of points

  ==============================================================================*/

#ifndef __LEAST_SQUARE_H
#define __LEAST_SQUARE_H

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/DMatrix.h"
#include "io/Fct.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file LeastSquare.h
 *  \brief Definition file for class LeastSquare.
 */

/*! \class LeastSquare
 * \ingroup Solver
 * \brief To compute a least square approximation.
 *
 * This class enables using approximation methods to mathematically define
 * a geometry.
 * 
 * The algorithms used in this class are largely inspired from the book:
 * An Introduction to NURBS, by David F. Rogers. Copyright (C) 2000 David F. Rogers,
 */

class LeastSquare
{

 public:

/** 
 * \brief Default constructor
 * \details The functions set(...) must be used to define data before computing the approximation
 * using run()
 */
   LeastSquare();

/** 
 * \brief Constructor of least square approximation using given basis functions
 * \details The least square approximation defines the function:
 *     a[0]*f[0] + a[1]*f[1] + ... + a[N]*f[N] 
 * @param [in] f Vector of references to functions (class Fct) defining basis of least square 
 *               approximation
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a Vector containing solution: Coefficients of basis functions
 */
   LeastSquare(vector<Fct *>& f, const Vect<real_t>& x, const Vect<real_t>& y, Vect<real_t>& a);

/**
 * \brief Constructor with given least square matrix
 * \details Matrix entries is defined by \c B(k,i) as the k-th basis function evaluated
 * at \c x(i) 
 * @param [in] B Rectangle matrix of Least square approximation
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a Vector containing solution: Coefficients of basis functions
 */
   LeastSquare(DMatrix<real_t>& B, const Vect<real_t>& y, Vect<real_t>& a);

/**
 * \brief Constructor for polynomial regression
 * \details Define N-degree polynomial to approximate in the sense of least squares
 *  of a given set of points in the plane. The resulting line has the equation:
 *  y = a[0] + a[1]*x + ... + a[N]*x^N
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [in] N Degree of approximation polynomial
 * @param [out] a Vector of coefficients of polynomial as defined here above
 */
   LeastSquare(const Vect<real_t>& x, const Vect<real_t>& y, size_t N, Vect<real_t>& a);

/**
 * \brief Constructor for linear regression
 * \details Define 1-degree polynomial to approximate in the sense of least squares
 * of a given set of points in the plane. The resulting line has the equation:
  * y = a0 + a1*x
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a0 Coefficient of constant term
 * @param [out] a1 Coefficient of first degree term
 */
   LeastSquare(const Vect<real_t>& x, const Vect<real_t>& y, real_t& a0, real_t& a1);

/// \brief Destructor
   ~LeastSquare();

/** \brief Set least square approximation using given basis functions
 * \details The least square approximation defines the function:
 * a[0]*f[0] + a[1]*f[1] + ... + a[N]*f[N] 
 * @param [in] f Vector of references to functions (class Fct) defining basis of least square 
 *               approximation
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a Vector containing solution: Coefficients of basis functions
 */
   void set(vector<Fct *>& f, const Vect<real_t>& x, const Vect<real_t>& y, Vect<real_t>& a);

/** \brief Set least square approximation using least square matrix
 *  \details Matrix entries is defined by \c B(k,i) as the k-th basis function evaluated
 * at \c x(i) 
 * @param [in] B Rectangle matrix of Least square approximation
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a Vector containing solution: Coefficients of basis functions
 */
   void set(DMatrix<real_t>& B, const Vect<real_t>& y, Vect<real_t>& a);

/** \brief Set least square approximation by polynomial regression
 * \details Define N-degree polynomial to approximate in the sense of least squares
 * of a given set of points in the plane. The resulting line has the equation:
  * y = a[0] + a[1]*x + ... + a[N]*x^N
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [in] N Degree of approximation polynomial
 * @param [out] a Vector of coefficients of polynomial as defined here above
 */
   void set(const Vect<real_t>& x, const Vect<real_t>& y, size_t N, Vect<real_t>& a);

/** \brief Set least square approximation by linear regression
 * \details Define 1-degree polynomial to approximate in the sense of least squares
 * of a given set of points in the plane. The resulting line has the equation:
  * y = a0 + a1*x
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a0 Coefficient of constant term
 * @param [out] a1 Coefficient of first degree term
 */
   void set(const Vect<real_t>& x, const Vect<real_t>& y, real_t& a0, real_t& a1);

/// \brief Compute least square approximation
   int run();

 private:

   int _opt, _alloc;
   size_t _np, _dim;
   vector<Fct *> *_fct;
   const Vect<real_t> *_x, *_y;
   Vect<real_t> _b, *_a;
   real_t *_a0, *_a1;
   DMatrix<real_t> _A, *_B;
};

} /* namespace OFELI */

#endif
