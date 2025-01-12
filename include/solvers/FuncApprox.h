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

                    Definition of class 'FuncApprox' for various
                          methods of function approximation

  ==============================================================================*/

#ifndef __FUNC_APPROX_H
#define __FUNC_APPROX_H

#include "OFELI_Config.h"
#include "linear_algebra/DMatrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file FuncApprox.h
 *  \brief Definition file for class FuncApprox.
 */

/*! \class FuncApprox
 * \ingroup Solver
 * \brief To set function approximation methods.
 *
 * \details This class enables using approximation methods for functions.
 * 
 * The algorithms used in this class are largely inspired from the book:
 * An Introduction to NURBS, by David F. Rogers. Copyright (C) 2000 David F. Rogers,
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define FUNC_APPROX_THRESHOLD 5.0e-6
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

class FuncApprox
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
enum Approx {
   LAGRANGE,
   LEAST_SQUARE,
   BSPLINE,
   BSPLINE_SURFACE,
   BEZIER,
   BEZIER_SURFACE,
   NURBS,
   NURBS_SURFACE
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

 public:

/// \brief Default constructor
/// \details The function setData can then be used
    FuncApprox();

/// \brief Destructor
    ~FuncApprox();

/** 
 * \brief Define Lagrange interpolation
 * \details This member function defines Lagrange interpolation data
 * @param [in] n Degree of interpolation polynomial (must be >= 1)
 * @param [in] x Abcissa of defining points
 * @param [in] y Values of points to interpolate
 * @param [out] f Function that will contain Lagrange interpolation polynomial once the
 *                function run() is invoked
 */
   void setLagrange(int n, const Vect<real_t>& x, const Vect<real_t>& y, Fct& f);

/** 
 * \brief Define least square approximation using given basis functions
 * \details The least square approximation defines the function:
 *     a[0]*f[0] + a[1]*f[1] + ... + a[N]*f[N] 
 * @param [in] f Vector of references to functions (class Fct) defining basis of least square 
 *               approximation
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a Vector containing solution: Coefficients of basis functions
 */
   void setLeastSquare(const vector<Fct *>& f, const Vect<real_t>& x, const Vect<real_t>& y, Vect<real_t>& a);

/**
 * \brief Define least square approximation with given least square matrix
 * \details Matrix entries is defined by \c B(k,i) as the k-th basis function evaluated
 * at \c x(i) 
 * @param [in] B Rectangle matrix of Least square approximation
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a Vector containing solution: Coefficients of basis functions
 */
   void setLeastSquare(DMatrix<real_t>& B, const Vect<real_t>& y, Vect<real_t>& a);

/**
 * \brief Define least square approximation using polynomial regression
 * \details Define N-degree polynomial to approximate in the sense of least squares
 *  of a given set of points in the plane. The resulting line has the equation:
 *  y = a[0] + a[1]*x + ... + a[N]*x^N
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [in] N Degree of approximation polynomial
 * @param [out] a Vector of coefficients of polynomial as defined here above
 */
   void setLeastSquare(const Vect<real_t>& x, const Vect<real_t>& y, size_t N, Vect<real_t>& a);

/**
 * \brief Constructor foDefine least square approximation using linear regression
 * \details Define 1-degree polynomial to approximate in the sense of least squares
 * of a given set of points in the plane. The resulting line has the equation:
  * y = a0 + a1*x
 * @param [in] x Vector defining x-coordinates of points
 * @param [in] y Vector defining y-coordinates of points
 * @param [out] a0 Coefficient of constant term
 * @param [out] a1 Coefficient of first degree term
 */
   void setLeastSquare(const Vect<real_t>& x, const Vect<real_t>& y, real_t& a0, real_t& a1);

/** 
 * \brief Compute the resulting least square fitting function
 * \details Once the function run() is invoked, this member function gives the analytical expression
 * of the resulting least square fitting function
 * @param [in] f Reference to approximating function using the least square method
 */
   void getLeastSquare(Fct& f);

/** 
 * \brief Define BSpline approximation
 * \details We set here BSpline approximation data. No computation are done yet.
 * @param [in] n Number of defining polygon vertices
 * @param [in] c Order of the B-spline basis function
 * @param [in] np Number of points to be calculated on the curve
 * @param [in] b Vector containing the defining polygon vertices
 * @param [out] p Vector containing the resulting curve points
 */
   void setBSpline(size_t n, size_t c, size_t np, const Vect<real_t>& b, Vect<real_t>& p);

/**
 * \brief Define BSpline surface modelling
 * @param [in] m One less than the number of polygon vertices in w direction
 * @param [in] n One less than the number of polygon vertices in u direction
 * @param [in] c Order of the B-spline basis function in w direction
 * @param [in] d Order of the B-spline basis function in u direction
 * @param [in] npu Number of parametric lines in the u direction
 * @param [in] npw Number of parametric lines in the w direction
 * @param [in] b Vector containing the defining polygon vertices
 * @param [out] p Vector containing the resulting curve points
 */
   void setBSplineSurface(size_t m, size_t n, size_t c, size_t d, size_t npu, size_t npw,
                          const Vect<real_t>& b, Vect<real_t>& p);

/**
 * \brief Define Bezier modelling
 * @param [in] n Number of defining polygon vertices
 * @param [in] nc Number of points to be calculated on the curve
 * @param [in] b Vector containing the defining polygon vertices
 * @param [out] p Vector containing the resulting curve points
 */
   void setBezier(size_t n, size_t nc, const Vect<real_t>& b, Vect<real_t>& p);

/**
 * \brief Define Bezier surface modelling
 * @param [in] m One less than the number of polygon vertices in w direction
 * @param [in] n One less than the number of polygon vertices in u direction
 * @param [in] npu Number of parametric lines in the u direction
 * @param [in] npw Number of parametric lines in the w direction
 * @param [in] b Vector containing the defining polygon vertices
 * @param [out] p Vector containing the resulting curve points
 */
   void setBezierSurface(size_t m, size_t n, size_t npu, size_t npw, const Vect<real_t>& b,
                         Vect<real_t>& p);

/**
 * \brief Define for Nurbs modelling
 * @param [in] n Number of defining polygon vertices
 * @param [in] c Order of the B-spline basis function
 * @param [in] np Number of points to be calculated on the curve
 * @param [in] b Vector containing the defining polygon vertices
 * @param [in] h Vector containing the homogeneous weighting factors
 * @param [out] p Vector containing the resulting curve points
 */
   void setNurbs(size_t n, size_t c, size_t np, const Vect<real_t>& b, const Vect<real_t>& h,
                 Vect<real_t>& p);

/**
 * \brief Define Nurbs surface modelling
 * @param [in] m Number of polygon vertices in w direction
 * @param [in] n Number of polygon vertices in u direction
 * @param [in] c Order of the B-spline basis function in w direction
 * @param [in] d Order of the B-spline basis function in u direction
 * @param [in] npu Number of parametric lines in the u direction
 * @param [in] npw Number of parametric lines in the w direction
 * @param [in] b Vector containing the defining polygon vertices
 * @param [out] p Vector containing the resulting curve points
 */
   void setNurbsSurface(size_t m, size_t n, size_t c, size_t d, size_t npu, size_t npw,
                        const Vect<real_t>& b, Vect<real_t>& p);

/**
 * \brief Run approximation process.
 * \details Results are stored in designated vector(s)
 * 
 * @return 0 if no problem occurred given in set...() functions. 
 */
   int run();

 private:

   Approx _ap;
   bool _alloc;
   size_t _n1, _n2, _c1, _c2, _np1, _np2;
   int _ls_opt, _degree;
   const vector<Fct *> *_fct;
   Fct *_ffct;
   Vect<real_t> *_p, _qq, _b;
   const Vect<real_t> *_h, *_q, *_x, *_y;
   real_t *_a0, *_a1;
   DMatrix<real_t> _A, *_B;

   real_t factrl(int n);
   void knot(size_t n, size_t c, vector<int>& x);
   real_t sumrbas(const vector<real_t>& N, const vector<real_t>& M);
   real_t BernsteinBasis(int n, int i, real_t t);
   void BSplineBasis(size_t n, size_t c, real_t t, const vector<int>& x, vector<real_t>& N);
   void RationalBasis(size_t c, real_t t, size_t n, const vector<int>& x, vector<real_t>& N);
   void SRationalBasis(size_t c, real_t t, size_t n, const vector<int>& x, vector<real_t>& N);
   int runLagrange();
   int runLeastSquare();
   int runBSpline();
   int runBSplineSurface();
   int runBezier();
   int runBezierSurface();
   int runNurbs();
   int runNurbsSurface();
   int rbsurf(size_t ibnum, vector<real_t>& bold, vector<real_t>& ni, vector<real_t>& mj,
              vector<real_t>& rsumij, vector<real_t>& savrsumij);

};

} /* namespace OFELI */
#endif