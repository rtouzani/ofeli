/*==============================================================================

                                    O  F  E  L  I

                             Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                    Definition of class 'GeoModel' for various
                          methods of geometry modelling

  ==============================================================================*/

#ifndef __GEO_MODEL_H
#define __GEO_MODEL_H

#include "OFELI_Config.h"
#include "linear_algebra/Vect.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file GeoModel.h
 *  \brief Definition file for class GeoModel.
 */

/*! \class GeoModel
 * \ingroup Solver
 * \brief To set geometry modelling.
 *
 * This class enables using approximation methods to mathematically define
 * a geometry.
 * 
 * The algorithms used in this class are largely inspired from the book:
 * An Introduction to NURBS, by David F. Rogers. Copyright (C) 2000 David F. Rogers,
 */

#define NURBS_THRESHOLD 5.0e-6

class GeoModel
{

 public:

/// \brief Default constructor
/// \details The function setData can then be used
    GeoModel();

/** \brief Constructor with given polygon points and solution vector
 *  \details This function is to be used if the default constructed was used
 *  @param [in] b Vector containing the defining polygon vertices
 *  @param [in,out] p Vector containing the resulting curve points
 */
   GeoModel(const Vect<real_t>& b, Vect<real_t>& p);

/** \brief Constructor with given data for nurbs
 *  \details This function is to be used if the default constructed was used
 *  @param [in] b Vector containing the defining polygon vertices
 *  @param [in] h Vector containing the homogeneous weighting factors
 *  @param [in,out] p Vector containing the resulting curve points
 */
   GeoModel(const Vect<real_t>& b, const Vect<real_t>& h, Vect<real_t>& p);

/// \brief Destructor
   ~GeoModel() { }

/** \brief Set vector data
 *  @param [in] b Vector containing the defining polygon vertices
 *  @param [in,out] p Vector containing the resulting curve points
 */
   void setData(const Vect<real_t>& b, Vect<real_t>& p);

/** \brief Set vector data for Nurbs
 *  @param [in] b Vector containing the defining polygon vertices
 *  @param [in] h Vector containing the homogeneous weighting factors
 *  @param [in,out] p Vector containing the resulting curve points
 */
   void setData(const Vect<real_t>& b, const Vect<real_t>& h, Vect<real_t>& p);

/** \brief Set parameters for BSpline modelling
 * @param [in] n Number of defining polygon vertices
 * @param [in] c Order of the B-spline basis function
 * @param [in] np Number of points to be calculated on the curve
 */
   void setBSplinePar(size_t n, size_t c, size_t np);

/** \brief Set parameters for BSplineS modelling
 * @param [in] m One less than the number of polygon vertices in w direction
 * @param [in] n One less than the number of polygon vertices in u direction
 * @param [in] c Order of the B-spline basis function in w direction
 * @param [in] d Order of the B-spline basis function in u direction
 * @param [in] npu Number of parametric lines in the u direction
 * @param [in] npw Number of parametric lines in the w direction
 */
   void setBSplineSurfacePar(size_t m, size_t n, size_t c, size_t d, size_t npu, size_t npw);

/** \brief Set parameters for BSpline modelling
 * @param [in] n Number of defining polygon vertices
 * @param [in] nc Number of points to be calculated on the curve
 */
  void setBezierPar(size_t n, size_t nc);

/** \brief Set parameters for BSpline modelling
 *  @param [in] m One less than the number of polygon vertices in w direction
 *  @param [in] n One less than the number of polygon vertices in u direction
 *  @param [in] npu Number of parametric lines in the u direction
 *  @param [in] npw Number of parametric lines in the w direction
 */
   void setBezierSurfacePar(size_t m, size_t n, size_t npu, size_t npw);

/** \brief Set parameters for Nurbs modelling
 *  @param [in] n Number of defining polygon vertices
 *  @param [in] c Order of the B-spline basis function
 *  @param [in] np Number of points to be calculated on the curve
 */
   void setNurbsPar(size_t n, size_t c, size_t np);

/** \brief Run bspline modelling
 * @remark The resulting vector of curve points is the vector \c p given
 * by the constructor or by setData
 */
   void BSpline();

/** \brief Run surface bspline modelling
 * @remark The resulting vector of curve points is the vector \c p given
 * by the constructor or by setData
 */
   void BSplineSurface();

/** \brief Run Bezier modelling
 * @remark The resulting vector of curve points is the vector \c p given
 * by the constructor or by setData
 */
   void Bezier();

/** \brief Run Surface Bezier modelling
 * @remark The resulting vector of curve points is the vector \c p given
 * by the constructor or by setData
 */
   void BezierSurface();

/** \brief Run Nurbs modelling
 * @remark The resulting vector of curve points is the vector \c p given
 * by the constructor or by setData
 */
   void Nurbs();

 private:

   size_t _nu, _nw, _cu, _cw, _npu, _npw;
   Vect<real_t> *_cp;
   const Vect<real_t> *_h, *_pv;

   real_t factrl(int n);
   void knot(size_t n, size_t c, vector<int>& x);
   real_t BernsteinBasis(int n, int i, real_t t);
   void BSplineBasis(size_t n, size_t c, real_t t, const vector<int>& x, vector<real_t>& N);
   void RationalBasis(size_t c, real_t t, size_t n, const vector<int>& x, vector<real_t>& N);

};

} /* namespace OFELI */

#endif
