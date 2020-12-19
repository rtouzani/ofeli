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

            Class to compute locations and weights for Gauss quadrature
                                on the line [-1,+1]

  ==============================================================================*/

#ifndef __GAUSS_H
#define __GAUSS_H

#include "OFELI_Config.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/Point.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Gauss.h
 *  \brief Definition file for struct Gauss.
 */

/*! \class Gauss
 *  \ingroup Util
 * \brief Calculate data for %Gauss integration.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

template<class T_, size_t N_> class LocalVect;

class Gauss 
{

 public:

/// \brief Default constructor
    Gauss() { _triang = false; }

/// \brief Constructor using number of Gauss points
/// @param [in] np Number of integration points
    Gauss(size_t np);

/// \brief Set number of integration points
    void setNbPoints(size_t np) { _np = np; }

/** \brief Choose integration on triangle (7-point formula)
 *  \details If this is not selected, Gauss integration formula on <tt>[-1,1]</tt> is calculated.
 *  @param [out] w Array of weights of integration points
 *  @param [out] x Array of coordinates of integration points
 */
    void setTriangle(LocalVect<real_t,7>&        w,
                     LocalVect<Point<real_t>,7>& x);

/// \brief Return coordinate of <tt>i</tt>-th Gauss-Legendre point.
    real_t x(size_t i) const { return _x[i-1]; }

/// \brief Return coordinates of points in the reference triangle
    const Point<real_t> &xt(size_t i) const { return _xx[i-1]; }
    
/// \brief Return weight of <tt>i</tt>-th Gauss-Legendre point.
    real_t w(size_t i) const { return _w[i-1]; }

 private:
    size_t        _np;
    real_t        _x[10], _w[10];
    Point<real_t> _xx[10];
    bool _triang;
    void legendre(real_t y, real_t& p, real_t& dp);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
