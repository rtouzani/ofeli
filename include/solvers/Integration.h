/*==============================================================================

                                    O  F  E  L  I

                          Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                    Class Integration for Numerical Integration

  ==============================================================================*/


#ifndef __INTEGRATION_H
#define __INTEGRATION_H

#include <functional>
using std::function;

#include "OFELI_Config.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Integration.h
 *  \brief Definition file for numerical integration class.
 */

/*! \enum IntegrationScheme
 * Choose numerical integration scheme
 */
enum IntegrationScheme {
   LEFT_RECTANGLE  = 0,   /*!< Left rectangle integration formula */
   RIGHT_RECTANGLE = 1,   /*!< Right rectangle integration formula   */
   MID_RECTANGLE   = 2,   /*!< Midpoint (central) rectangle formula       */
   TRAPEZOIDAL     = 3,   /*!< Trapezoidal rule     */
   SIMPSON         = 4,   /*!< Simpson formula     */
   GAUSS_LEGENDRE  = 5,   /*!< Gauss-Legendre quadrature formulae            */
};

/*! \class Integration
 *  \ingroup Solver
 *  \brief Class for numerical integration methods
 *  \details Class NumInt defines and stores numerical integration data
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Integration
{

 public:

/// \brief Default constructor.
    Integration();

/** \brief Constructor 
 *  @param [in] low Lower value of integration interval
 *  @param [in] high Upper value of integration interval
 *  @param [in] f Function to integrate
 *  @param [in] s Integration scheme. To choose among enumerated values:
 *  <ul>
 *    <li><tt>LEFT_RECTANGLE</tt>: 
 *    <li><tt>RIGHT_RECTANGLE</tt>: 
 *    <li><tt>MID_RECTANGLE</tt>: 
 *    <li><tt>TRAPEZOIDAL</tt>: 
 *    <li><tt>SIMPSON</tt>: 
 *    <li><tt>GAUSS_LEGENDRE</tt>: 
 *  </ul>
 *  @param [in] error
 */
    Integration(real_t                          low,
                real_t                          high,
                function<real_t(real_t)> const& f,
                IntegrationScheme               s,
                real_t                          error);

/// \brief Destructor
    ~Integration() { }

/** \brief Define function to integrate numerically
 *  @param [in] f Function to integrate
 */
    void setFunction(function<real_t(real_t)> const& f);

/** \brief Set time inegration scheme
 *  @param [in] s Scheme to choose among enumerated values:
 *  <ul>
 *    <li><tt>LEFT_RECTANGLE</tt>: 
 *    <li><tt>RIGHT_RECTANGLE</tt>: 
 *    <li><tt>MID_RECTANGLE</tt>: 
 *    <li><tt>TRAPEZOIDAL</tt>: 
 *    <li><tt>SIMPSON</tt>: 
 *    <li><tt>GAUSS_LEGENDRE</tt>: 
 *  </ul>
 */
    void setScheme(IntegrationScheme s) { _scheme = s; }

/** \brief Define integration domain as a quadrilateral
 *  \details
 *  @param [in] x1 x-coordinate of first vertex of triangle
 *  @param [in] y1 y-coordinate of first vertex of triangle
 *  @param [in] x2 x-coordinate of second vertex of triangle
 *  @param [in] y2 y-coordinate of second vertex of triangle
 *  @param [in] x3 x-coordinate of third vertex of triangle
 *  @param [in] y3 y-coordinate of third vertex of triangle
 */
    void setTriangle(real_t x1,
                     real_t y1,
                     real_t x2,
                     real_t y2,
                     real_t x3,
                     real_t y3);

/** \brief Define integration domain as a quadrilateral
 *  \details
 *  @param [in] x1 x-coordinate of first vertex of quadrilateral
 *  @param [in] y1 y-coordinate of first vertex of quadrilateral
 *  @param [in] x2 x-coordinate of second vertex of quadrilateral
 *  @param [in] y2 y-coordinate of second vertex of quadrilateral
 *  @param [in] x3 x-coordinate of third vertex of quadrilateral
 *  @param [in] y3 y-coordinate of third vertex of quadrilateral
 *  @param [in] x4 x-coordinate of fourth vertex of quadrilateral
 *  @param [in] y4 y-coordinate of fourth vertex of quadrilateral
 */
    void setQuadrilateral(real_t x1,
                          real_t y1,
                          real_t x2,
                          real_t y2,
                          real_t x3,
                          real_t y3,
                          real_t x4,
                          real_t y4);

/** \brief Run numerical integration
 *  @return Computed approximate value of integral
 */
    real_t run();

 private:

    IntegrationScheme _scheme;
    real_t            _low, _high, _err;
    Point<real_t>     _x1, _x2, _x3, _x4;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
