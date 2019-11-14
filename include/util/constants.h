/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                                  Some Constants

  ==============================================================================*/

#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <float.h>

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file constants.h
 *  \ingroup Util
 *  \brief File that contains some widely used constants.
 */

/*! \def OFELI_E
 *  \ingroup Util
 * Value of \a e or \a exp (with 28 digits)
 */
#ifndef OFELI_E
#define OFELI_E         2.71828182845904523536028747135
#endif

/*! \def OFELI_PI
 *  \ingroup Util
 * Value of \a Pi (with 28 digits)
 */
#ifndef OFELI_PI
#define OFELI_PI        3.14159265358979323846264338328
#endif

/*! \def OFELI_THIRD
 *  \ingroup Util
 * Value of \a 1/3 (with 28 digits)
 */
#ifndef OFELI_THIRD
#define OFELI_THIRD     0.33333333333333333333333333333
#endif

/*! \def OFELI_SIXTH
 *  \ingroup Util
 * Value of \a 1/6 (with 28 digits)
 */
#ifndef OFELI_SIXTH
#define OFELI_SIXTH     0.16666666666666666666666666667
#endif

/*! \def OFELI_TWELVETH
 *  \ingroup Util
 * Value of \a 1/12 (with 28 digits)
 */
#ifndef OFELI_TWELVETH
#define OFELI_TWELVETH  0.08333333333333333333333333333
#endif

/*! \def OFELI_SQRT2
 *  \ingroup Util
 * Value of \a sqrt(2) (with 28 digits)
 */
#ifndef OFELI_SQRT2
#define OFELI_SQRT2     1.41421356237309504880168872421
#endif

/*! \def OFELI_SQRT3
 *  \ingroup Util
 * Value of \a sqrt(3) (with 28 digits)
 */
#ifndef OFELI_SQRT3
#define OFELI_SQRT3     1.73205080756887729352744634151
#endif

/*! \def OFELI_ONEOVERPI
 *  \ingroup Util
 * Value of \a 1/Pi (with 28 digits)
 */
#ifndef OFELI_ONEOVERPI
#define OFELI_ONEOVERPI 0.31830988618379067153776752675
#endif

/*! \def OFELI_GAUSS2
 *  \ingroup Util
 * Value of \a 1/sqrt(3) (with 32 digits)
 */
#ifndef OFELI_GAUSS2
#define OFELI_GAUSS2    0.57735026918962576450914878050196
#endif

/*! \def OFELI_EPSMCH
 *  \ingroup Util
 * Value of Machine Epsilon
 */
#ifndef OFELI_EPSMCH
#define OFELI_EPSMCH    DBL_EPSILON
#endif

/*! \def OFELI_TOLERANCE
 *  \ingroup Util
 * Default tolerance for an iterative process
 * = OFELI_EPSMCH * 10000
 */
#define OFELI_TOLERANCE OFELI_EPSMCH*10000

/*! \def VLG
 *  \ingroup Util
 * Very large number: A real number for penalty
 */
#ifndef VLG
#define VLG  1.e10
#endif

/*! \def OFELI_IMAG
 *  \ingroup Util
 * = Unit imaginary number (\c i)
 */
#define OFELI_IMAG std::complex<double>(0.,1.);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
