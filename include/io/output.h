/*==============================================================================

                                    O  F  E  L  I

                          Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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
                          Some utility functions for output
  ==============================================================================*/

#ifndef __OUTPUT_H
#define __OUTPUT_H

#include "OFELI_Config.h"
#include <iostream>
using std::ostream;

#include <iomanip>
#include <complex>
#include <string>
#include <list>

#include <vector>
using std::vector;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file output.h
 *  \brief File that contains some output utility functions.
 */


/** \fn ostream & operator<<(ostream& s, const std::complex<double> &x)
 *  \brief Output a complex number.
 * \ingroup Util
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
inline ostream& operator<<(ostream&         s,
                           const complex_t& x)
{
   if (x.imag()<0)
      s << x.real() << " - " << -x.imag() << " I";
   else
      s << x.real() << " + " << x.imag() << " I";
   return s;
}


/** \fn ostream & operator<<(ostream& s, const std::string &c)
 *  \brief Output a string.
 *  \ingroup Util
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
inline ostream& operator<<(ostream&           s,
                           const std::string& c)
{
    for (size_t i=0; i<c.length(); i++)
       s << c[i];
    return s;
}


/** \fn ostream & operator<<(ostream& s, const vector<T_> &v)
 *  \brief Output a vector instance.
 *  \ingroup Util
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template<class T_>
ostream&operator<<(ostream&          s,
                   const vector<T_>& v)
{
   for (size_t i=0; i<v.size(); i++)
      s << i << ": " << v[i] << std::endl;
   return s;
}


/** \fn ostream & operator<<(ostream& s, const std::list<T_> &l)
 *  \brief Output a vector instance.
 *  \ingroup Util
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template<class T_>
ostream&operator<<(ostream&             s,
                   const std::list<T_>& l)
{
   typename std::list<T_>::const_iterator it = l.begin();
   while (it != l.end()) {
      s << *it << std::endl;
      it++;
   }
   return s;
}


/** \fn ostream & operator<<(ostream& s, const std::pair<T_,T_> &a)
 *  \brief Output a pair instance.
 *  \ingroup Util
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
template<class T_>
inline ostream& operator<<(ostream&                s,
                           const std::pair<T_,T_>& a)
{
   s << "(" << a.first << "," << a.second << ") ";
   return s;
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
