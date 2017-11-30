/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

      Definition of class 'OFELIException' for handling exceptions in OFELI

  ==============================================================================*/

#ifndef __OFELI_EXCEPTION_H
#define __OFELI_EXCEPTION_H

/*!
 * \file OFELIException.h
 * \brief 
 */

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


/*!
 * \class OFELIException
 * \ingroup 
 * \brief To handle exceptions in OFELI
 * \details This class enables using exceptions in programs using OFELI
 */
class OFELIException : public std::runtime_error
{

 public:

/// \brief This form will be used most often in a throw.
    OFELIException(const string& description) : runtime_error(description)
    { };

/// \brief Throw with no error message.
    OFELIException() : std::runtime_error("Exception thrown in OFELI:  ") { }; 
};

} /* namespace OFELI */

#endif
