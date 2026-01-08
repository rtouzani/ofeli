/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

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

                       Definition of for class 'Funct'

  ==============================================================================*/

#ifndef __FUNCT_H
#define __FUNCT_H

#include <stdlib.h>
#include <math.h>
#include <string>
using std::string;

#include "OFELI_Config.h"
#include "io/exprtk.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
extern exprtk::parser<real_t> theParser;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file FunctTab.h
 *  \brief Definition file for class FunctTab.
 */
 

/*! \class FunctTab
 * \ingroup Util
 * \brief A simple class to define functions by a tabulation.
 * \details Functions must have 1, 2, 3 or at most 4 variables.
 * \warning Data in the file must be listed in the following order:
 * \verbatim
      for x=x_0,...,x_I
          for y=y_0,...,y_J
              for z=z_0,...,z_K
                  read v(x,y,z)
   \endverbatim
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */


class FunctTab
{

 public:

/// \brief Default constructor
    FunctTab() { }

/// \brief Constructor for a function of one variable
/// @param [in] v Name of the variable
    FunctTab(string v)
    {
       _vn = v;
    }

/// \brief Constructor for a function of two variables
/// @param [in] v1 Name of the first variable
/// @param [in] v2 Name of the second variable
    FunctTab(string v1,
             string v2)
    {
        _vn = v1 + string(",") + v2;
    }

/// \brief Constructor for a function of three variables
/// @param [in] v1 Name of the first variable
/// @param [in] v2 Name of the second variable
/// @param [in] v3 Name of the third variable
    FunctTab(string v1,
             string v2,
             string v3)
    {
       _vn = v1 + string(",") + v2 + string(",") + v3;
    }

/// \brief Constructor for a function of four variables
/// @param [in] v1 Name of the first variable
/// @param [in] v2 Name of the second variable
/// @param [in] v3 Name of the third variable
/// @param [in] v4 Name of the fourth variable
    FunctTab(string v1,
             string v2,
             string v3,
             string v4)
    {
       _vn = v1 + string(",") + v2 + string(",") + v3 + string(",") + v4;
    }

/// \brief Destructor
    ~FunctTab() { }

/// \brief Operator () to evaluate the function with one variable <tt>x</tt>
    real_t operator()(real_t x) const { return theParser.Eval(x); }

/// \brief Operator () to evaluate the function with two variables <tt>x</tt>, <tt>y</tt>
    real_t operator()(real_t x,
                      real_t y) const { return theParser.Eval(x,y); }

/// \brief Operator () to evaluate the function with three variables <tt>x</tt>, <tt>y</tt>, <tt>z</tt>
    real_t operator()(real_t x,
                      real_t y,
                      real_t z) const { return theParser.Eval(x,y,z); }

/// \brief Operator () to evaluate the function with four variables <tt>x</tt>, <tt>y</tt>,
/// <tt>z</tt>
    real_t operator()(real_t x,
                      real_t y,
                      real_t z,
                      real_t t) const { return theParser.Eval(x,y,z,t); }

/** \brief Operator =.
 *  \details Define the function by an algebraic expression
 *  @param [in] e Algebraic expression defining the function.
 */
    void operator=(string e) { theParser.Parse(e.c_str(),_vn.c_str()); }


 private:

   string _vn;
};

/*! @} End of Doxygen Groups */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
} /* namespace OFELI */

#endif
