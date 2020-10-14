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

                     Some additional functions to use with exprtk 

  ==============================================================================*/

#include "io/exprtk_adds.h"
#include <math.h>

exprtk::parser<real_t> theParser;

namespace OFELI {

void add_constants(exprtk::symbol_table<real_t>& symbol_table)
{
   real_t pi=3.14159265358979323846264338328, e=2.71828182845904523536028747135;
   symbol_table.add_constant("pi",pi);
   symbol_table.add_constant("e",e);
   symbol_table.add_constants();
}

}
