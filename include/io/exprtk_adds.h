/*==============================================================================

                                     I  R  I  S

               An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2018 - 2020 Rachid Touzani

    This file is part of IRIS.

    IRIS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IRIS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================

             A utility function to add some constants to expressions

  ==============================================================================*/

#ifndef __EXPRTK_CONST_H
#define __EXPRTK_CONST_H

#include "exprtk.hpp"
#include "OFELI_Config.h"

namespace OFELI {

void add_constants(exprtk::symbol_table<real_t>& symbol_table);

}

#endif
