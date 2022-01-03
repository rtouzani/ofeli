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

        Definition of Material properties associated to phase changes

  ==============================================================================*/


#include "OFELI_Config.h"
#include "equations/therm/PhaseChange.h"

namespace OFELI {


void PhaseChange::setMaterial(Material &m, int code)
{
   _material = &m;
   _code = code;
   _material->setCode(_code);
   _name = _material->getName(_code);
}


int PhaseChange::E2T(double &H, double &T, double &gamma)
{
   int ret = -1;
   if (_name=="Aluminium") {
      gamma = _material->Density()*_material->SpecificHeat();
      T = H/gamma;
      ret = 0;
   }
   else
      EnthalpyToTemperature(H,T,gamma);
   return ret;
}

} /* namespace OFELI */
