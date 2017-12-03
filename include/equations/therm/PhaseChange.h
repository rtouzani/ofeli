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

         Definition of Material properties associated to phase change

  ==============================================================================*/


#ifndef __PHASE_CHANGE_H
#define __PHASE_CHANGE_H

#include "mesh/Material.h"

#include <string>
using std::string;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file PhaseChange.h
 *  \brief Definition file for class PhaseChange and its parent abstract class.
 */

/*! \class PhaseChange
 *  \ingroup Therm
 *  \brief This class enables defining phase change laws for a given material.
 *  \details These laws are predefined for a certain number of materials.
 *  The user can set himself a specific behavior for his own materials
 *  by defining a class that inherits from PhaseChange. The derived
 *  class must has at least the member function
 *  
 *  int EnthalpyToTemperature(real_t &H, real_t &T, real_t &gamma)
 *
 */

class PhaseChange
{

  public:

/// Destructor
    virtual ~PhaseChange() { }

/** \brief Calculate temperature from enthalpy.
 *  \details This member function is to be called in any equation class that needs phase change laws.
 *  @param [in] H Enthalpy value
 *  @param [out] T Calculated temperature value
 *  @param [out] gamma Maximal slope of the curve <tt>H -> T</tt>
 */
    int E2T(real_t& H,
            real_t& T,
            real_t& gamma);

/** \brief Virtual function to calculate temperature from enthalpy.
 *  \details This member function must be implemented in any derived class
 *  in order to define user's own material laws.
 *  @param [in] H Enthalpy value
 *  @param [out] T Calculated temperature value
 *  @param [out] gamma Maximal slope of the curve <tt>H -> T</tt>
 */
    virtual int EnthalpyToTemperature(real_t& H,
                                      real_t& T,
                                      real_t& gamma)
    {
       gamma = _material->Density()*_material->SpecificHeat();
       T = H/gamma;
       return 0;
    }

/// \brief Choose Material instance and material code
    void setMaterial(Material& m,
                     int       code);

/// \brief Return reference to Material instance
    Material &getMaterial() const { return *_material; }

  protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    Material *_material;
    int       _code;
    string    _name;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
};

template<class C_>
void setPhaseChangeLaw(C_&     c,
                       real_t& H,
                       real_t& T,
                       real_t& gamma);

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
