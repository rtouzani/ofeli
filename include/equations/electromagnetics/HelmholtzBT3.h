/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                         Definition of class HelmholtzBT3
               for Helmholtz Equation in a bounded domain using
                         3-node triangular finite element

  ==============================================================================*/


#ifndef __HelmholtzBT3_H
#define __HelmholtzBT3_H


#include "equations/electromagnetics/Equa_Electromagnetics.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "io/UserData.h"

using std::complex;

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file HelmholtzBT3.h
 *  \brief Definition file for class HelmholtzBT3.
 *
 */

/*! \class HelmholtzBT3
 *  \ingroup Electromagnetics
 *  \brief Builds finite element arrays for Helmholtz equations in a bounded media
 *  using 3-Node triangles.
 *
 *  \details Problem being formulated in time harmonics, the solution is complex valued.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class HelmholtzBT3 : virtual public Equa_Electromagnetics<complex_t,3,3,2,2>
{

 public:

/// \brief Default Constructor
    HelmholtzBT3()
    {
       _tr = NULL; _ln = NULL;
    }

/// \brief Constructor using element data
    HelmholtzBT3(Element* el);

/// \brief Constructor using side data
    HelmholtzBT3(Side* sd);

/// \brief Destructor
    ~HelmholtzBT3();

/// \brief Add element Left-Hand Side
    void LHS(real_t wave_nb);

/// \brief Add element Right-Hand Side using a UserData instance
    void BoundaryRHS(UserData<complex_t>& ud);

 private:

   Triang3      *_tr;
   Line2        *_ln;
   real_t        _lx[3], _ly[3];
   Point<real_t> _x[3];
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
