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

                         Definition of class HelmholtzBT3
               for Helmholtz Equation in a bounded domain using
                         3-node triangular finite element

  ==============================================================================*/


#ifndef __HelmholtzBT3_H
#define __HelmholtzBT3_H


#include "equations/electromagnetics/Equa_Electromagnetics.h"


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
 *  \details Problem being formulated in time harmonics, the solution is complex-valued
 *  but stored in 2-degree of freedom real-valued vector.
 *  Therefore, mesh must be defined with 2 degrees of freedom per node
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class HelmholtzBT3 : virtual public Equa_Electromagnetics<3,6,2,4>
{

 public:

/// \brief Default Constructor
    HelmholtzBT3();

/** \brief Constructor using mesh data
 *  @param [in] ms Mesh instance
 */
    HelmholtzBT3(Mesh& ms);

/** \brief Constructor using mesh and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance containing solution
 */
    HelmholtzBT3(Mesh&         ms,
                 Vect<real_t>& u);

/// \brief Destructor
    ~HelmholtzBT3();

/// \brief Builds system of equations
    void build();

/// \brief Add element Left-Hand Side
    void LHS();

/// \brief Add element Right-Hand Side
    void BodyRHS(Vect<real_t>& f);

/// \brief Add side Right-Hand Side
    void BoundaryRHS(Vect<real_t>& f);

 private:
   void set(const Element* el);
   void set(const Side* sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
