/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                         Definition of class EC2DB2T3
     for Eddy Current Problems in Two-Dimensions with a vector magnetic field
                           using the 3-Node triangle
                      The Problem is solved in a bounded domain
                      by using a transparent boundary condition

  ==============================================================================*/


#ifndef __EC2DB2T3_H
#define __EC2DB2T3_H


#include "equations/electromagnetics/Equa_Electromagnetics.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/SkMatrix.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


class EC2DB2T3 : virtual public Equa_Electromagnetics<3,6,2,4>
{

 public:

/// \brief Default Constructor
    EC2DB2T3() { }

/// \brief Constructor using element data
    EC2DB2T3(Element* el);

/// \brief Constructor using one side data
    EC2DB2T3(Side* sd);

/// \brief Destructor
    ~EC2DB2T3() { }

/// \brief Compute Contribution to Right-Hand Side
/// @param [in] coef Coefficient to multiply by the RHS before adding [Default: <tt>1</tt>]
    void RHS(real_t coef=1.);

/// \brief Compute Finite Element Diagonal Block
    void EMatr();

/** \brief Compute constant to multiply by potential
 *  @param [in] u Vector containing potential at nodes
 *  @param [in] I Prescribed total current
 */
    complex_t Constant(const LocalVect<real_t,6>& u,
                       complex_t                  I);

/// \brief Compute magnetic pressure in element
/// @param [in] u Vector containing potential at nodes
    real_t MagneticPressure(const LocalVect<real_t,6>& u);

 private:

   void set(const Side* sd);
   void set(const Element* el);
   LocalVect<Point<real_t>,3> _N;
   size_t   _ns;
   real_t   _area, _det, _length;
   Point<real_t> _center;
   LocalVect<Point<real_t>,3> _x;
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
