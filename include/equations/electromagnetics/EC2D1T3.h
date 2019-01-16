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

                         Definition of class EC2D1T3
     for Eddy Current Problems in Two-Dimensions with a scalar magnetic field
                           using the 3-Node triangle

  ==============================================================================*/


#ifndef __EC2D1T3_H
#define __EC2D1T3_H


#include "equations/electromagnetics/Equa_Electromagnetics.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file EC2D1T3.h
 *  \brief Definition file for class EC2D1T3.
 */

/*! \class EC2D1T3
 *  \ingroup Electromagnetics
 *  \brief Eddy current problems in 2-D domains using solenoidal approximation.
 *
 *  Builds finite element arrays for time harmonic eddy current problems in 2-D domains
 *  with solenoidal configurations (Magnetic field has only one nonzero component).
 *  Magnetic field is constant in the vacuum, and then zero in the outer vacuum.\n
 *  Uses 3-Node triangles.
 *
 *  The unknown is the time-harmonic magnetic induction (complex valued).
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class EC2D1T3 : public Equa_Electromagnetics<complex_t,3,3,2,2>
{

 public :

/// \brief Default constructor
    EC2D1T3();

/// \brief Constructor using element data
/// @param [in] el Pointer to Element instance
    EC2D1T3(const Element* el);

/// \brief Constructor using side data
/// @param [in] side Pointer to Side instance
    EC2D1T3(const Side* side);

/** \brief Constructor using element and previous time data
 *  @param [in] el Pointer to Element instance
 *  @param [in] u Solution at previous iteration
 *  @param [in] time Time value [Default: <tt>0</tt>]
 */
    EC2D1T3(const Element*         el,
            const Vect<complex_t>& u,
            const real_t&          time=0.);

/** \brief Constructor using side and previous time data
 *  @param [in] sd Pointer to Side instance
 *  @param [in] u Solution at previous iteration
 *  @param [in] time Time value [Default: <tt>0</tt>]
 */
    EC2D1T3(const Side*            sd,
            const Vect<complex_t>& u,
            const real_t&          time=0.);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    EC2D1T3(Mesh&            ms,
            Vect<complex_t>& u,
            real_t           omega,
            real_t           volt);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/// \brief Destructor
    ~EC2D1T3();

/** \brief Add magnetic contribution to matrix
 *  @param [in] omega Angular frequency
 *  @param [in] coef Coefficient to multiply by [Default: <tt>1</tt>]
 */
    void Magnetic(real_t omega,
                  real_t coef=1.);

/// \brief Add electric contribution to matrix
/// @param [in] coef Coefficient to multiply by [Default: <tt>1</tt>]
    void Electric(real_t coef=1.);

/// \brief Compute Joule density in element
    real_t Joule();

/// \brief Add element integral contribution
    complex_t IntegMF();

/** \brief Compute integral of normal derivative on edge
 *  @param [in] h Vect instance containing magnetic field at element nodes
 *  @param [in] opt Vector <tt>h</tt> is local (<tt>LOCAL_ARRAY</tt>) with size <tt>3</tt> 
 *  or global (<tt>GLOBAL_ARRAY</tt>) with size = Number of nodes
 *  [Default: <tt>GLOBAL_ARRAY</tt>].
 *  @note This member function is to be called within each element, it detects
 *  boundary sides as the ones with nonzero code
 */
    complex_t IntegND(const Vect<complex_t>& h,
                            int              opt=GLOBAL_ARRAY);

/// \brief Add contribution to vacuum area calculation
    real_t VacuumArea();

 protected:

    void set(const Element* el);
    void set(const Side* el);

 private:
    real_t _omega, _volt;

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
