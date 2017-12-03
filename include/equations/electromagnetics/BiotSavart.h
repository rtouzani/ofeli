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

                        Definition of for class 'BiotSavart'

  ==============================================================================*/

#ifndef __BIOT_SAVART_H
#define __BIOT_SAVART_H

#include <stdlib.h>
#include <math.h>

#include "OFELI_Config.h"
#include <iostream>
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/LocalVect.h"
#include "mesh/MeshExtract.h"
#include "linear_algebra/Point.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file BiotSavart.h
 *  \brief Definition file for class BiotSavart.
 */

/*! \class BiotSavart
 * \ingroup Electromagnetics
 * \brief Class to compute the magnetic induction from the current density
 * using the Biot-Savart formula.
 * \details Given a current density vector given at elements, a collection of sides
 * of edges (piecewise constant),
 * this class enables computing the magnetic induction vector (continuous and
 * piecewise linear) using the Ampere equation. This magnetic induction is
 * obtained by using the Biot-Savart formula which can be either a volume,
 * surface or line formula depending on the nature of the current density vector.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */
class BiotSavart
{

 public:

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor.
    BiotSavart();

/// \brief Constructor using mesh data
/// @param [in] ms Mesh instance
    BiotSavart(Mesh& ms);

/** \brief Constructor using mesh and vector of real current density
 *  \details The current density is assumed piecewise constant
 *  @param [in] ms Mesh instance
 *  @param [in] J Sidewise vector of current density (<tt>J</tt> is a real valued vector),
 *  in the case of a surface supported current
 *  @param [in] B Nodewise vector that contains, once the member function run is used,
 *  the magnetic induction
 *  @param [in] code Only sides with given \a code support current [Default: <tt>0</tt>]
 */
    BiotSavart(      Mesh&         ms,
               const Vect<real_t>& J,
                     Vect<real_t>& B,
                     int           code=0);

/** \brief Constructor using mesh and vector of complex current density
 *  \details The current density is assumed piecewise constant
 *  @param [in] ms Mesh instance
 *  @param [in] J Sidewise vector of current density (<tt>J</tt> is a complex valued vector),
 *  in the case of a surface supported current
 *  @param [in] B Nodewise vector that contains, once the member function run is used,
 *  the magnetic induction
 *  @param [in] code Only sides with given <tt>code</tt> support current [Default: <tt>0</tt>]
 */
    BiotSavart(      Mesh&            ms,
               const Vect<complex_t>& J,
                     Vect<complex_t>& B,
                     int              code=0);

/// \brief Destructor
    ~BiotSavart();

//-------------------------------   MODIFIERS  ---------------------------------

/** \brief Set (real) current density given at elements
 *  \details The current density is assumed piecewise constant and real valued.
 *  This function can be used in the case of the volume Biot-Savart formula.
 *  @param [in] J Current density vector (Vect instance) and real entries
 */
    void setCurrentDensity(const Vect<real_t>& J);

/** \brief Set (real) current density given at elements
 *  \details The current density is assumed piecewise constant and complex valued.
 *  This function can be used in the case of the volume Biot-Savart formula.
 *  @param [in] J Current density vector (Vect instance) of complex entries
 */
    void setCurrentDensity(const Vect<complex_t>& J);

/** \brief Transmit (real) magnetic induction vector given at nodes
 *  @param [out] B Magnetic induction vector (Vect instance) and real entries
 */
    void setMagneticInduction(Vect<real_t>& B);
   
/** \brief Transmit (complex) magnetic induction vector given at nodes
 *  @param [out] B Magnetic induction vector (Vect instance) and complex entries
 */    
    void setMagneticInduction(Vect<complex_t>& B);

/// Choose code of faces or edges at which current density is given
    void selectCode(int code) { _code = code; }

/// \brief Set the magnetic permeability coefficient
/// @param [in] mu Magnetic permeability
    void setPermeability(real_t mu);

/** \brief Choose to compute the magnetic induction at boundary nodes only
 *  \details By default the magnetic induction is computed (using the function run)
 *  at all mesh nodes
 *  @note This function has no effect for surface of line Biot-Savart formula
 */
    void setBoundary() { _bound = true; }

/** \brief Compute the real magnetic induction at a given point using the volume Biot-Savart
 *  formula
 *  \details This function computes a real valued magnetic induction for a real valued
 *  current density field
 *  @param [in] x Coordinates of point at which the magnetic induction is computed
 *  @return Value of the magnetic induction at <tt>x</tt>
 */
    Point<real_t> getB3(Point<real_t> x);
   
/** \brief Compute the real magnetic induction at a given point using the surface Biot-Savart formula
 *  \details This function computes a real valued magnetic induction for a real valued
 *  current density field
 *  @param [in] x Coordinates of point at which the magnetic induction is computed
 *  @return Value of the magnetic induction at <tt>x</tt>
 */
    Point<real_t> getB2(Point<real_t> x);
   
/** \brief Compute the real magnetic induction at a given point using the line Biot-Savart formula
 *  \details This function computes a real valued magnetic induction for a real valued
 *  current density field
 *  @param [in] x Coordinates of point at which the magnetic induction is computed
 *  @return Value of the magnetic induction at <tt>x</tt>
 */
    Point<real_t> getB1(Point<real_t> x);
   
/** \brief Compute the complex magnetic induction at a given point using the volume Biot-Savart formula
 *  \details This function computes a complex valued magnetic induction for a complex valued
 *  current density field
 *  @param [in] x Coordinates of point at which the magnetic induction is computed
 *  @return Value of the magnetic induction at <tt>x</tt>
 */
    Point<complex_t> getBC3(Point<real_t> x);

/** \brief Compute the complex magnetic induction at a given point using the surface Biot-Savart formula
 *  \details This function computes a complex valued magnetic induction for a complex valued
 *  current density field
 *  @param [in] x Coordinates of point at which the magnetic induction is computed
 *  @return Value of the magnetic induction at <tt>x</tt>
 */
    Point<complex_t> getBC2(Point<real_t> x);

/** \brief Compute the complex magnetic induction at a given point using the line Biot-Savart formula
 *  \details This function computes a complex valued magnetic induction for a complex valued
 *  current density field
 *  @param [in] x Coordinates of point at which the magnetic induction is computed
 *  @return Value of the magnetic induction at <tt>x</tt>
 */
    Point<complex_t> getBC1(Point<real_t> x);

/** \brief Run the calculation by the Biot-Savart formula
 *  \details This function computes the magnetic induction, which is stored in the
 *  vector <tt>B</tt> given in the constructor
 */
    int run();

//-----------------------------   INSPECTORS  ----------------------------------


 private:

   bool _C;
   int _type, _code, _bound;
   real_t _mu;
   Mesh *_theMesh;
   SideList *_theSideList;
   EdgeList *_theEdgeList;
   const Vect<real_t> *_J;
   const Vect<complex_t> *_JC;
   Vect<real_t> *_B, _meas;
   Vect<complex_t> *_BC;
   Vect<Point<real_t> > _center;
   void setVolume();
   void setSurface();
   void setLine();
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
