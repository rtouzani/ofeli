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

                           Definition of for class 'PostProcess'

  ==============================================================================*/

#ifndef __RECONSTRUCTION_H
#define __RECONSTRUCTION_H

#include <stdlib.h>
#include <math.h>

#include <iostream>
using std::ostream;
using std::endl;
using std::cerr;

#include <iomanip>
using std::setw;

#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"

namespace OFELI {

/*! \file Reconstruction.h
 *  \brief Definition file for class Reconstruction.
 */

/*! \class Reconstruction
 * \ingroup Solver
 * \brief To perform various reconstruction operations.
 *
 * \details This class enables various reconstruction operations like smoothing, projections, ...
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

template<class T_> class Vect;

class Reconstruction
{

 public:

//----------------------------   BASIC OPERATIONS   ----------------------------

/// \brief Default constructor
    Reconstruction() { }

/// \brief Constructor using a refrence to a Mesh instance
    Reconstruction(const Mesh& ms) { _theMesh = &ms; }

/// \brief Destructor
    ~Reconstruction() { }

/// Provide Mesh instance
    void setMesh(Mesh& ms) { _theMesh = &ms; }

/** \brief Smooth an elementwise field to obtain a nodewise field by L<sup>2</sup> projection
 *  @param [in] u Vect instance that contains field to smooth
 *  @param [out] v Vect instance that contains on output smoothed field
 */
    void P0toP1(const Vect<real_t>& u,
                Vect<real_t>&       v);

/** \brief Smooth an Discontinuous P1 field to obtain a nodewise (Continuous P<sub>1</sub>)
 *  field by L<sup>2</sup> projection
 *  @param [in] u Vect instance that contains field to smooth
 *  @param [out] v Vect instance that contains on output smoothed field
 *  @warning This function is valid for P<sub>1</sub> triangles (2-D) only.
 */
    void DP1toP1(const Vect<real_t>& u,
                 Vect<real_t>&       v);

    friend class Mesh;

 private:

    const Mesh   *_theMesh;
    Vect<real_t> _M;
};

} /* namespace OFELI */

#endif
