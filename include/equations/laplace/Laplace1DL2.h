/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                      Definition of class Laplace1DL2
              for 1-D Laplace equation using 2-node line element (P1)

  ==============================================================================*/


#ifndef __LAPLACE_1DL2_H
#define __LAPLACE_1DL2_H

#include <float.h>
#include <stdlib.h>
#include <math.h>

#include "OFELI_Config.h"
#include "equations/laplace/Equa_Laplace.h"
#include "mesh/Element.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Laplace1DL2.h
 *  \brief Definition file for class Laplace1DL2.
 */

/*! \class Laplace1DL2
 *  \ingroup Laplace
 *  \brief To build element equation for a 1-D elliptic equation
 *  using the 2-Node line element (<tt>P<sub>1</sub></tt>).
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Laplace1DL2 : public Equa_Laplace<2,2,1,1> {

 public :

    using Equa::setInput;

/// \brief Default constructor
    Laplace1DL2();

/** Constructor using mesh instance and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Vect instance that contains, after execution of \b run() the solution
 */
    Laplace1DL2(Mesh&         ms,
                Vect<real_t>& u);

/// Constructor using mesh instance
/// @param [in] ms Mesh instance
    Laplace1DL2(Mesh& ms);

/// \brief Destructor
    ~Laplace1DL2() { }

/// \brief Add finite element matrix to left hand side
    void LHS();

/** \brief Build global stiffness and mass matrices for the eigen system
 *  @param [in] opt Flag to choose a lumped mass matrix (0) or consistent (1) [Default: <tt>0</tt>]
 */
    void buildEigen(int opt=0);

/// \brief Add Right-Hand %Side Contribution
/// @param [in] f Vector containing the source given function at mesh nodes
    void BodyRHS(const Vect<real_t>& f);

/** \brief Add Neumann contribution to Right-Hand %Side
 *  @param [in] f Vector with size the total number of nodes. The first entry stands for
 *  the force at the first node (Neumann condition) and the last entry is the force at
 *  the last node (Neumann condition)
 */
    void BoundaryRHS(const Vect<real_t>& f);

/** \brief Set Dirichlet boundary data
 *  @param [in] f Value to assign
 *  @param [in] lr Option to choose location of the value (<tt>-1</tt>: Left end,
 *  <tt>1</tt>: Right end)
 */
    void setBoundaryCondition(real_t f,
                              int    lr);

/** \brief Set Traction data
 *  @param [in] f Value of traction (Neumann boundary condition)
 *  @param [in] lr Option to choose location of the traction (<tt>-1</tt>: Left end, 
 *  <tt>1</tt>: Right end)
 */
    void setTraction(real_t f,
                     int    lr);

    
 private:
    real_t _lsf, _rsf, _lbc, _rbc;
    bool   _is_lbc, _is_rbc;
    void set(const Element* el);
    void set(const Side* sd) { }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
