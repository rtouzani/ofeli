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

                             Definition of class Beam3DL2
                       for Beam element with 6 d.o.f. per node
                           using Mindlin - Reissner Theory

  ==============================================================================*/


#ifndef __BEAM3DL2_H
#define __BEAM3DL2_H

#include "equations/solid/Equa_Solid.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Beam3DL2.h
 *  \brief Definition file for class Beam3DL2.
 */

/*! \class Beam3DL2
 *  \ingroup Solid
 *  \brief To build element equations for 3-D beam equations using 2-node lines
 *
 *  \details This class enables building finite element arrays for 3-D beam elements using
 *  6 degrees of freedom per node and 2-Node line elements.
 */

class Beam3DL2 : virtual public Equa_Solid<2,12,1,6> {

public :

   using Equation<2,12,1,6>::_x;

/// \brief Default Constructor
    Beam3DL2() {
       _bending = _axial = _torsion = _shear = true;
       _reduced_integration = false;
    }

/** \brief Constructor using mesh and constant beam properties
    @param [in] ms Mesh instance
    @param [in] A Section area of the beam
    @param [in] I1 first (x) momentum of inertia
    @param [in] I2 second (y) momentum of inertia
 */
    Beam3DL2(Mesh&  ms,
             real_t A,
             real_t I1,
             real_t I2);

/** \brief Constructor using a Mesh instance
 *  @param [in] ms Reference to Mesh instance
 */
    Beam3DL2(Mesh& ms);

/** \brief Constructor using a Mesh instance and solution vector
 *  @param [in] ms Reference to Mesh instance
 *  @param [in,out] u Solution vector
 */
    Beam3DL2(Mesh&         ms,
             Vect<real_t>& u);

/// \brief Destructor
    ~Beam3DL2() { }

/** \brief Set constant beam properties
    @param [in] A Section area of the beam
    @param [in] I1 first (x) momentum of inertia
    @param [in] I2 second (y) momentum of inertia
 */
    void set(real_t A,
             real_t I1,
             real_t I2);

/** \brief Set nonconstant beam properties
 *  @param [in] A Vector containing section areas of the beam (for each element)
 *  @param [in] I1 Vector containing first (x) momentum of inertia (for each element)
 *  @param [in] I2 Vector containing second (y) momentum of inertia (for each element)
 */
    void set(const Vect<real_t>& A,
             const Vect<real_t>& I1,
             const Vect<real_t>& I2);

/** \brief Get vector of displacements at nodes
 *  @param [out] d Vector containing three components for each node that are <tt>x</tt>,
 *  <tt>y</tt> and <tt>z</tt> displacements.
 */
    void getDisp(Vect<real_t>& d);

/// \brief Add element lumped Mass contribution to element matrix after multiplication by <tt>coef</tt>
    void LMass(real_t coef=1.);

/// \brief Add element consistent Mass contribution to RHS after multiplication by 
/// <tt>coef</tt> (not implemented)
    void Mass(real_t coef=1.) { cerr << "Error: Beam3DL2::Mass not implemented" << endl; }

/// \brief Add element stiffness to element matrix
    void Stiffness(real_t coef=1.);

/// \brief Add contributions for loads
    void Load(const Vect<real_t>& f);

// Return Stress in beam
//    real_t Stress();

/// \brief Set bending contribution to stiffness
    void setBending() { _bending = true; }

/// \brief Set axial contribution to stiffness
    void setAxial() { _axial = true; }

/// \brief Set shear contribution to stiffness
    void setShear() { _shear = true; }

/// \brief Set torsion contribution to stiffness
    void setTorsion() { _torsion = true; }

/// \brief Set no bending contribution
    void setNoBending() { _bending = false; }

/// \brief Set no axial contribution
    void setNoAxial() { _axial = false; }

/// \brief Set no shear contribution
    void setNoShear() { _shear = false; }

/// \brief Set no torsion contribution
    void setNoTorsion() { _torsion = false; }

/// \brief Set reduced integration
    void setReducedIntegration() { _reduced_integration = true; }

/** \brief Return axial force in element
 *  @param [out] f Vector containing axial force in each element. This vector is
 *                 resized in the function
 */
    void AxialForce(Vect<real_t>& f);

/** \brief Return shear force in element
 *  @param [out] sh Vector containing shear forces (2 components) in each element. This vector is
 *               resized in the function
 */
    void ShearForce(Vect<real_t>& sh);

/** \brief Return bending moment in element
 *  @param [out] m Vector containing bending moments (2 components) in each element. This vector is
 *                 resized in the function
 */
    void BendingMoment(Vect<real_t>& m);

/** \brief Return twisting moments
 *  @param [out] m Vector containing twisting moment in each element. This vector is
 *                 resized in the function
 */
    void TwistingMoment(Vect<real_t>& m);

/// \brief Build the linear system of equations
    void build();

/** \brief Build global stiffness and mass matrices for the eigen system
 *  \details Case where the mass matrix is lumped
 *  @param [in] K Stiffness matrix
 *  @param [in] M Vector containing diagonal mass matrix
 */
    void buildEigen(SkSMatrix<real_t>& K,
                    Vect<real_t>&      M);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   void BodyRHS(const Vect<real_t> &f) { }
   void BoundaryRHS(const Vect<real_t>& f) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
   
 private:

   real_t       _ae, _I1e, _I2e, _h, _mu;
   Vect<real_t> _section, _I1, _I2;
   bool         _bending, _axial, _torsion, _shear, _reduced_integration;
   void set(const Element *el);
   void set(const Side *sd) { }
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
