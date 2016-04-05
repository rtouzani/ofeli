/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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
#include "io/UserData.h"

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

class Beam3DL2 : virtual public Equa_Solid<real_t,2,12,1,1> {

public :

/// \brief Default Constructor
    Beam3DL2() {
      _bending = _axial = _torsion = _shear = true;
      _reduced_integration = false;
    }

/** \brief Constructor using element data
    @param [in] el Pointer to Element
    @param [in] A Section area of the beam
    @param [in] I1 first (x) momentum of inertia
    @param [in] I2 second (y) momentum of inertia
 */
    Beam3DL2(Element* el,
             real_t   A,
             real_t   I1,
             real_t   I2);

/** \brief Constructor for dynamic problems
    @param [in] el Pointer to Element
    @param [in] A Section area of the beam
    @param [in] I1 first (x) momentum of inertia
    @param [in] I2 second (y) momentum of inertia
    @param [in] u Vector containing previous solution (at previous time step)
    @param [in] time Current time value
 */
    Beam3DL2(      Element*      el,
                   real_t        A,
                   real_t        I1,
                   real_t        I2,
             const Vect<real_t>& u,
             const real_t&       time=0);

/** \brief Constructor to determine displacements
 *  \details The unknowns consist in planar and rotational degrees of freedom.
 *  This member function construct a 3-D node vector that gives the displacement vector
 *  at each node.
 *  @param [in] ms Mesh instance
 *  @param [in] u Vector containing the solution vector
 *  @param [out] d Vector containing three components for each node that are <tt>x</tt>,
 *  <tt>y</tt> and <tt>z</tt> displacements.
 */
    Beam3DL2(      Mesh& ms,
             const Vect<real_t>& u,
                   Vect<real_t>& d);

/// \brief Destructor
    ~Beam3DL2() { }

/// \brief Add element lumped Mass contribution to matrix after multiplication by <tt>coef</tt>
    void LMassToLHS(real_t coef=1.);

/// \brief Add element lumped Mass contribution to RHS after multiplication by <tt>coef</tt>
    void LMassToRHS(real_t coef=1.);

/// \brief Add element consistent Mass contribution to matrix after multiplication by 
/// <tt>coef</tt> (not implemented)
    void MassToLHS(real_t coef=1.) { coef = 0; cerr << "Error: Beam3DL2::Mass2LHS not implemented" << endl; }

/// \brief Add element consistent Mass contribution to RHS after multiplication by 
/// <tt>coef</tt> (not implemented)
    void MassToRHS(real_t coef=1.) { coef = 0; cerr << "Error: Beam3DL2::Mass2RHS not implemented" << endl; }

/// \brief Add element stiffness to left hand side
    void Stiffness(real_t coef=1.);

// Add body right-hand side term to right hand side
//    void BodyRHS(UserData<real_t> &ud);

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

/// \brief Return axial force in element
    real_t AxialForce() const;

/// \brief Return shear force in element
    Point<real_t> ShearForce() const;

/// \brief Return bending moment in element
    Point<real_t> BendingMoment() const;

/// \brief Return twisting moment in element
    real_t TwistingMoment() const;

/** \brief Build global stiffness and mass matrices for the eigen system
 *  \details Case where the mass matrix is lumped
 *  @param [in] K Stiffness matrix
 *  @param [in] M Vector containing diagonal mass matrix
 */
    void buildEigen(SkSMatrix<real_t>& K,
                    Vect<real_t>&      M);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
   void BodyRHS(const Vect<real_t> &f, int opt=GLOBAL_ARRAY) { }
   void BoundaryRHS(const Vect<real_t>& f) { }
   void BoundaryRHS(UserData<real_t>& ud) { }
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
   
private:

   real_t        _section, _I1, _I2, _h, _mu;
   Point<real_t> _x[2];
   bool          _bending, _axial, _torsion, _shear, _reduced_integration;
   void set(const Element *el);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
