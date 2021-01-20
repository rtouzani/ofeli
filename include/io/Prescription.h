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

         Definition of class 'Prescription' to read data for prescription

  ==============================================================================*/

#ifndef __PRESCRIPTION_H
#define __PRESCRIPTION_H

#include <string>
#include <fstream>
#include <vector>
#include "OFELI_Config.h"
#include "io/Fct.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Prescription.h
 *  \brief Definition file for class Prescription.
 */

template <class T_> class Vect;
class Mesh;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct PrescriptionPar {
   size_t dof;
   EqDataType type;
   int code;
   real_t x, y, z, t;
   std::string fct;
   bool bx, by, bz, bt;
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** \class Prescription
 *  \ingroup IO
 *  \brief To prescribe various types of data by an algebraic expression.
 *  Data may consist in boundary conditions, forces, tractions, fluxes, initial condition.
 *  All these data types can be defined through an enumerated variable.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Prescription
{

 public:

/// Default constructor
    Prescription();

/** \brief Constructor that gives an instance of class Mesh and the data file name.
 *  \details It reads parameters in Prescription Format from this file.
 *  @param [in] mesh Mesh instance
 *  @param [in] file Name of Prescription file
 */
    Prescription(Mesh&              mesh,
                 const std::string& file);

/// Destructor
    ~Prescription();

/** Read data in the given file and stores in a Vect instance for a chosen DOF.
 *  The input value type determines the type of data to read.
 *  @param [in] type Type of data to seek. To choose among the enumerated values:
 *  <ul>
 *    <li><tt>BOUNDARY_CONDITION</tt>: Read values for (Dirichlet) boundary conditions
 *    <li><tt>BOUNDARY_FORCE</tt>: Read values for boundary force (Neumann boundary condition).\n
 *       The values <tt>TRACTION</tt> and <tt>FLUX</tt> have the same effect.
 *    <li><tt>BODY_FORCE</tt>: Read values for body (or volume) forces.\n
 *       The value <tt>SOURCE</tt> has the same effect.
 *    <li><tt>POINT_FORCE</tt>: Read values for pointwise forces
 *    <li><tt>CONTACT_DISTANCE</tt>: Read values for contact distance (for contact mechanics)
 *    <li><tt>INITIAL_FIELD</tt>: Read values for initial solution
 *    <li><tt>SOLUTION</tt>: Read values for a solution vector
 *  @param [in,out] v Vect instance that is instantiated on input and filled on output
 *  @param [in] time Value of time for which data is read [Default: <tt>0</tt>].
 *  @param [in] dof DOF to store (Default is <tt>0</tt>: All DOFs are chosen).
 */
    int get(EqDataType    type,
            Vect<real_t>& v,
            real_t        time=0,
            size_t        dof=0);

 private:

   std::ifstream *_if;
   Mesh *_theMesh;
   Vect<real_t> *_v;
   real_t _time;
   bool pforce, initial, bc, force, flux;
   std::string _file;
   std::vector<PrescriptionPar> _p;
   Fct _theFct;
   int Type(int type);
   void get_point_force(size_t k);
   void get_point_force(size_t k, size_t dof);
   void get_boundary_condition(size_t k, size_t dof);
   void get_boundary_condition(size_t k);
   void get_vector(size_t k, size_t dof);
   void get_vector(size_t k);
   void get_boundary_force(size_t k, size_t dof);
   void get_boundary_force(size_t k);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
