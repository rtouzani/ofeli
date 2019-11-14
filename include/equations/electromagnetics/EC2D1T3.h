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

   using Equation<complex_t,3,3,2,2>::_A;
   using Equation<complex_t,3,3,2,2>::_b;

/// \brief Default constructor
    EC2D1T3();

/** \brief Constructor using mesh
 *  @param [in] ms Mesh instance
 */
    EC2D1T3(Mesh& ms);

/** \brief Constructor using mesh and solution vector
 *  @param [in] ms Mesh instance
 *  @param [in,out] u Reference to solution vector instance
 */
    EC2D1T3(Mesh&            ms,
            Vect<complex_t>& u);

/// \brief Destructor
    ~EC2D1T3();

/** \brief Define data for equation
 *  @param [in] omega Angular frequency
 *  @param [in] volt Voltage
 */
    void setData(real_t omega,
                 real_t volt);

/** \brief Build the linear system of equations
 *  \details Before using this function, one must have properly selected 
 *  appropriate options for:
 *  <ul>
 *     <li>The choice of a steady state or transient analysis. By default, the analysis is stationary
 *     <li>In the case of transient analysis, the choice of a time integration scheme
 *         and a lumped or consistent capacity matrix. If transient analysis is chosen, the lumped
 *         capacity matrix option is chosen by default, and the implicit Euler scheme is used
 *         by default for time integration.
 *  </ul>
 */
    void build()
    {
       mesh_elements(*_theMesh) {
          set(the_element);
          Magnetic(_omega,1.);
          Electric();
          AbsEqua<complex_t>::_A->Assembly(The_element,eMat.get());
       }

       mesh_nodes(*_theMesh) {
          int m = the_node->getDOF(1);
          if (the_node->getCode(1)==1) {
             _A->set(m,m,(*_A)(m,m)*VLG);
             (*_b)[m-1] = (*_A)(m,m)*_current;
          }
          if (the_node->getCode(1)==2) {
             _A->set(m,m,(*_A)(m,m)*VLG);
             (*_b)[m-1] = 0;
          }
       }
    }

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
 *  @param [in] h Vect instance containing magnetic field at nodes
 *  @note This member function is to be called within each element, it detects
 *  boundary sides as the ones with nonzero code
 */
    complex_t IntegND(const Vect<complex_t>& h);

/// \brief Add contribution to vacuum area calculation
    real_t VacuumArea();

 protected:

    void set(const Element* el);
    void set(const Side* el);

 private:
    real_t _omega, _volt, _current;

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
