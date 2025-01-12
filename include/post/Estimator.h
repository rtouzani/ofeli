/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

             Definition of class 'Estimator' for error post-estimation

  ==============================================================================*/


#ifndef __ESTIMATOR_H
#define __ESTIMATOR_H

#include "OFELI_Config.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"

namespace OFELI {

/*! \file Estimator.h
 *  \brief Definition file for class Estimator.
 */

/*! \class Estimator
 *  \ingroup Equation
 *  \brief To calculate an a posteriori estimator of the solution.
 *
 *  \details This class enables calculating an estimator of a solution in order to
 *  evaluate reliability. Estimation uses the so-called Zienkiewicz-Zhu estimator.
 *
 * \author Rachid Touzani
 * \copyright GNU Lesser Public License
 */

class Estimator
{

 public:

/*! \enum EstimatorType
 * Enumerate variable that selects an error estimator for mesh adaptation purposes
 */
enum EstimatorType {
   ESTIM_ZZ      =  0,    /*!< Zhu-Zienckiewicz elementwise estimator     */
   ESTIM_ND_JUMP =  1     /*!< Normal derivative jump sidewise estimator  */
};

/// \brief Default Constructor
    Estimator() { }

/** \brief Constructor using finite element mesh
 *  @param [in] m Mesh instance
 */
    Estimator(Mesh& m);

/// \brief Destructor
    ~Estimator() { }

/** \brief Select type of a posteriori estimator
 *  @param [in] t Type of estimator. It has to be chosen among the enumerated values:
 *  <ul>
 *     <li><tt>ESTIM_ZZ</tt>: The Zhu-Zienckiewicz estimator (Default value)
 *     <li><tt>ESTIM_ND_JUMP</tt>: An estimator based on the jump of normal derivatives of
 *         the solution across mesh sides
 *  </ul>
 */
    void setType(EstimatorType t=ESTIM_ZZ);

/** \brief Provide solution vector in order to determine error index.
 *  @param [in] u Vector containing solution at mesh nodes
 */
    void setSolution(const Vect<real_t>& u);

/** \brief Get vector containing elementwise error index
 *  @param [in,out] e Vector that contains once the member function setError is invoked
 *  a posteriori estimator at each element
 */
    void getElementWiseIndex(Vect<real_t>& e);

/** \brief Get vector containing nodewise error index
 *  @param [in,out] e Vector that contains once the member function setError is invoked
 *  a posteriori estimator at each node
 */
    void getNodeWiseIndex(Vect<real_t>& e);

/** \brief Get vector containing sidewise error index
 *  @param [in,out] e Vector that contains once the member function setError is invoked
 *  a posteriori estimator at each side
 */
    void getSideWiseIndex(Vect<real_t>& e);

/// \brief Return averaged error.
    real_t getAverage() const { return _average; }

/// \brief Return a reference to the finite element mesh.
    Mesh& getMesh() const { return *_mesh; }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    friend ostream& operator<<(ostream& s, const Estimator& r);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

private:

    size_t _nb_el, _nb_sd, _nb_nd, _nb_dof;
   Mesh *_mesh;
   real_t _average;
   Vect<real_t> _el_I, _sd_I, _nd_I;
   Vect<Point<real_t> > _N;
   EstimatorType _est_type;
   void elementT3_ZZ(const Vect<real_t>&   u,
                     Vect<Point<real_t> >& b);
   void elementT3_ND_JUMP(const Vect<real_t>& u);
};


/// \fn ostream& operator<<(ostream& s, const Estimator &r)
/// \ingroup Equation
/// \brief Output estimator vector in output stream
ostream& operator<<(ostream&         s,
                    const Estimator& r);

} /* namespace OFELI */

#endif
