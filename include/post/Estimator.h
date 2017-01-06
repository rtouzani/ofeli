/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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
 */

class Estimator
{

 public:

/*! \enum ESTIMATOR_TYPE
 * Enumerate variable that selects an error estimator for mesh adaptation purposes
 */
enum EstimatorType {
   ESTIM_ZZ      =  0,    /*!< Zhu-Zienckiewicz elementwise estimator     */
   ESTIM_ND_JUMP =  1     /*!< Normal derivative jump sidewise estimator  */
};

/// \brief Default Constructor
    Estimator() { }

/** \brief Constructor using finite element mesh and elementwise estimator
 *  @param [in,out] e Vector that contains once the member function setError is invoked
 *  a posteriori estimator at each element
 */
    Estimator(Mesh&         m,
              Vect<real_t>& e);

/// \brief Destructor
    ~Estimator() { }

    void setType(EstimatorType t=ESTIM_ZZ);

/// \brief Calculate error using Vect solution vector \a u.
    void setError(const Vect<real_t>& u);

/// \brief Return averaged error.
    real_t getAverage() const { return _average; }

/// \brief Return a reference to the finite element mesh.
    Mesh& getMesh() const { return *_mesh; }

    friend ostream& operator<<(ostream& s, const Estimator& r);

private:

   size_t _nb_el, _nb_sd, _nb_nd;
   Mesh *_mesh;
   real_t _average;
   Vect<real_t> *_est;
   Vect<Point<real_t> > _N;
   EstimatorType _est_type;
   void elementT3_ZZ(const Vect<real_t>&   u,
                     Vect<real_t>&         M,
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
