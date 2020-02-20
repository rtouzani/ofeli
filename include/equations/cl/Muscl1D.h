/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                          Definition of Class Muscl1D
                    Class for Muscl finite volumes for 1-D problems

  ==============================================================================*/

#ifndef __MUSCL1D_H
#define __MUSCL1D_H

#include "equations/cl/Muscl.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

/*! \file Muscl1D.h
 *  \brief Definition file for class Muscl1D.
 */

template<class T_,size_t N_> class LocalVect;
template<class T_> class Vect;


/*! \class Muscl1D
 *  \ingroup ConservationLaws
 *  \brief Class for 1-D hyperbolic solvers with Muscl scheme.
 *
 * \author S. Clain, V. Clauzon
 * \copyright GNU Lesser Public License
 */


class Muscl1D : public Muscl {

 public:

   using Muscl::_nb_sides;
   using Muscl::_nb_elements;

/// \brief Constructor using mesh instance
    Muscl1D(Mesh& m);

/// \brief Destructor
    ~Muscl1D() { }

/// \brief Return mean length
    real_t getMeanLength() const { return _MeanLength; }

/// \brief Return maximal length
    real_t getMaximumLength() const { return _MaximumLength; }

/// \brief Return mimal length
    real_t getMinimumLength() const { return _MinimumLength; }

/// \brief Return mean length
    real_t getTauLim() const { return _taulim; }

/// \brief Output mesh information
    void print_mesh_stat();

 protected:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    Vect<real_t> _BBl, _BQl;
    Vect<real_t> _Lrate, _Rrate;
    void Initialize();

/// Implementation of the virtual function defined in Muscl.h
    void FirstOrder(const Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof);
    void SecondOrder(const Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof);
    void MultiSlopeM(const Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof)
    {
       SecondOrder(U,LU,RU,dof);
    }

    void MultiSlopeQ(Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof) 
    {
       SecondOrder(U,LU,RU,dof);
    }

    void GradientQ(Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof)
    {
       SecondOrder(U,LU,RU,dof);
    }

    void GradientM(Vect<real_t> &U, Vect<real_t> &LU, Vect<real_t> &RU, size_t dof)
    {
       SecondOrder(U,LU,RU,dof);
    }

//  mesh characteristics
    real_t _MinimumLength, _MaximumLength, _MeanLength, _taulim;
    void mesh_analyze();
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
