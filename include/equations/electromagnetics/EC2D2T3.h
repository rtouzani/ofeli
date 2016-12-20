/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Libry

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

                         Definition of class EC2D2T3
     for Eddy Current Problems in Two-Dimensions with a vector magnetic field
                           using the 3-Node triangle

  ==============================================================================*/


#ifndef __EC2D2T3_H
#define __EC2D2T3_H


#include "equations/electromagnetics/Equa_Electromagnetics.h"
#include "linear_algebra/SpMatrix.h"
#include "linear_algebra/SkMatrix.h"

#define GAUSS2 0.577350269189625764509148780501957

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */


class EC2D2T3 : virtual public Equa_Electromagnetics<real_t,3,6,2,4>
{

 public :

/// \brief Default Constructor
    EC2D2T3() {
    }

/// \brief Constructor using element data
    EC2D2T3(const Element* el);

/// \brief Constructor using one side data
    EC2D2T3(const Side* sd);

/// \brief Constructor using two side data
    EC2D2T3(const Side* sd1,
            const Side* sd2);

/// \brief Destructor
    ~EC2D2T3()
    {  }

/// \brief Compute Contribution to Right-Hand Side
    void RHS(real_t coef=1.);

/// \brief Compute Finite Element Diagonal Block
    void FEBlock(real_t omega);

/// \brief Compute boundary element blocks
    void BEBlocks(size_t            n1,
                  size_t            n2,
                  SpMatrix<real_t>& L,
                  SpMatrix<real_t>& U,
                  SpMatrix<real_t>& D);

/// \brief Compute constant to multiply by potential
    complex_t Constant(      real_t        omega,
                       const Vect<real_t>& u,
                             complex_t&    I);

/// \brief Compute magnetic pressure in element
    real_t MagneticPressure(const Vect<real_t>& u);

 private:

   real_t        _ll1, _ll2, _a3, _area;
   LocalVect<Point<real_t>,3> _N;
   Point<real_t> _N1, _N2, _M1, _M2;
   size_t        _ns, _nt, _i1, _j1, _i2, _j2;
   Triang3       _tr;
   size_t _log_det(const Point<real_t> &ck, const Point<real_t> &cl);
   complex_t _ablog(size_t det, complex_t a, complex_t b, real_t t);
   void _Lkl(const Point<real_t> &uk, const Point<real_t> &vk,
             const Point<real_t> &ul, const Point<real_t> &vl,
             real_t &d1, real_t &d2);
   void _Dkl(const Point<real_t> &uk, const Point<real_t> &vk,
             const Point<real_t> &ul, const Point<real_t> &vl, real_t &d);
   void set(const Element *el);
   void set(const Side *sd);
};

/*! @} End of Doxygen Groups */
} /* namespace OFELI */

#endif
