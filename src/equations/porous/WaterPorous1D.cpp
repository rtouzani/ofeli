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

                      Implementation of class WaterPorous1D

  ==============================================================================*/

#include "equations/porous/WaterPorous1D.h"
#include "shape_functions/Line2.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

WaterPorous1D::WaterPorous1D()
              : _nb_nodes(0), _nb_elements(0), _verb(1)
{
   _nb_dof = 1;
   _cw = 1.e-6;
   _Kx = 1.;
   _mu = 1.;
}


WaterPorous1D::WaterPorous1D(Mesh&  ms,
                             size_t verb)
              : _nb_nodes(ms.getNbNodes()), _nb_elements(ms.getNbElements()),
                _verb(verb)
{
   _theMesh = &ms;
   _nb_dof = 1;
   _cw = 1.e-6;
   _phi.setSize(_nb_nodes);
   _phi = 1.;
   _Kx = 1.;
   _mu = 1.;
   if (_verb > 3) {
      cout << "Number of nodes:\t\t" << _nb_nodes << endl;
      cout << "Number of equations:\t\t" << _theMesh->getNbEq() << endl;
   }
}


WaterPorous1D::~WaterPorous1D() { }


void WaterPorous1D::setCoef(real_t cw,
                            real_t phi,
                            real_t rho,
                            real_t K,
                            real_t mu)
{
   _cw = cw;
   _phi.setSize(_nb_nodes);
   _phi = phi;
   _density = rho;
   _Mw = 1./mu;
   _Kx.setSize(_nb_elements);
   _Kx = K;
}


void WaterPorous1D::set(const Element *el)
{
   Init(el);
   eA0 = 0;
   eA1 = 0;
   eRHS = 0;
   Line2 ln(el);
   _length = ln.getLength();
   _Ke = _Kx(el->n());
   _dSh(1) = ln.DSh(1), _dSh(2) = ln.DSh(2);
}


void WaterPorous1D::Mass()
{
   real_t c=0.5*_length*_cw*_density;
   eA1(1,1) = c*_phi(1);
   eA1(2,2) = c*_phi(1);
}


void WaterPorous1D::Mobility()
{
   real_t c=_length*_Mw*_density;
   for (size_t i=1; i<=2; i++) {
      for (size_t j=1; j<=2; j++)
         eA0(i,j) += c*(_Ke*_dSh(i).x*_dSh(j).x);
   }
}


void WaterPorous1D::BodyRHS(const Vect<real_t>& bf,
                            int                 opt)
{
   real_t c=0.5*_length;
   if (opt==GLOBAL_ARRAY) {
      eRHS(1) += c*bf((*_theElement)(1)->n());
      eRHS(2) += c*bf((*_theElement)(2)->n());
   }
   else {
      eRHS(1) += c*bf(1);
      eRHS(2) += c*bf(2);
   }
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */
