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

                      Implementation of class WaterPorous1D

  ==============================================================================*/

#include "equations/porous/WaterPorous1D.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

WaterPorous1D::WaterPorous1D()
{
   _nb_dof = 1;
   _cw = 1.e-6;
   _Kx = 1.;
   _mu = 1.;
}


WaterPorous1D::WaterPorous1D(Mesh&  ms)
              : Equation<real_t,2,2,1,1>(ms)
{
   _nb_dof = 1;
   _cw = 1.e-6;
   _phi.setSize(_nb_nodes);
   _phi = 1.;
   _Kx = 1.;
   _mu = 1.;
   if (Verbosity>3) {
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
   _Kx.setSize(_nb_el);
   _Kx = K;
}


void WaterPorous1D::set(const Element *el)
{
   _theElement = el, _theSide = nullptr;
   eA0 = 0;
   eA1 = 0;
   eRHS = 0;
   Line2 ln(el);
   _el_geo.length = ln.getLength();
   _Ke = _Kx(el->n());
   _dSh = ln.DSh();
}


void WaterPorous1D::Mass()
{
   real_t c=0.5*_el_geo.length*_cw*_density;
   eA1(1,1) = c*_phi(1);
   eA1(2,2) = c*_phi(1);
}


void WaterPorous1D::Mobility()
{
   real_t c=_el_geo.length*_Mw*_density;
   for (size_t i=1; i<=2; i++) {
      for (size_t j=1; j<=2; j++)
         eA0(i,j) += c*(_Ke*_dSh[i-1].x*_dSh[j-1].x);
   }
}


void WaterPorous1D::BodyRHS(const Vect<real_t>& bf)
{
   eRHS(1) += 0.5*_el_geo.length*bf((*_theElement)(1)->n());
   eRHS(2) += 0.5*_el_geo.length*bf((*_theElement)(2)->n());
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */
