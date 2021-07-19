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

                      Implementation of class WaterPorous2D

  ==============================================================================*/

#include "equations/porous/WaterPorous2D.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

WaterPorous2D::WaterPorous2D()
{
   _nb_dof = 1;
   _cw = 1.e-6;
   _Kx = _Ky = 1.;
   _mu = 1.;
}


WaterPorous2D::WaterPorous2D(Mesh&  ms)
              : Equation<3,3,2,2>(ms)
{
   _theMesh = &ms;
   _nb_dof = 1;
   _cw = 1.e-6;
   _phi.setSize(_nb_nodes);
   _phi = 1.;
   _Kx = _Ky = 1.;
   _mu = 1.;
   if (Verbosity>3) {
      cout << "Number of nodes:\t\t" << _nb_nodes << endl;
      cout << "Number of equations:\t\t" << _theMesh->getNbEq() << endl;
   }
}


WaterPorous2D::~WaterPorous2D() { }


void WaterPorous2D::setCoef(real_t cw,
                            real_t phi,
                            real_t rho,
                            real_t Kx,
                            real_t Ky,
                            real_t mu)
{
   _cw = cw;
   _phi.setSize(_nb_nodes);
   _phi = phi;
   _density = rho;
   _Mw = 1./mu;
   _Kx.setSize(_nb_el);
   _Ky.setSize(_nb_el);
   _Kx = Kx;
   _Ky = Ky;
}


void WaterPorous2D::set(const Element *el)
{
  _theElement = el, _theSide = nullptr;
   eA0 = 0;
   eA1 = 0;
   eRHS = 0;
   Triang3 tr(el);
   _el_geo.area = tr.getArea();
   _Kxe = _Kx(el->n()), _Kye = _Ky(el->n());
   _dSh = tr.DSh();
}


void WaterPorous2D::set(const Side *sd)
{
   _theSide = sd, _theElement = nullptr;
   sRHS = 0;
   _el_geo.length = Line2(sd).getLength();
}


void WaterPorous2D::Mass()
{
   real_t c=OFELI_THIRD*_el_geo.area*_cw*_density;
   for (size_t i=1; i<=3; i++)
      eA1(i,i) = c*_phi(i);
}


void WaterPorous2D::Mobility()
{
   real_t c=_el_geo.area*_Mw*_density;
   for (size_t i=1; i<=3; i++) {
      for (size_t j=1; j<=3; j++) {
         eA0(i,j) += c*(_Kxe*_dSh[i-1].x*_dSh[j-1].x +
                        _Kye*_dSh[i-1].y*_dSh[j-1].y);
      }
   }
}


void WaterPorous2D::BodyRHS(const Vect<real_t>& bf)
{
   for (size_t i=1; i<=3; i++)
      eRHS(i) += OFELI_THIRD*_el_geo.area*bf((*_theElement)(i)->n());
}


void WaterPorous2D::BoundaryRHS(const Vect<real_t>& sf)
{
   sRHS(1) += 0.5*_el_geo.length*sf((*_theSide)(1)->n());
   sRHS(2) += 0.5*_el_geo.length*sf((*_theSide)(2)->n());
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */
