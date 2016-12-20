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

                      Implementation of class WaterPorous2D

  ==============================================================================*/

#include "equations/porous/WaterPorous2D.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"

namespace OFELI {
/*!
 *  \addtogroup OFELI
 *  @{
 */

WaterPorous2D::WaterPorous2D()
              : _nb_nodes(0), _nb_elements(0), _max_iter(50), _verb(1), _toler(1.e-5)
{
   _nb_dof = 1;
   _verb = 1;
   _max_iter = 50;
   _toler = 1.e-5;
   _cw = 1.e-6;
   _Kx = _Ky = 1.;
   _mu = 1.;
}


WaterPorous2D::WaterPorous2D(Mesh&  ms,
                             size_t verb)
              : _nb_nodes(ms.getNbNodes()), _nb_elements(ms.getNbElements()), _max_iter(50), _verb(verb)
{
   _theMesh = &ms;
   _nb_dof = 1;
   _nb_nodes = ms.getNbNodes();
   _nb_elements = ms.getNbElements();
   _max_iter = 50;
   _toler = 1.e-6;
   _cw = 1.e-6;
   _phi.setSize(_nb_nodes);
   _phi = 1.;
   _Kx = _Ky = 1.;
   _mu = 1.;
   if (_verb > 3) {
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
   _Kx.setSize(_nb_elements);
   _Ky.setSize(_nb_elements);
   _Kx = Kx;
   _Ky = Ky;
}


void WaterPorous2D::set(const Element *el)
{
   Init(el);
   eA0 = 0;
   eA1 = 0;
   eRHS = 0;
   Triang3 tr(el);
   _area = tr.getArea();
   _Kxe = _Kx(el->n()), _Kye = _Ky(el->n());
   _dSh(1) = tr.DSh(1), _dSh(2) = tr.DSh(2), _dSh(3) = tr.DSh(3);
}


void WaterPorous2D::set(const Side *sd)
{
   Init(sd);
   sRHS = 0;
   _length = Line2(sd).getLength();
}


void WaterPorous2D::Mass()
{
   real_t c=OFELI_THIRD*_area*_cw*_density;
   for (size_t i=1; i<=3; i++)
      eA1(i,i) = c*_phi(i);
}


void WaterPorous2D::Mobility()
{
   real_t c=_area*_Mw*_density;
   for (size_t i=1; i<=3; i++) {
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += c*(_Kxe*_dSh(i).x*_dSh(j).x +
                        _Kye*_dSh(i).y*_dSh(j).y);
   }
}


void WaterPorous2D::BodyRHS(const Vect<real_t>& bf,
                            int                 opt)
{
   real_t c=OFELI_THIRD*_area;
   if (opt==GLOBAL_ARRAY) {
      for (size_t i=1; i<=3; i++)
         eRHS(i) += c*bf((*_theElement)(i)->n());
   }
   else {
      for (size_t i=1; i<=3; i++)
         eRHS(i) += c*bf(i);
   }
}


void WaterPorous2D::BoundaryRHS(const Vect<real_t>& sf,
                                int                 opt)
{
   real_t c=0.5*_length;
   if (opt==GLOBAL_ARRAY) {
      for (size_t i=1; i<=2; i++)
         sRHS(i) += c*sf((*_theSide)(i)->n());
   }
   else {
      for (size_t i=1; i<=2; i++)
         sRHS(i) += c*sf(i);
   }
}

/*! @} End of Doxygen Groups */
} /* namespace OFELI */
