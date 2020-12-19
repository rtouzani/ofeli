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

                         Implementation of class Muscl1D

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/


#include "equations/cl/Muscl1D.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {


Muscl1D::Muscl1D(Mesh& m) : Muscl(m)
{
   _theMesh = &m;
   _BBl.setSize(_nb_elements,2);        // ||B_j - B_i||
   _BQl.setSize(_nb_elements,2);
   _Lrate.setSize(_nb_sides);           // Ratio surface/length of left side
   _Rrate.setSize(_nb_sides);           // Ratio surface/length of right side
   Initialize();
   mesh_analyze();
   _method = FIRST_ORDER_METHOD;
   _limiter = MINMOD_LIMITER;
}


void Muscl1D::Initialize()
{
   MESH_SD {
      the_element = the_side->getNeighborElement(1);
      _Lrate.set(side_label,1./Line2(the_element).getLength());
      _Rrate.set(side_label,0.);  // default value if on boundary
      if (the_side->isOnBoundary()==false)
         _Rrate.set(side_label,1./Line2(the_side->getNeighborElement(2)).getLength());
   }

   MESH_EL {
      size_t n = element_label;
      if (The_element.isOnBoundary() == false) {
         Side *sd1=The_element.getPtrSide(1), *sd2=The_element.getPtrSide(2);
         Element *el1=sd1->getNeighborElement(1), *el2=sd2->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);
         Line2 ln(the_element), ln1(el1), ln2(el2);

//       Compute BiBj: distance center to center
         _BBl(n,1) = fabs(ln1.getCenter().x-ln.getCenter().x);
         _BBl(n,2) = fabs(ln2.getCenter().x-ln.getCenter().x);
//       Compute BiQj
         _BQl(n,1) = fabs(The_element(1)->getX()-ln.getCenter().x);
         _BQl(n,2) = fabs(The_element(2)->getY()-ln.getCenter().x);
      }
   }
}


void Muscl1D::FirstOrder(const Vect<real_t>& U,
                         Vect<real_t>&       LU,
                         Vect<real_t>&       RU,
                         size_t              dof)
{
   Element *Lel, *Rel;
   MESH_SD {
      size_t s = side_label;
      Lel = The_side.getNeighborElement(1);
      (The_side.isOnBoundary()) ? Rel = Lel : Rel = The_side.getNeighborElement(2);
      LU(s,dof) = U(Lel->n(),dof);
      RU(s,dof) = U(Rel->n(),dof);
   }
}


void Muscl1D::SecondOrder(const Vect<real_t>& U,
                          Vect<real_t>&       LU,
                          Vect<real_t>&       RU,
                          size_t              dof)
{
   MESH_EL {
      size_t n = element_label;
      Side *sd1 = The_element.getPtrSide(1); size_t ns1 = sd1->n();
      Side *sd2 = The_element.getPtrSide(2); size_t ns2 = sd2->n();
      real_t val1=U(n,dof), val2=U(n,dof);
      if (The_element.isOnBoundary()) {
         if (sd1->isOnBoundary())
            RU(ns1,dof) = U(n,dof);
         if (sd2->isOnBoundary())
            RU(ns2,dof) = U(n,dof);
      }
      else {  // Inner element: Second order method
         Element *el1=sd1->getNeighborElement(1), *el2=sd2->getNeighborElement(1);
         if (el1==the_element)
            el1 = sd1->getNeighborElement(2);
         el2 = sd2->getNeighborElement(1);
         if (el2==the_element)
            el2 = sd2->getNeighborElement(2);

//       compute slope_plus
         real_t SlopeQForward1 = (U(el1->n(),dof) - U(n,dof))/_BBl(n,1),
                SlopeQForward2 = (U(el2->n(),dof) - U(n,dof))/_BBl(n,2);

//       compute slope_minus
         real_t SlopeQBackward1 = -SlopeQForward2, SlopeQBackward2 = -SlopeQForward1;

//       compute limited slope
         real_t SlopeQ1 = minmod(SlopeQBackward1,SlopeQForward1)*SlopeQForward1;
         real_t SlopeQ2 = minmod(SlopeQBackward2,SlopeQForward2)*SlopeQForward2;
         val1 = SlopeQ1*_BQl(n,1) + U(n,dof);
         val2 = SlopeQ2*_BQl(n,2) + U(n,dof);
      }
      (sd1->getNeighborElement(1)==the_element) ? LU(ns1,dof) = val1 : RU(ns1,dof) = val1;
      (sd2->getNeighborElement(1)==the_element) ? LU(ns2,dof) = val2 : RU(ns2,dof) = val2;
   }
}


void Muscl1D::mesh_analyze()
{
   size_t nb_el = 0;
   real_t aux, tau;
   _MinimumLength = 1./OFELI_EPSMCH;
   _MaximumLength = _MeanLength = 0.;
   _taulim = 2.;
   MESH_EL {
      nb_el++;
      aux = Line2(the_element).getLength();
      if (_MinimumLength>aux)
         _MinimumLength = aux;
      if (_MaximumLength<aux)
         _MaximumLength = aux;
      _MeanLength += aux;
      if (!the_element->getPtrSide(1)->isOnBoundary() && !the_element->getPtrSide(1)->isOnBoundary()) {
         size_t n = element_label;
         tau = _BBl(n,1)/_BQl(n,1);
         if (tau<_taulim)
            _taulim = tau;
         tau = _BBl(n,2)/_BQl(n,2);
         if (tau<_taulim)
            _taulim = tau;
      }
   }
   _MeanLength /= nb_el;
}

} /* namespace OFELI */
