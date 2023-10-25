/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                     Implementation of class Pres2DT3

  ==============================================================================*/


#include "equations/acoustics/Pres2DT3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/LocalMatrix_impl.h"
#include "linear_algebra/LocalVect_impl.h"
#include "equations/Equation_impl.h"

namespace OFELI {


Pres2DT3::Pres2DT3()
{
   _equation_name = "Acoustics";
   _finite_element = "2-D, 3-Node Triangles (P1)";
}


Pres2DT3::Pres2DT3(Mesh& ms)
         : Equation<3,3,2,2>(ms)
{
   _equation_name = "Acoustics";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
}


Pres2DT3::Pres2DT3(Mesh&         ms,
                   Vect<real_t>& u)
       : Equation<3,3,2,2>(ms,u)
{
   _equation_name = "Acoustics";
   _finite_element = "2-D, 3-Node Triangles (P1)";
   setMatrixType(SPARSE);
   setSolver(CG_SOLVER,DILU_PREC);
}


Pres2DT3::~Pres2DT3() { }


void Pres2DT3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   _dSh = tr.DSh();
   _el_geo.size = 2*tr.getCircumRadius();
   ElementNodeCoordinates();
   ElementNodeVector(*_u,_eu);
   eA0 = 0, eA1 = 0, eA2 = 0;
   eRHS = 0;
}


void Pres2DT3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _el_geo.center = ln.getCenter();
   _el_geo.length = ln.getLength();
   sMat = 0, sRHS = 0;
}


void Pres2DT3::setInput(EqDataType    opt,
                        Vect<real_t>& u)
{
   Equa::setInput(opt,u);
}


void Pres2DT3::LMass(real_t coef)
{
   real_t c=OFELI_THIRD*_el_geo.area/(_speed*_speed);
   eA2(1,1) += c;
   eA2(2,2) += c;
   eA2(3,3) += c;
   eMat += eA2;
}


void Pres2DT3::Mass(real_t coef)
{
   real_t c=OFELI_SIXTH*_el_geo.area*coef/(_speed*_speed);
   real_t d=0.5*c;
   eA2(1,1) += c; eA2(2,2) += c; eA2(3,3) += c;
   eA2(1,2) += d; eA2(2,1) += d; eA2(1,3) += d;
   eA2(3,1) += d; eA2(2,3) += d; eA2(3,2) += d;
}


void Pres2DT3::Diffusion(real_t coef)
{
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += coef*_el_geo.area*(_dSh[i-1],_dSh[j-1]);
   eMat += eA0;
}


void Pres2DT3::BodyRHS(real_t bf)
{
   real_t c=OFELI_THIRD*_el_geo.area*bf;
   eRHS(1) += c;
   eRHS(2) += c;
   eRHS(3) += c;
}


void Pres2DT3::BodyRHS(const Vect<real_t>& f)
{
   real_t c = OFELI_THIRD*_el_geo.area;
   for (size_t i=1; i<=3; i++)
      eRHS(i) += c*f((*_theElement)(i)->n());
}


void Pres2DT3::BoundaryRHS(real_t flux)
{
   real_t c = 0.5*_el_geo.length*flux;
   for (size_t i=1; i<=2; i++)
       sRHS(i) += c;
}


void Pres2DT3::BoundaryRHS(const Vect<real_t>& f)
{
   if (_theSide->getCode(1)>0) {
      real_t c = 0.5*_el_geo.length;
      if (f.getDOFType()==NODE_DOF) {
         sRHS(1) += c*f((*_theSide)(1)->n());
         sRHS(2) += c*f((*_theSide)(2)->n());
      }
      else if (f.getDOFType()==SIDE_DOF) {
         sRHS(1) += c*f(_theSide->n());
         sRHS(2) += c*f(_theSide->n());
      }
   }
}


Point<real_t> &Pres2DT3::Flux() const
{
   _f = _eu(1)*_dSh[0] + _eu(2)*_dSh[1] + _eu(3)*_dSh[2];
   return _f;
}


void Pres2DT3::Grad(Vect<Point<real_t> >& g)
{
   MESH_EL {
      set(the_element);
      g(element_label) = _eu(1)*_dSh[0] + _eu(2)*_dSh[1] + _eu(3)*_dSh[2];
   }
}

} /* namespace OFELI */