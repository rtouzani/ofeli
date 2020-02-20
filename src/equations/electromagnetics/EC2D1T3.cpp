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

                        Implementation of class EC2D1T3

  ==============================================================================*/


#include "equations/electromagnetics/EC2D1T3.h"
#include "shape_functions/Triang3.h"
#include "shape_functions/Line2.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {


EC2D1T3::EC2D1T3()
{
   _equation_name = "Eddy Current";
   _finite_element = "2-D, 3-Node Triangles (P1)";
}


EC2D1T3::EC2D1T3(Mesh& ms)
        : Equation<complex_t,3,3,2,2>(ms)
{
   _omega = 1.;
   _volt = 1.;
}


EC2D1T3::EC2D1T3(Mesh&            ms,
                 Vect<complex_t>& u)
        : Equation<complex_t,3,3,2,2>(ms,u)
{ }


EC2D1T3::~EC2D1T3() { }


void EC2D1T3::setData(real_t omega,
                      real_t volt)
{
   _omega = omega;
   _volt = volt;
}


void EC2D1T3::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   Triang3 tr(_theElement);
   _el_geo.area = tr.getArea();
   _el_geo.center = tr.getCenter();
   ElementNodeCoordinates();
   _dSh = tr.DSh();
   ElementNodeVector(*_u,_eu);
   eMat = complex_t(0.);
   eRHS = complex_t(0.);
}


void EC2D1T3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _el_geo.center = ln.getCenter();
   _el_geo.length = ln.getLength();
   sMat = complex_t(0.);
   sRHS = complex_t(0.);
}

/*
int EC2D1T3::run()
{
   MESH_EL {
     set(theElement);
     Magnetic(1.,_omega);
     Electric();
     _A->Assembly(theElement,eMat.get());
     _b.Assembly(theElement,eRHS.get());
   }
   Prec<complex_t> p(*_Am,ILU_PREC);
   int nb_it = BiCGStab(*_A,p,_b,_uu,_max_it,_toler,0);
   if (_verbose>1)
      cout << "Nb. of iterations for electromagnetics: " << nb_it << endl;
   insertBC(_uu,*_u);

      eqH.NormalizedMF(T1,h);
      eqH.UpdateMF(h);
   return nb_it;
}*/


void EC2D1T3::Magnetic(real_t coef,
                       real_t omega)
{
   complex_t c = 0.5*coef*_mu*omega*_el_geo.area*OFELI_SIXTH*OFELI_IMAG;
   complex_t d = coef*_mu*omega*_el_geo.area*OFELI_SIXTH*OFELI_IMAG;
   eMat(1,1) += d; eMat(2,2) += d; eMat(3,3) += d;
   eMat(1,2) += c; eMat(2,1) += c; eMat(1,3) += c;
   eMat(3,1) += c; eMat(2,3) += c; eMat(3,2) += c;
}


void EC2D1T3::Electric(real_t coef)
{
   static real_t d = coef*_rho*_el_geo.area;
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eMat(i,j) += d*(_dSh[i-1],_dSh[j-1]);
}


complex_t EC2D1T3::IntegMF()
{
   static real_t c = _mu*_el_geo.area*OFELI_SIXTH;
   complex_t x=0;
   for (size_t i=1; i<=3; i++)
      x += _eu(i)*c;
   return x;
}


complex_t EC2D1T3::IntegND(const Vect<complex_t>& h)
{
   int c1 = (*_theElement)(1)->getCode(1),
       c2 = (*_theElement)(2)->getCode(1),
       c3 = (*_theElement)(3)->getCode(1);
   complex_t dhx, dhy, xx=0;
   complex_t h1 = h((*_theElement)(1)->n()),
             h2 = h((*_theElement)(2)->n()),
             h3 = h((*_theElement)(3)->n());
   size_t j=0, i1[2], i2[2];
   if (c1==1 && c2==1)
      i1[j] = 1, i2[j] = 2;
   if (c2==1 && c3==1)
      i1[j] = 2, i2[j] = 3;
   if (c3==1 && c1==1)
      i1[j] = 3, i2[j] = 1;
   if (j) {
      dhx = h1*_dSh[0].x + h2*_dSh[1].x + h3*_dSh[2].x;
      dhy = h1*_dSh[0].y + h2*_dSh[1].y + h3*_dSh[2].y;
      for (size_t i=0; i<j; i++) {
         real_t n1 = _x[i2[i]-1].y - _x[i1[i]-1].y;
         real_t n2 = _x[i1[i]-1].x - _x[i2[i]-1].x;
         xx += _rho*(n1*dhx + n2*dhy);
      }
   }
   return xx;
}


real_t EC2D1T3::VacuumArea()
{
   return 0.25*((_x[0].x+_x[1].x)*(_x[1].y-_x[0].y)+(_x[0].y+_x[1].y)*(_x[0].x-_x[1].x));
}

/*
void EC2D1T3::UpdateMF()
{
   complex_t a(0,omega);
   complex_t yy;
   MeshElements(ms0) {
      set(theElement);
      yy += IntegND(b0);
   }
   complex_t xx=0;
   MeshElements(ms1) {
      set(theElement);
      real_t c = _mu*_area*OFELI_SIXTH;
      for (size_t i=1; i<=3; i++)
         xx += c*b1(theElement->getNodeLabel(i));
   }

// Calculate area of vacuum
   real_t aa=0.;
   MeshSides(ms0) {
      set(theSide);
      aa += VacuumArea();
   }
   MeshSides(ms1) {
      set(theSide);
      aa -= VacuumArea();
   }
   size_t i;
   for (i=0; i<b0.getSize(); i++)
      b0[i] *= volt/(a*(xx+MU0*aa)+yy);
   for (i=0; i<b1.getSize(); i++)
      b1[i] *= volt/(a*(xx+MU0*aa)+yy);
}*/


real_t EC2D1T3::Joule()
{
   Point<real_t> dhr = _eu(1).real()*_dSh[0] + _eu(2).real()*_dSh[1] + _eu(3).real()*_dSh[2], 
                 dhi = _eu(1).imag()*_dSh[0] + _eu(2).imag()*_dSh[1] + _eu(3).imag()*_dSh[2];
   return 0.5*_rho/(dhr*dhr+dhi*dhi);
}

} /* namespace OFELI */
