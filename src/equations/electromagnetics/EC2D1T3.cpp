/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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
#include "io/Fct.h"

namespace OFELI {


EC2D1T3::EC2D1T3()
{
   _equation_name = "Eddy Current";
   _finite_element = "2-D, 3-Node Triangles (P1)";
}


EC2D1T3::EC2D1T3(Mesh& ms)
        : Equation<3,6,2,4>(ms)
{
   if (Equa::_nb_dof!=2)
      throw OFELIException("In EC2D1T3::EC2D1T3(..): Nodes must have "
                            "2 degrees of freedom because of complex-valued formulation.");
   _omega = 1.;
   _volt = 1.;
}


EC2D1T3::EC2D1T3(Mesh&         ms,
                 Vect<real_t>& u)
        : Equation<3,6,2,4>(ms,u)
{
   if (Equa::_nb_dof!=2)
      throw OFELIException("In EC2D1T3::EC2D1T3(..): Nodes must have "
                            "2 degrees of freedom because of complex-valued formulation.");
}


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
   if (_omega_set)
      _omega = _omega_fct(_el_geo.center,0.);
   if (_Mu_set)
      _Mu = _Mu_fct(_el_geo.center,0.);
   if (_sigma_set)
      _sigma = _sigma_fct(_el_geo.center,0.);
   eMat = 0.;
   eRHS = 0.;
}


void EC2D1T3::set(const Side* sd)
{
   _theElement = nullptr, _theSide = sd;
   Line2 ln(sd);
   SideNodeCoordinates();
   _el_geo.center = ln.getCenter();
   _el_geo.length = ln.getLength();
   sMat = 0.;
   sRHS = 0.;
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


void EC2D1T3::Magnetic(real_t coef)
{
   real_t c = 0.5*coef*_Mu*_omega*_el_geo.area*OFELI_SIXTH;
   for (size_t i=1; i<=3; ++i) {
      for (size_t j=1; j<=3; ++j) {
         eMat(2*i,2*j-1) += c;
         eMat(2*i,2*j  ) += c;
      }
      eMat(2*i,2*i-1) += c;
      eMat(2*i,2*i  ) += c;
   }
}


void EC2D1T3::Electric(real_t coef)
{
   real_t d = coef*_el_geo.area/_sigma;
   for (size_t i=1; i<=3; ++i) {
      for (size_t j=1; j<=3; ++j) {
         eMat(2*i-1,2*j-1) += d*(_dSh[i-1],_dSh[j-1]);
         eMat(2*i-1,2*j  ) += d*(_dSh[i-1],_dSh[j-1]);
         eMat(2*i  ,2*j-1) += d*(_dSh[i-1],_dSh[j-1]);
         eMat(2*i  ,2*j  ) += d*(_dSh[i-1],_dSh[j-1]);
      }
   }
}


void EC2D1T3::IntegMF(real_t& vr, real_t& vi)
{
   vr = vi = 0.;
   static real_t c = _Mu*_el_geo.area*OFELI_SIXTH;
   for (size_t i=1; i<=3; ++i) {
      vr += _eu(2*i-1)*c;
      vi += _eu(2*i  )*c;
   }
}


void EC2D1T3::IntegND(const Vect<real_t>& h,
                      real_t&             vr,
                      real_t&             vi)
{
   int c1 = (*_theElement)(1)->getCode(1),
       c2 = (*_theElement)(2)->getCode(1),
       c3 = (*_theElement)(3)->getCode(1);
   Point<real_t> dhr, dhi;
   vr = vi = 0.;
   real_t h1r = h(2*(*_theElement)(1)->n()-1),
          h2r = h(2*(*_theElement)(2)->n()-1),
          h3r = h(2*(*_theElement)(3)->n()-1);
   real_t h1i = h(2*(*_theElement)(1)->n()),
          h2i = h(2*(*_theElement)(2)->n()),
          h3i = h(2*(*_theElement)(3)->n());
   size_t j=0, i1[2], i2[2];
   if (c1==1 && c2==1)
      i1[j] = 1, i2[j] = 2;
   if (c2==1 && c3==1)
      i1[j] = 2, i2[j] = 3;
   if (c3==1 && c1==1)
      i1[j] = 3, i2[j] = 1;
   if (j) {
      dhr = h1r*_dSh[0] + h2r*_dSh[1] + h3r*_dSh[2];
      dhi = h1i*_dSh[0] + h2i*_dSh[1] + h3i*_dSh[2];
      for (size_t i=0; i<j; i++) {
         real_t n1 = _x[i2[i]-1].y - _x[i1[i]-1].y;
         real_t n2 = _x[i1[i]-1].x - _x[i2[i]-1].x;
         vr += (n1*dhr.x + n2*dhr.y)/_sigma;
         vi += (n1*dhi.x + n2*dhi.y)/_sigma;
      }
   }
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
      real_t c = _Mu*_area*OFELI_SIXTH;
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
   Point<real_t> dhr = _eu(1)*_dSh[0] + _eu(3)*_dSh[1] + _eu(5)*_dSh[2], 
                 dhi = _eu(2)*_dSh[0] + _eu(4)*_dSh[1] + _eu(6)*_dSh[2];
   return 0.5/(dhr*dhr+dhi*dhi)/_sigma;
}

} /* namespace OFELI */
