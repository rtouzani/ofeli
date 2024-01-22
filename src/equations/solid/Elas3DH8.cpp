/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                         Implementation of class ELAS3DH8
  for 3-D Linear Elasticity Equations using 8-node hexahedral Q1 finite element

  ==============================================================================*/


#include "equations/solid/Elas3DH8.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Elas3DH8::Elas3DH8(Mesh& ms)
         : Equation<8,24,4,12>(ms), _hexa(nullptr), _quad(nullptr)
{ }


Elas3DH8::~Elas3DH8()
{
   if (_hexa!=nullptr)
      delete _hexa; 
   if (_quad!=nullptr)
      delete _quad;
}


void Elas3DH8::set(const Element* el)
{
   _theElement = el, _theSide = nullptr;
   setMaterial();
   _hexa = new Hexa8(el);
   if (_quad != nullptr)
      delete _quad, _quad = nullptr;
   _hexa->atGauss(2,_dSh,_wg);
   _hexa->atGauss(2,_sh,_wg);
   _el_geo.center = _hexa->getCenter();
   ElementNodeCoordinates();
   _xl[0].x = _xl[3].x = _xl[4].x = _xl[7].x = -1.;
   _xl[1].x = _xl[2].x = _xl[5].x = _xl[6].x =  1.;
   _xl[0].y = _xl[1].y = _xl[4].y = _xl[5].y = -1.;
   _xl[2].y = _xl[3].y = _xl[6].y = _xl[7].y =  1.;
   _xl[0].z = _xl[1].z = _xl[2].z = _xl[3].z = -1.;
   _xl[4].z = _xl[5].z = _xl[6].z = _xl[7].z =  1.;
   ElementNodeVector(*_u,_eu);
   if (Equa::_bf!=nullptr)
      ElementNodeVector(*_bf,_ebf);
   if (_rho_set)
      _rho = _rho_fct(_el_geo.center,_TimeInt.time);
   if (_young_set)
      _young = _young_fct(_el_geo.center,_TimeInt.time);
   if (_poisson_set)
      _poisson = _poisson_fct(_el_geo.center,_TimeInt.time);
   eA0 = 0; eA1 = 0; eA2 = 0;
   eRHS = 0;
}


void Elas3DH8::set(const Side* sd)
{
   _theSide = sd, _theElement = nullptr;
   _quad = new Quad4(sd);
   if (_hexa != nullptr)
      delete _hexa, _hexa = nullptr;
   SideNodeCoordinates();
   SideNodeVector(*_u,_su);
   if (Equa::_sf!=nullptr)
      SideVector(*_sf,_ssf);
   sA0 = 0;
   sRHS = 0;
}


void Elas3DH8::LMass(real_t coef)
{
   for (size_t i=0; i<8; ++i) {
      _hexa->setLocal(_xl[i]);
      real_t c = _rho*coef*_quad->getDet();
      eA2(3*i+1,3*i+1) += c;
      eA2(3*i+2,3*i+2) += c;
      eA2(3*i+3,3*i+3) += c;
   }
}


void Elas3DH8::Deviator(real_t coef)
{
   for (size_t k=0; k<8; ++k) {
      real_t c1 = coef*_wg[k]*_G;
      real_t c2 = 2*c1;
      for (size_t j=0; j<8; ++j) {
         real_t db11 = c2*_dSh[8*j+k].x, db22 = c2*_dSh[8*j+k].y;
         real_t db33 = c2*_dSh[8*j+k].z, db42 = c1*_dSh[8*j+k].z;
         real_t db43 = c1*_dSh[8*j+k].y, db53 = c1*_dSh[8*j+k].x;
         for (size_t i=0; i<8; ++i) {
            eA0(3*i+1,3*j+1) += _dSh[8*i+k].x*db11 + _dSh[8*i+k].z*db42 + _dSh[8*i+k].y*db43;
            eA0(3*i+1,3*j+2) += _dSh[8*i+k].y*db53;
            eA0(3*i+1,3*j+3) += _dSh[8*i+k].z*db53;
            eA0(3*i+2,3*j+1) += _dSh[8*i+k].x*db43;
            eA0(3*i+2,3*j+2) += _dSh[8*i+k].y*db22 + _dSh[8*i+k].z*db42 + _dSh[8*i+k].x*db53;
            eA0(3*i+2,3*j+3) += _dSh[8*i+k].z*db43;
            eA0(3*i+3,3*j+1) += _dSh[8*i+k].x*db42;
            eA0(3*i+3,3*j+2) += _dSh[8*i+k].y*db42;
            eA0(3*i+3,3*j+3) += _dSh[8*i+k].z*db33 + _dSh[8*i+k].y*db43 + _dSh[8*i+k].x*db53;
         }
      }
   }
   eMat += eA0;
}


void Elas3DH8::Dilatation(real_t coef)
{
   for (size_t k=0; k<8; k++) {
      real_t c = coef*_wg[k]*_lambda;
      for (size_t j=0; j<8; ++j) {
         real_t db11 = c*_dSh[8*j+k].x, db12 = c*_dSh[8*j+k].y, db13 = c*_dSh[8*j+k].z;
         for (size_t i=0; i<8; ++i) {
            eA0(3*i+1,3*j+1) += _dSh[8*i+k].x*db11;
            eA0(3*i+1,3*j+2) += _dSh[8*i+k].x*db12;
            eA0(3*i+1,3*j+3) += _dSh[8*i+k].x*db13;
            eA0(3*i+2,3*j+1) += _dSh[8*i+k].y*db11;
            eA0(3*i+2,3*j+2) += _dSh[8*i+k].y*db12;
            eA0(3*i+2,3*j+3) += _dSh[8*i+k].y*db13;
            eA0(3*i+3,3*j+1) += _dSh[8*i+k].z*db11;
            eA0(3*i+3,3*j+2) += _dSh[8*i+k].z*db12;
            eA0(3*i+3,3*j+3) += _dSh[8*i+k].z*db13;
         }
      }
   }
   eMat += eA0;
}


void Elas3DH8::BodyRHS(const Vect<real_t>& f)
{
   for (size_t k=0; k<8; ++k) {
      real_t fx = 0., fy = 0., fz = 0.;
      for (size_t j=0; j<8; ++j) {
         fx += _sh[8*j+k]*f((*_theElement)(j+1)->n(),1);
         fy += _sh[8*j+k]*f((*_theElement)(j+1)->n(),2);
         fz += _sh[8*j+k]*f((*_theElement)(j+1)->n(),3);
      }
      for (size_t i=0; i<8; ++i) {
         eRHS(3*i+1) += _wg[k]*fx*_sh[8*i+k];
         eRHS(3*i+2) += _wg[k]*fy*_sh[8*i+k];
         eRHS(3*i+3) += _wg[k]*fz*_sh[8*i+k];
      }
   }
}


void Elas3DH8::BodyRHS()
{
   for (size_t k=0; k<8; ++k) {
      real_t fx = 0., fy = 0., fz = 0.;
      for (size_t j=0; j<8; ++j) {
         fx += _sh[8*j+k]*_ebf(3*j  );
         fy += _sh[8*j+k]*_ebf(3*j+1);
         fz += _sh[8*j+k]*_ebf(3*j+2);
      }
      for (size_t i=0; i<8; ++i) {
         eRHS(3*i+1) += _wg[k]*fx*_sh[8*i+k];
         eRHS(3*i+2) += _wg[k]*fy*_sh[8*i+k];
         eRHS(3*i+3) += _wg[k]*fz*_sh[8*i+k];
      }
   }
}


void Elas3DH8::BoundaryRHS(const Vect<real_t>& f)
{
   for (size_t i=0; i<4; ++i) {
      for (size_t k=1; k<=3; k++) {
         if (_theSide->getCode(k+1) != CONTACT_BC)
            sRHS(3*i+k) += _el_geo.area*f(_theSide->n(),k);
      }
   }
}


void Elas3DH8::BoundaryRHS()
{
   for (size_t i=0; i<4; ++i) {
      for (size_t k=1; k<=3; ++k) {
         if (_theSide->getCode(k) != CONTACT_BC)
            sRHS(3*i+k) += _el_geo.area*_ssf[k-1];
      }
   }
}

} /* namespace OFELI */