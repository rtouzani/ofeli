/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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
                         Implementation of class Beam3DL2
  ==============================================================================*/


#include "equations/solid/Beam3DL2.h"

namespace OFELI {

Beam3DL2::Beam3DL2(Element* el,
                   real_t   A,
                   real_t   I1,
                   real_t   I2)
{
   set(el);
   _I1 = I1;
   _I2 = I2;
   _section = A;
}


Beam3DL2::Beam3DL2(      Element*      el,
                         real_t        A,
                         real_t        I1,
                         real_t        I2,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   _time = time;
   set(el);
   _section = A;
   _I1 = I1;
   _I2 = I2;
   ElementVector(u);
}


Beam3DL2::Beam3DL2(      Mesh&         ms,
                   const Vect<real_t>& u,
                         Vect<real_t>& d)
{
   mesh_nodes(ms) {
      size_t n = node_label;
      d(n,1) = u(n,1);
      d(n,2) = u(n,2);
      d(n,3) = u(n,3);
   }
}


void Beam3DL2::set(const Element* el)
{
   _nb_dof = 6;
   Init(el);
   setMaterial();
   ElementNodeCoordinates();
   _x[0] = (*_theElement)(1)->getCoord();
   _x[1] = (*_theElement)(2)->getCoord();
   eMat = 0;
   eRHS = 0;
   _mu = 0.5*_E/(1+_nu);
   _h = sqrt((_x[1].x - _x[0].x)*(_x[1].x - _x[0].x)
           + (_x[1].y - _x[0].y)*(_x[1].y - _x[0].y)
           + (_x[1].z - _x[0].z)*(_x[1].z - _x[0].z));
   _bending = _shear = _axial = _torsion = true;
   _reduced_integration = false;
}


void Beam3DL2::LMassToLHS(real_t coef)
{
   real_t c = 0.5*coef*_rho*_h;
   eMat( 1, 1) += 0.5*c;
   eMat( 2, 2) += 0.5*c;
   eMat( 3, 3) -= 0.5*c;
   eMat( 4, 4) -= 0.5*c*_I1;
   eMat( 5, 5) -= 0.5*c*_I2;
   eMat( 6, 6) += 0.5*c*(_I1+_I2);
   eMat( 7, 7) += 0.5*c;
   eMat( 8, 8) += 0.5*c;
   eMat( 9, 9) -= 0.5*c;
   eMat(10,10) -= 0.5*c*_I1;
   eMat(11,11) -= 0.5*c*_I2;
   eMat(12,12) += 0.5*c*(_I1+_I2);
}


void Beam3DL2::LMassToRHS(real_t coef)
{
   real_t c = 0.5*coef*_rho*_h;
   eRHS( 1) += c * ePrev( 1);
   eRHS( 2) += c * ePrev( 2);
   eRHS( 3) += c * ePrev( 3);
   eRHS( 4) += c * ePrev( 4);
   eRHS( 5) += c * ePrev( 5);
   eRHS( 6) += c * ePrev( 6);
   eRHS( 7) += c * ePrev( 7);
   eRHS( 8) += c * ePrev( 8);
   eRHS( 9) += c * ePrev( 9);
   eRHS(10) += c * ePrev(10);
   eRHS(11) += c * ePrev(11);
   eRHS(12) += c * ePrev(12);
}


void Beam3DL2::Stiffness(real_t coef)
{
   real_t c1, c2, c3;

// Bending
   if (_bending) {
      c1 = coef*_E*_I1/_h;
      c2 = coef*_E*_I2/_h;
      eMat( 4, 4) += c1;
      eMat( 5, 5) += c2;
      eMat(10,10) += c1;
      eMat(11,11) += c2;
      eMat( 4,10) -= c1;
      eMat( 5,11) -= c2;
   }

// Shear
   if (_shear) {
      c1 = coef*_mu*_section/_h;
      c2 = coef*_mu*_section*_h*0.5;
      c3 = coef*_mu*_section*_h*OFELI_THIRD;
      eMat( 1, 1) += c1;
      eMat( 2, 2) += c1;
      eMat( 7, 7) += c1;
      eMat( 8, 8) += c1;
      eMat( 1, 7) -= c1;
      eMat( 2, 8) -= c1;

      if (_reduced_integration) {
         eMat( 4, 4) += c2;
         eMat( 5, 5) -= c2;
         eMat(10,10) += c2;
         eMat(11,11) -= c2;
      }
      else {
         eMat( 4, 4) += c3;
         eMat( 5, 5) -= c3;
         eMat(10,10) += c3;
         eMat(11,11) -= c3;
         eMat( 4,10) += c3;
         eMat( 5,11) += c3;
      }

      eMat( 2, 4) -= c2;
      eMat( 1, 5) += c2;
      eMat( 2,10) -= c2;
      eMat( 1,11) += c2;
      eMat( 5, 7) -= c2;
      eMat( 4, 8) += c2;
      eMat( 7,11) -= c2;
      eMat( 6,10) += c2;
   }

// Axial
   if (_axial) {
      c1 = coef*_E*_section/_h;
      eMat( 3, 3) += c1;
      eMat( 3, 9) -= c1;
      eMat( 9, 9) += c1;
   }

// Torsional
   if (_torsion) {
      c1 = coef*_mu*(_I1+_I2)/_h;
      eMat( 6, 6) += c1;
      eMat( 6,12) -= c1;
      eMat(12,12) += c1;
   }

// Symmetrize matrix
   eMat.Symmetrize();
}


void Beam3DL2::Load(const Vect<real_t> &f)
{
   eRHS(1) += f(_theElement->getNodeLabel(1),1);
   eRHS(2) += f(_theElement->getNodeLabel(1),2);
   eRHS(3) += f(_theElement->getNodeLabel(1),3);
   eRHS(4) += f(_theElement->getNodeLabel(1),4);
   eRHS(5) += f(_theElement->getNodeLabel(1),5);
   eRHS(6) += f(_theElement->getNodeLabel(1),6);
}


real_t Beam3DL2::AxialForce() const
{
   return _E*_section*(ePrev(9)-ePrev(3))/_h;
}


Point<real_t> Beam3DL2::ShearForce() const
{
   real_t a = _mu*_section*((ePrev(7)-ePrev(1))/_h - 0.5*(ePrev(11)+ePrev(5)));
   real_t b = _mu*_section*((ePrev(8)-ePrev(2))/_h - 0.5*(ePrev(10)+ePrev(4)));
   return Point<real_t>(a,b);
}


Point<real_t> Beam3DL2::BendingMoment() const
{
   Point<real_t> m;
   real_t c = -_E*(ePrev(10)-ePrev(4))/_h;
   m.x = c*_I2;
   m.y = c*_I1;
   return m;
}


real_t Beam3DL2::TwistingMoment() const
{
   return _mu*(_I1+_I2)*(ePrev(12)-ePrev(6))/_h;
}


void Beam3DL2::buildEigen(SkSMatrix<real_t>& K,
                          Vect<real_t>&      M)
{
   real_t c, c1, c2, c3;
   MESH_EL {
      set(theElement);
      if (_bending) {
         c1 = _E*_I1/_h;
         c2 = _E*_I2/_h;
         eMat( 4, 4) += c1; eMat( 5, 5) += c2;
         eMat(10,10) += c1; eMat(11,11) += c2;
         eMat( 4,10) -= c1; eMat( 5,11) -= c2;
      }
      if (_shear) {
         c1 = _mu*_section/_h;
         c2 = _mu*_section*_h*0.5;
         c3 = _mu*_section*_h*OFELI_THIRD;
         eMat( 1, 1) += c1; eMat( 2, 2) += c1;
         eMat( 7, 7) += c1; eMat( 8, 8) += c1;
         eMat( 1, 7) -= c1; eMat( 2, 8) -= c1;
         if (_reduced_integration) {
            eMat( 4, 4) += c2; eMat( 5, 5) -= c2;
            eMat(10,10) += c2; eMat(11,11) -= c2;
         }
         else {
            eMat( 4, 4) += c3; eMat( 5, 5) -= c3;
            eMat(10,10) += c3; eMat(11,11) -= c3;
            eMat( 4,10) += c3; eMat( 5,11) += c3;
         }
         eMat( 2, 4) -= c2; eMat( 1, 5) += c2;
         eMat( 2,10) -= c2; eMat( 1,11) += c2;
         eMat( 5, 7) -= c2; eMat( 4, 8) += c2;
         eMat( 7,11) -= c2; eMat( 6,10) += c2;
      }
      if (_axial) {
         c1 = _E*_section/_h;
         eMat( 3, 3) += c1; eMat( 3, 9) -= c1; eMat( 9, 9) += c1;
      }
      if (_torsion) {
         c1 = _mu*(_I1+_I2)/_h;
         eMat( 6, 6) += c1; eMat( 6,12) -= c1; eMat(12,12) += c1;
      }
      eMat.Symmetrize();
      c = 0.5*_rho*_h;
      eRHS( 1) += 0.5*c;     eRHS( 2) += 0.5*c;
      eRHS( 3) -= 0.5*c;     eRHS( 4) -= 0.5*c*_I1;
      eRHS( 5) -= 0.5*c*_I2; eRHS( 6) += 0.5*c*(_I1+_I2);
      eRHS( 7) += 0.5*c;     eRHS( 8) += 0.5*c;
      eRHS( 9) -= 0.5*c;     eRHS(10) -= 0.5*c*_I1;
      eRHS(11) -= 0.5*c*_I2; eRHS(12) += 0.5*c*(_I1+_I2);
      ElementAssembly(K);
      ElementAssembly(M);
   }
}

} /* namespace OFELI */
