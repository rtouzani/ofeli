/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

Beam3DL2::Beam3DL2(Mesh& ms)
         : Equation<real_t,2,12,1,6>(ms)
{
   _equation_name = "Beam equation";
   _finite_element = "2-D, 2-Node lines (P1)";
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
}


Beam3DL2::Beam3DL2(Mesh&         ms,
                   Vect<real_t>& u)
         : Equation<real_t,2,12,1,6>(ms,u)
{
   _equation_name = "Beam equation";
   _finite_element = "2-D, 2-Node lines (P1)";
   setMatrixType(SKYLINE|SYMMETRIC);
   setSolver(DIRECT_SOLVER);
}


void Beam3DL2::set(real_t A,
                   real_t I1,
                   real_t I2)
{
   _section.resize(_nb_el);
   _I1.resize(_nb_el);
   _I2.resize(_nb_el);
   for (size_t n=1; n<=_nb_el; ++n) {
      _section(n) = A;
      _I1(n) = I1;
      _I2(n) = I2;
   }
}


void Beam3DL2::set(const Vect<real_t>& A,
                   const Vect<real_t>& I1,
                   const Vect<real_t>& I2)
{
   _section.resize(_nb_el);
   _I1.resize(_nb_el);
   _I2.resize(_nb_el);
   for (size_t n=1; n<=_nb_el; ++n) {
      _section(n) = A(n);
      _I1(n) = I1(n);
      _I2(n) = I2(n);
   }
}


void Beam3DL2::set(const Element* el)
{
   _theElement = el;
   setMaterial();
   ElementNodeCoordinates();
   _x[0] = (*_theElement)(1)->getCoord();
   _x[1] = (*_theElement)(2)->getCoord();
   eA0 = 0, eA1 = 0, eA2 = 0;
   eRHS = 0;
   _mu = 0.5*_E/(1+_nu);
   _h = sqrt((_x[1].x - _x[0].x)*(_x[1].x - _x[0].x)
           + (_x[1].y - _x[0].y)*(_x[1].y - _x[0].y)
           + (_x[1].z - _x[0].z)*(_x[1].z - _x[0].z));
   _bending = _shear = _axial = _torsion = true;
   _I1e = _I1(_theElement->n());
   _I2e = _I2(_theElement->n());
   _ae = _section(_theElement->n());
   _reduced_integration = false;
   ElementVector(*_u);
}


void Beam3DL2::getDisp(Vect<real_t>& d)
{
   d.setSize(_nb_nodes,3);
   mesh_nodes(*_theMesh) {
      size_t n = node_label;
      d(n,1) = (*_u)(n,1);
      d(n,2) = (*_u)(n,2);
      d(n,3) = (*_u)(n,3);
   }
}


void Beam3DL2::build()
{
   AbsEqua<real_t>::_A->clear();
   mesh_elements(*_theMesh) {
      set(the_element);
      if (_analysis==TRANSIENT) {
         ElementVector(*AbsEqua<real_t>::_u);
         if (_terms&LUMPED_MASS)
            LMass();
         if (_terms&MASS)
            Mass();
      }
      Stiffness();
      if (_pf!=nullptr)
         *_b += *_pf;
      AbsEqua<real_t>::_A->Assembly(The_element,eMat.get());
      AbsEqua<real_t>::_b->Assembly(The_element,eRHS.get());
   }
}


void Beam3DL2::LMass(real_t coef)
{
   real_t c = 0.5*coef*_rho*_h;
   eA2( 1, 1) += 0.5*c;
   eA2( 2, 2) += 0.5*c;
   eA2( 3, 3) -= 0.5*c;
   eA2( 4, 4) -= 0.5*c*_I1e;
   eA2( 5, 5) -= 0.5*c*_I2e;
   eA2( 6, 6) += 0.5*c*(_I1e+_I2e);
   eA2( 7, 7) += 0.5*c;
   eA2( 8, 8) += 0.5*c;
   eA2( 9, 9) -= 0.5*c;
   eA2(10,10) -= 0.5*c*_I1e;
   eA2(11,11) -= 0.5*c*_I2e;
   eA2(12,12) += 0.5*c*(_I1e+_I2e);
}


void Beam3DL2::Stiffness(real_t coef)
{

// Bending
   if (_bending) {
      real_t c1 = coef*_E*_I1e/_h, c2 = coef*_E*_I2e/_h;
      eA0( 4, 4) += c1;
      eA0( 5, 5) += c2;
      eA0(10,10) += c1;
      eA0(11,11) += c2;
      eA0( 4,10) -= c1;
      eA0( 5,11) -= c2;
   }

// Shear
   if (_shear) {
      real_t c1 = coef*_mu*_ae/_h, c2 = coef*_mu*_ae*_h*0.5, c3 = coef*_mu*_ae*_h*OFELI_THIRD;
      eA0( 1, 1) += c1;
      eA0( 2, 2) += c1;
      eA0( 7, 7) += c1;
      eA0( 8, 8) += c1;
      eA0( 1, 7) -= c1;
      eA0( 2, 8) -= c1;

      if (_reduced_integration) {
         eA0( 4, 4) += c2;
         eA0( 5, 5) -= c2;
         eA0(10,10) += c2;
         eA0(11,11) -= c2;
      }
      else {
         eA0( 4, 4) += c3;
         eA0( 5, 5) -= c3;
         eA0(10,10) += c3;
         eA0(11,11) -= c3;
         eA0( 4,10) += c3;
         eA0( 5,11) += c3;
      }

      eA0( 2, 4) -= c2;
      eA0( 1, 5) += c2;
      eA0( 2,10) -= c2;
      eA0( 1,11) += c2;
      eA0( 5, 7) -= c2;
      eA0( 4, 8) += c2;
      eA0( 7,11) -= c2;
      eA0( 6,10) += c2;
   }

// Axial
   if (_axial) {
      real_t c1 = coef*_E*_ae/_h;
      eA0( 3, 3) += c1;
      eA0( 3, 9) -= c1;
      eA0( 9, 9) += c1;
   }

// Torsional
   if (_torsion) {
      real_t c1 = coef*_mu*(_I1e+_I2e)/_h;
      eA0( 6, 6) += c1;
      eA0( 6,12) -= c1;
      eA0(12,12) += c1;
   }

// Symmetrize matrix
   eA0.Symmetrize();
   eMat = eA0;
}


void Beam3DL2::Load(const Vect<real_t> &f)
{
  for (size_t i=1; i<=6; ++i)
     eRHS(i) += f((*_theElement)(1)->n(),i);
}


void Beam3DL2::AxialForce(Vect<real_t>& f)
{
   f.setSize(_nb_el);
   mesh_elements(*_theMesh) {
     set(the_element);
     f(element_label) = _E*_ae*(_eu(9)-_eu(3))/_h;
   }
}


void Beam3DL2::ShearForce(Vect<real_t>& sh)
{
   sh.setSize(_nb_el,2);
   mesh_elements(*_theMesh) {
     set(the_element);
     sh(element_label,1) = _mu*_ae*((_eu(7)-_eu(1))/_h - 0.5*(_eu(11)+_eu(5)));
     sh(element_label,2) = _mu*_ae*((_eu(8)-_eu(2))/_h - 0.5*(_eu(10)+_eu(4)));
   }
}


void Beam3DL2::BendingMoment(Vect<real_t>& m)
{
   m.setSize(_nb_el,2);
   mesh_elements(*_theMesh) {
      set(the_element);
      real_t c = -_E*(_eu(10)-_eu(4))/_h;
      m(element_label,1) = c*_I2e;
      m(element_label,2) = c*_I1e;
   }
}


void Beam3DL2::TwistingMoment(Vect<real_t>& m)
{
   m.setSize(_nb_el);
   mesh_elements(*_theMesh) {
      set(the_element);
      m(element_label) = _mu*(_I1e+_I2e)*(_eu(12)-_eu(6))/_h;
   }
}


void Beam3DL2::buildEigen(SkSMatrix<real_t>& K,
                          Vect<real_t>&      M)
{
   real_t c, c1, c2, c3;
   mesh_elements(*_theMesh) {
      set(the_element);
      if (_bending) {
         c1 = _E*_I1e/_h;
         c2 = _E*_I2e/_h;
         eA0( 4, 4) += c1; eA0( 5, 5) += c2;
         eA0(10,10) += c1; eA0(11,11) += c2;
         eA0( 4,10) -= c1; eA0( 5,11) -= c2;
      }
      if (_shear) {
         c1 = _mu*_ae/_h;
         c2 = _mu*_ae*_h*0.5;
         c3 = _mu*_ae*_h*OFELI_THIRD;
         eA0( 1, 1) += c1; eA0( 2, 2) += c1;
         eA0( 7, 7) += c1; eA0( 8, 8) += c1;
         eA0( 1, 7) -= c1; eA0( 2, 8) -= c1;
         if (_reduced_integration) {
            eA0( 4, 4) += c2; eA0( 5, 5) -= c2;
            eA0(10,10) += c2; eA0(11,11) -= c2;
         }
         else {
            eA0( 4, 4) += c3; eA0( 5, 5) -= c3;
            eA0(10,10) += c3; eA0(11,11) -= c3;
            eA0( 4,10) += c3; eA0( 5,11) += c3;
         }
         eA0( 2, 4) -= c2; eA0( 1, 5) += c2;
         eA0( 2,10) -= c2; eA0( 1,11) += c2;
         eA0( 5, 7) -= c2; eA0( 4, 8) += c2;
         eA0( 7,11) -= c2; eA0( 6,10) += c2;
      }
      if (_axial) {
         c1 = _E*_ae/_h;
         eA0( 3, 3) += c1; eA0( 3, 9) -= c1; eA0( 9, 9) += c1;
      }
      if (_torsion) {
         c1 = _mu*(_I1e+_I2e)/_h;
         eA0( 6, 6) += c1; eA0( 6,12) -= c1; eA0(12,12) += c1;
      }
      eMat.Symmetrize();
      c = 0.5*_rho*_h;
      eRHS( 1) += 0.5*c;      eRHS( 2) += 0.5*c;
      eRHS( 3) -= 0.5*c;      eRHS( 4) -= 0.5*c*_I1e;
      eRHS( 5) -= 0.5*c*_I2e; eRHS( 6) += 0.5*c*(_I1e+_I2e);
      eRHS( 7) += 0.5*c;      eRHS( 8) += 0.5*c;
      eRHS( 9) -= 0.5*c;      eRHS(10) -= 0.5*c*_I1e;
      eRHS(11) -= 0.5*c*_I2e; eRHS(12) += 0.5*c*(_I1e+_I2e);
      ElementAssembly(K);
      ElementAssembly(M);
   }
}

} /* namespace OFELI */
