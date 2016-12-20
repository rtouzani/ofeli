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

                         Implementation of class ELAS3DH8
  for 3-D Linear Elasticity Equations using 8-node hexahedral Q1 finite element

  ==============================================================================*/


#include "equations/solid/Elas3DH8.h"

namespace OFELI {

Elas3DH8::Elas3DH8(const Element* el)
{
   set(el);
}


Elas3DH8::Elas3DH8(const Side* sd)
{
   set(sd);
}


Elas3DH8::Elas3DH8(const Element*      el,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   set(el);
   _time = time;
   ElementVector(u);
}


Elas3DH8::Elas3DH8(const Side*         sd,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   set(sd);
   _time = time;
   SideVector(u);
}


Elas3DH8::~Elas3DH8()
{
   if (_hexa) { 
      delete _hexa; 
      _hexa = NULL;
   }
   if (_quad) { 
      delete _quad;
      _quad = NULL;
   }
}


void Elas3DH8::set(const Element* el)
{
   _nb_dof = 3;
   Init(el);
   setMaterial();
   _hexa = new Hexa8(el);
   _quad = NULL;
   ElementNodeCoordinates();
   _xl[0].x = _xl[3].x = _xl[4].x = _xl[7].x = -1.;
   _xl[1].x = _xl[2].x = _xl[5].x = _xl[6].x =  1.;
   _xl[0].y = _xl[1].y = _xl[4].y = _xl[5].y = -1.;
   _xl[2].y = _xl[3].y = _xl[6].y = _xl[7].y =  1.;
   _xl[0].z = _xl[1].z = _xl[2].z = _xl[3].z = -1.;
   _xl[4].z = _xl[5].z = _xl[6].z = _xl[7].z =  1.;
   _lambda = _nu*_E/((1+_nu)*(1-2*_nu));
   _G = 0.5*_E/(1+_nu);
   _hexa->atGauss2(_dsh,_w);
   Gauss g(2);
   _xg[0] = g.x(1); _xg[1] = g.x(2);
   _wg[0] = g.w(1); _wg[1] = g.w(2);
   eMat = 0;
   eRHS = 0;
}


void Elas3DH8::set(const Side* sd)
{
   _nb_dof = 3;
   Init(sd);
   _quad = new Quad4(sd);
   _hexa = NULL;
   SideNodeCoordinates();
   _area = _theSide->getMeasure();
   Gauss g(2);
   _xg[0] = g.x(1); _xg[1] = g.x(2);
   _wg[0] = g.w(1); _wg[1] = g.w(2);
   sMat = 0;
   sRHS = 0;
}


void Elas3DH8::LMassToLHS(real_t coef)
{
   for (size_t i=1; i<=8; i++) {
      _hexa->setLocal(_xl[i-1]);
      real_t c = _rho*coef*_quad->getDet();
      eMat(3*i-2,3*i-2) += c;
      eMat(3*i-1,3*i-1) += c;
      eMat(3*i  ,3*i  ) += c;
   }
}


void Elas3DH8::LMassToRHS(real_t coef)
{
   for (size_t i=1; i<=8; i++) {
      _hexa->setLocal(_xl[i-1]);
      real_t c = _rho*coef*_hexa->getDet();
      eRHS(3*i-2) += c*ePrev(3*i-2);
      eRHS(3*i-1) += c*ePrev(3*i-1);
      eRHS(3*i  ) += c*ePrev(3*i  );
   }
}


void Elas3DH8::Deviator(real_t coef)
{
   real_t db11, db22, db33, db42, db43, db53;
   for (size_t k=0; k<8; k++) {
      real_t c1 = coef*_w[k]*_G;
      real_t c2 = 2*c1;
      for (size_t j=1; j<=8; j++) {
         db11 = c2*_dsh(j,k).x;
         db22 = c2*_dsh(j,k).y;
         db33 = c2*_dsh(j,k).z;
         db42 = c1*_dsh(j,k).z;
         db43 = c1*_dsh(j,k).y;
         db53 = c1*_dsh(j,k).x;
         for (size_t i=1; i<=8; i++) {
            eMat(3*i-2,3*j-2) += _dsh(i,k).x*db11 + _dsh(i,k).z*db42 + _dsh(i,k).y*db43;
            eMat(3*i-2,3*j-1) += _dsh(i,k).y*db53;
            eMat(3*i-2,3*j  ) += _dsh(i,k).z*db53;
            eMat(3*i-1,3*j-2) += _dsh(i,k).x*db43;
            eMat(3*i-1,3*j-1) += _dsh(i,k).y*db22 + _dsh(i,k).z*db42 + _dsh(i,k).x*db53;
            eMat(3*i-1,3*j  ) += _dsh(i,k).z*db43;
            eMat(3*i  ,3*j-2) += _dsh(i,k).x*db42;
            eMat(3*i  ,3*j-1) += _dsh(i,k).y*db42;
            eMat(3*i  ,3*j  ) += _dsh(i,k).z*db33 + _dsh(i,k).y*db43 + _dsh(i,k).x*db53;
         }
      }
   }
}


void Elas3DH8::DeviatorToRHS(real_t coef)
{
   real_t db11, db22, db33, db42, db43, db53;
   for (size_t k=0; k<8; k++) {
      real_t c1 = coef*_w[k]*_G;
      real_t c2 = 2*c1;
      for (size_t j=1; j<=8; j++) {
         db11 = c2*_dsh(j,k).x;
         db22 = c2*_dsh(j,k).y;
         db33 = c2*_dsh(j,k).z;
         db42 = c1*_dsh(j,k).z;
         db43 = c1*_dsh(j,k).y;
         db53 = c1*_dsh(j,k).x;
         for (size_t i=1; i<=8; i++) {
            eRHS(3*i-2) -= (_dsh(i,k).x*db11+_dsh(i,k).z*db42+_dsh(i,k).y*db43) * ePrev(3*j-2)
                         + (_dsh(i,k).y*db53) * ePrev(3*j-1)
                         + (_dsh(i,k).z*db53) * ePrev(3*j  );
            eRHS(3*i-1) -= (_dsh(i,k).x*db43) * ePrev(3*j-2)
                         + (_dsh(i,k).y*db22+_dsh(i,k).z*db42+_dsh(i,k).x*db53) * ePrev(3*j-1)
                         + (_dsh(i,k).z*db43) * ePrev(3*j  );
            eRHS(3*i  ) -= (_dsh(i,k).x*db42) * ePrev(3*j-2)
                         + (_dsh(i,k).y*db42) * ePrev(3*j-1)
                         + (_dsh(i,k).z*db33+_dsh(i,k).y*db43+_dsh(i,k).x*db53) * ePrev(3*j  );
         }
      }
   }
}


void Elas3DH8::DilatationToRHS(real_t coef)
{
   real_t db11, db12, db13;
   for (size_t k=0; k<8; k++) {
      real_t c = coef*_w[k]*_lambda;
      for (size_t j=1; j<=8; j++) {
         db11 = c*_dsh(j,k).x;
         db12 = c*_dsh(j,k).y;
         db13 = c*_dsh(j,k).z;
         for (size_t i=1; i<=8; i++) {
            eRHS(3*i-2) -= _dsh(i,k).x*db11 * ePrev(3*j-2)
                         + _dsh(i,k).x*db12 * ePrev(3*j-1)
                         + _dsh(i,k).x*db13 * ePrev(3*j  );
            eRHS(3*i-1) -= _dsh(i,k).y*db11 * ePrev(3*j-2)
                         + _dsh(i,k).y*db12 * ePrev(3*j-1)
                         + _dsh(i,k).y*db13 * ePrev(3*j  );
            eRHS(3*i  ) -= _dsh(i,k).z*db11 * ePrev(3*j-2)
                         + _dsh(i,k).z*db12 * ePrev(3*j-1)
                         + _dsh(i,k).z*db13 * ePrev(3*j  );
         }
      }
   }
}


void Elas3DH8::Dilatation(real_t coef)
{
   real_t db11, db12, db13;
   for (size_t k=0; k<8; k++) {
      real_t c = coef*_w[k]*_lambda;
      for (size_t j=1; j<=8; j++) {
         db11 = c*_dsh(j,k).x;
         db12 = c*_dsh(j,k).y;
         db13 = c*_dsh(j,k).z;
         for (size_t i=1; i<=8; i++) {
            eMat(3*i-2,3*j-2) += _dsh(i,k).x*db11;
            eMat(3*i-2,3*j-1) += _dsh(i,k).x*db12;
            eMat(3*i-2,3*j  ) += _dsh(i,k).x*db13;
            eMat(3*i-1,3*j-2) += _dsh(i,k).y*db11;
            eMat(3*i-1,3*j-1) += _dsh(i,k).y*db12;
            eMat(3*i-1,3*j  ) += _dsh(i,k).y*db13;
            eMat(3*i  ,3*j-2) += _dsh(i,k).z*db11;
            eMat(3*i  ,3*j-1) += _dsh(i,k).z*db12;
            eMat(3*i  ,3*j  ) += _dsh(i,k).z*db13;
         }
      }
   }
}


void Elas3DH8::BodyRHS(UserData<real_t>& ud)
{
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++)
         for (size_t m=0; m<2; m++) {
            Point<real_t> g(_xg[k],_xg[l],_xg[m]);
            _hexa->setLocal(g);
            Point<real_t> x = _hexa->getLocalPoint();
            real_t fx = ud.BodyForce(x, _time, 1),
                   fy = ud.BodyForce(x, _time, 2),
                   fz = ud.BodyForce(x, _time, 3);
            real_t c = _wg[k]*_wg[l]*_wg[m]*_hexa->getDet();
            for (size_t i=1; i<=8; i++) {
               eRHS(3*i-2) += c*fx*_hexa->Sh(i);
               eRHS(3*i-1) += c*fy*_hexa->Sh(i);
               eRHS(3*i  ) += c*fz*_hexa->Sh(i);
         }
      }
}


void Elas3DH8::BodyRHS(const Vect<real_t>& bf,
                             int           opt)
{
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++)
         for (size_t m=0; m<2; m++) {
            Point<real_t> g(_xg[k],_xg[l],_xg[m]);
            _hexa->setLocal(g);
            real_t c = _wg[k]*_wg[l]*_wg[m]*_hexa->getDet();
            real_t fx = 0., fy = 0., fz = 0.;
            if (opt==LOCAL_ARRAY) {
               for (size_t j=1; j<=8; j++) {
                  fx += _hexa->Sh(j)*bf(3*j-2);
                  fy += _hexa->Sh(j)*bf(3*j-1);
                  fz += _hexa->Sh(j)*bf(3*j  );
               }
            }
            else {
               for (size_t j=1; j<=8; j++) {
                  fx += _hexa->Sh(j)*bf(3*_theElement->getNodeLabel(j)-2);
                  fy += _hexa->Sh(j)*bf(3*_theElement->getNodeLabel(j)-1);
                  fz += _hexa->Sh(j)*bf(3*_theElement->getNodeLabel(j)  );
               }
            }
            for (size_t i=1; i<=8; i++) {
               eRHS(3*i-2) += c*fx*_hexa->Sh(i);
               eRHS(3*i-1) += c*fy*_hexa->Sh(i);
               eRHS(3*i  ) += c*fz*_hexa->Sh(i);
            }
         }
}


void Elas3DH8::BoundaryRHS(const Vect<real_t>& f)
{
   for (size_t i=1; i<=4; i++) {
      for (size_t k=1; k<=3; k++) {
         if (_theSide->getCode(k) != CONTACT) {
            real_t c = _area*f(_theSide->n(),k);
            sRHS(3*(i-1)+k) += c;
         }
      }
   }
}

} /* namespace OFELI */
