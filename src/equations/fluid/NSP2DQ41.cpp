/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

                        Implementation of class NSP2DQ41
             for 2-D Navier-Stokes equations using penalty formulation and
                               Q1/P0 finite element

  ==============================================================================*/


#include "equations/fluid/NSP2DQ41.h"

namespace OFELI {

NSP2DQ41::NSP2DQ41(const Element* el)
{
   _nb_dof = 2;
   Init(el);
   setMaterial();
   _quad = new Quad4(_theElement);
   _ln = NULL;
   _gauss[0] = -0.577350269189626;
   _gauss[1] = -_gauss[0];
   _cgx = _cgy = 0.;
   _xl[0] = _yl[0] = _yl[1] = _xl[3] = -1.;
   _xl[1] = _xl[2] = _yl[2] = _yl[3] =  1.;
   for (size_t i=0; i<4; i++)
      _x[i] = _theElement->getPtrNode(i+1)->getCoord();
   eMat = 0;
   eRHS = 0;
}


NSP2DQ41::NSP2DQ41(const Element*      el,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   _nb_dof = 2;
   Init(el);
   _time = time;
   setMaterial();
   _quad = new Quad4(el);
   _ln = NULL;
   _gauss[0] = -0.577350269189626;
   _gauss[1] = -_gauss[0];
   _cgx = _cgy = 0.;
   _xl[0] = _yl[0] = _yl[1] = _xl[3] = -1.;
   _xl[1] = _xl[2] = _yl[2] = _yl[3] =  1.;
   ElementVector(u);
   for (size_t i=0; i<4; i++)
      _x[i] = (*_theElement)(i+1)->getCoord();
   eMat = 0;
   eRHS = 0;
}


NSP2DQ41::NSP2DQ41(const Side* sd)
{
   _nb_dof = 2;
   Init(sd);
   _ln = new Line2(_theSide);
   _quad = NULL;
   _gauss[0] = -0.577350269189626;
   _gauss[1] = -_gauss[0];
   _cgx = _cgy = 0.;
   _xl[0] = _yl[0] = _yl[1] = _xl[3] = -1.;
   _xl[1] = _xl[2] = _yl[2] = _yl[3] =  1.;
   SideNodeCoordinates();
   sMat = 0;
   sRHS = 0;
}


NSP2DQ41::NSP2DQ41(const Side*         sd,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   _nb_dof = 2;
   Init(sd);
   _time = time;
   _ln = new Line2(_theSide);
   _quad = NULL;
   _gauss[0] = -0.577350269189626;
   _gauss[1] = -_gauss[0];
   _cgx = _cgy = 0.;
   SideVector(u);
   SideNodeCoordinates();
   sMat = 0;
   sRHS = 0;
}


NSP2DQ41::~NSP2DQ41()
{
   if (_quad) { delete _quad; _quad = NULL; }
   if (_ln) { delete _ln; _ln = NULL; }
}


void NSP2DQ41::LMass(real_t coef)
{
   for (size_t i=1; i<=4; i++) {
      _quad->setLocal(Point<real_t> (_xl[i-1],_yl[i-1]));
      real_t c = _dens*coef*_quad->getDet();
      eMat(2*i-1,2*i-1) += c;
      eMat(2*i  ,2*i  ) += c;
      eRHS(2*i-1) += c*ePrev(2*i-1);
      eRHS(2*i  ) += c*ePrev(2*i  );
   }
}


void NSP2DQ41::Mass(real_t coef)
{
   coef = 1;
   std::cerr << "Sorry, consistent mass matrix is not implemented !\n";
}


void NSP2DQ41::Viscous(real_t coef)
{
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++) {
         Point<real_t>  g(_gauss[k],_gauss[l]);
         _quad->setLocal(g);
         real_t c = coef*_visc*_quad->getDet();
         for (size_t i=1; i<=4; i++) {
            Point<real_t> a = c*_quad->DSh(i);
            for (size_t j=1; j<=4; j++) {
               eMat(2*i-1,2*j-1) += 2*a.x*_quad->DSh(j).x + a.y*_quad->DSh(j).y;
               eMat(2*i-1,2*j  ) += a.y*_quad->DSh(j).x;
               eMat(2*i  ,2*j-1) += a.x*_quad->DSh(j).y;
               eMat(2*i  ,2*j  ) += 2*a.y*_quad->DSh(j).y + a.x*_quad->DSh(j).x;
            }
         }
      }
}


void NSP2DQ41::RHS_Viscous(real_t coef)
{
   coef = 1;
   std::cerr << "Sorry, not yet implemented !\n";
}


void NSP2DQ41::Penal(real_t coef)
{
   _quad->setLocal(Point<real_t> (0.));
   Point<real_t>  c = 4.*_quad->getDet()*coef;
   for (size_t i=1; i<=4; i++) {
      Point<real_t>  a = c*_quad->DSh(i);;
      for (size_t j=1; j<=4; j++) {
         eMat(2*i-1,2*j-1) += a.x*_quad->DSh(j).x;
         eMat(2*i-1,2*j  ) += a.x*_quad->DSh(j).y;
         eMat(2*i  ,2*j-1) += a.y*_quad->DSh(j).x;
         eMat(2*i  ,2*j  ) += a.y*_quad->DSh(j).y;
      }
   }
}


void NSP2DQ41::RHS_Convection(real_t coef)
{
   static real_t ug[2];
   Point<real_t>  dug[2];
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++) {
         _quad->setLocal(Point<real_t> (_gauss[k],_gauss[l]));
         for (size_t j=0; j<2; j++) {
            ug[j] = dug[j].x = dug[j].y = 0.;
            for (size_t i=1; i<=4; i++) {
               size_t ii = 2*i+j-1;
               ug[j]    += _quad->Sh(i)*ePrev(ii);
               dug[j].x += _quad->DSh(i).x*ePrev(ii);
               dug[j].y += _quad->DSh(i).y*ePrev(ii);
            }
         }

         real_t t1 = ug[0]*dug[0].x + ug[1]*dug[0].y;
         real_t t2 = ug[0]*dug[1].x + ug[1]*dug[1].y;
         for (size_t k=1; k<=4; k++) {
            real_t c = coef * _dens * _quad->getDet() * _quad->Sh(k);
            eRHS(2*k-1) -= c*t1;
            eRHS(2*k  ) -= c*t2;
         }
      }
}


void NSP2DQ41::LHS1_Convection(real_t coef)
{
   static real_t ug[2];
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++) {
         _quad->setLocal(Point<real_t> (_gauss[k],_gauss[l]));
         for (size_t j=0; j<2; j++) {
            ug[j] = 0.;
            for (size_t i=1; i<=4; i++)
               ug[j] += _quad->Sh(i)*ePrev(2*i+j-1);
         }
         for (size_t ii=1; ii<=4; ii++)
            for (size_t jj=1; jj<=4; jj++) {
               real_t c = coef*_dens*_quad->getDet()*(ug[0]*_quad->DSh(jj).x + ug[1]*_quad->DSh(jj).y);
               eMat(2*ii-1,2*jj-1) += c*_quad->Sh(ii);
               eMat(2*ii  ,2*jj  ) += c*_quad->Sh(ii);
            }
      }
}


void NSP2DQ41::LHS2_Convection(real_t coef)
{
   Point<real_t>  dug[2];
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++) {
         _quad->setLocal(Point<real_t> (_gauss[k],_gauss[l]));
         for (size_t j=0; j<2; j++) {
            dug[j].x = dug[j].y = 0.;
            for (size_t i=1; i<=4; i++) {
               dug[j].x += _quad->DSh(i).x*ePrev(2*i+j-1)*coef;
               dug[j].y += _quad->DSh(i).y*ePrev(2*i+j-1)*coef;
            }
         }
         for (size_t ii=1; ii<=4; ii++)
            for (size_t jj=1; jj<=4; jj++) {
               real_t c = _dens * _quad->getDet() * _quad->Sh(ii) * _quad->Sh(jj);
               eMat(2*ii-1,2*jj-1) += c*dug[0].x;
               eMat(2*ii-1,2*jj  ) += c*dug[0].y;
               eMat(2*ii  ,2*jj-1) += c*dug[1].x;
               eMat(2*ii  ,2*jj  ) += c*dug[1].y;
            }
      }
}


void NSP2DQ41::BodyRHS(UserData<real_t>& ud)
{
   for (size_t k=0; k<2; k++)
      for (size_t l=0; l<2; l++) {
         _quad->setLocal(Point<real_t> (_gauss[k],_gauss[l]));
         Point<real_t> x(_quad->getLocalPoint());
         real_t fx = ud.BodyForce(x, _time, 1);
         real_t fy = ud.BodyForce(x, _time, 2);
         real_t c = _quad->getDet();
         for (size_t i=1; i<=4; i++) {
            eRHS(2*i-1) += c*fx*_quad->Sh(i);
            eRHS(2*i  ) += c*fy*_quad->Sh(i);
         }
      }
}


void NSP2DQ41::BoundaryRHS(UserData<real_t>& ud)
{
   for (size_t k=0; k<2; k++) {
      real_t g = _gauss[k];
      real_t fx = ud.SurfaceForce(_ln->getLocalPoint(g), _theSide->getCode(1), _time, 1);
      real_t fy = ud.SurfaceForce(_ln->getLocalPoint(g), _theSide->getCode(2), _time, 2);
      real_t c = 0.5*_ln->getLength();
      for (size_t i=1; i<=2; i++) {
         sRHS(2*i-1) += c*fx*_ln->Sh(i,g);
         sRHS(2*i  ) += c*fy*_ln->Sh(i,g);
      }
   }
}


void NSP2DQ41::Periodic(real_t coef)
{
   for (size_t i=1; i<=2; i++) {
      real_t c = 0.5*_ln->getLength()*coef;
      if (_theSide->getCode(1) == PERIODIC_A)
         sMat(2*i-1,2*i-1) += c;
      else if (_theSide->getCode(1) == PERIODIC_B)
         sMat(2*i-1,2*i-1) -= c;
      if (_theSide->getCode(2) == PERIODIC_A)
         sMat(2*i  ,2*i  ) += c;
      else if (_theSide->getCode(2) == PERIODIC_B)
         sMat(2*i  ,2*i  ) -= c;
   }
}


real_t NSP2DQ41::Pressure(real_t coef)
{
   real_t p = 0.;
   for (size_t k=0; k<2; k++) {
      for (size_t l=0; l<2; l++) {
         _quad->setLocal(Point<real_t> (_gauss[k],_gauss[l]));
         for (size_t i=1; i<=4; i++)
             p += _quad->getDet() * (_quad->DSh(i).x * ePrev(2*i-1)
                                   + _quad->DSh(i).y * ePrev(2*i  ));
      }
   }
   return -p*coef;
}


} /* namespace OFELI */
