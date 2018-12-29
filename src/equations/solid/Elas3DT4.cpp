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

                         Implementation of class ELAS3DT4
  for 3-D Linear Elasticity Equations using 4-node tetrahedral P1 finite element

  ==============================================================================*/


#include "equations/solid/Elas3DT4.h"
#include "shape_functions/Tetra4.h"
#include "shape_functions/Triang3.h"

namespace OFELI {

Elas3DT4::Elas3DT4(const Element* el)
{
   set(el);
}


Elas3DT4::Elas3DT4(const Side* sd)
{
   set(sd);
}


Elas3DT4::Elas3DT4(const Element*      el,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   set(el);
   ElementVector(u);
}


Elas3DT4::Elas3DT4(const Side*         sd,
                   const Vect<real_t>& u,
                   const real_t&       time)
{
   set(sd);
   SideVector(u);
}


Elas3DT4::~Elas3DT4() { }


void Elas3DT4::set(const Element* el)
{
   _nb_dof = 3;
   Init(el);
   setMaterial();
   Tetra4 tetra(el);
   _det = tetra.getDet();
   _center = tetra.getCenter();
   ElementNodeCoordinates();
   _dSh(1) = tetra.DSh(1);
   _dSh(2) = tetra.DSh(2);
   _dSh(3) = tetra.DSh(3);
   _dSh(4) = tetra.DSh(4);
   _lambda = _nu*_E/((1+_nu)*(1-2*_nu));
   _G = 0.5*_E/(1+_nu);
   eMat = 0;
   eRHS = 0;
}


void Elas3DT4::set(const Side* sd)
{
   _nb_dof = 3;
   Init(sd);
   Triang3 tr(sd);
   _area = tr.getArea();
   SideNodeCoordinates();
   sMat = 0;
   sRHS = 0;
}


void Elas3DT4::Media(real_t E,
                     real_t nu,
                     real_t rho)
{
   _E = E;
   _nu = nu;
   _rho = rho;
}


void Elas3DT4::LMassToLHS(real_t coef)
{
   real_t c = 0.5*coef*OFELI_TWELVETH*_det*_rho;
   for (size_t i=1; i<=4; i++) {
      eMat(3*i-2,3*i-2) += c;
      eMat(3*i-1,3*i-1) += c;
      eMat(3*i  ,3*i  ) += c;
   }
}


void Elas3DT4::LMassToRHS(real_t coef)
{
   for (size_t i=1; i<=4; i++)
      eRHS(i) += 0.5*coef*OFELI_TWELVETH*_det*_rho*ePrev(i);
}


void Elas3DT4::Deviator(real_t coef)
{
   real_t c=_G*_det*coef;
   for (size_t j=1; j<=4; j++) {
      Point<real_t> db=c*_dSh(j);
      for (size_t i=1; i<=4; i++) {
         eMat(3*i-2,3*j-2) += 2*_dSh(i).x*db.x + _dSh(i).z*db.z + _dSh(i).y*db.y;
         eMat(3*i-2,3*j-1) += _dSh(i).y*db.x;
         eMat(3*i-2,3*j  ) += _dSh(i).z*db.x;
         eMat(3*i-1,3*j-2) += _dSh(i).x*db.y;
         eMat(3*i-1,3*j-1) += 2*_dSh(i).y*db.y + _dSh(i).z*db.z + _dSh(i).x*db.x;
         eMat(3*i-1,3*j  ) += _dSh(i).z*db.y;
         eMat(3*i  ,3*j-2) += _dSh(i).x*db.z;
         eMat(3*i  ,3*j-1) += _dSh(i).y*db.z;
         eMat(3*i  ,3*j  ) += 2*_dSh(i).z*db.z + _dSh(i).y*db.y + _dSh(i).x*db.x;
      }
   }
}


void Elas3DT4::DeviatorToRHS(real_t coef)
{
   real_t c1=_G*_det*coef, c2=2*c1;
   for (size_t j=1; j<=4; j++) {
      real_t db11 = c2*_dSh(j).x;
      real_t db22 = c2*_dSh(j).y;
      real_t db33 = c2*_dSh(j).z;
      real_t db42 = c1*_dSh(j).z;
      real_t db43 = c1*_dSh(j).y;
      real_t db53 = c1*_dSh(j).x;
      for (size_t i=1; i<=4; i++) {
         eRHS(3*i-2) -= (_dSh(i).x*db11+_dSh(i).z*db42+_dSh(i).y*db43)*ePrev(3*j-2)
                       + _dSh(i).y*db53*ePrev(3*j-1) + _dSh(i).z*db53*ePrev(3*j);
         eRHS(3*i-1) -= _dSh(i).x*db43*ePrev(3*j-2) + _dSh(i).z*db43*ePrev(3*j)
                      + (_dSh(i).y*db22+_dSh(i).z*db42+_dSh(i).x*db53)*ePrev(3*j-1);
         eRHS(3*i  ) -= _dSh(i).x*db42*ePrev(3*j-2) + _dSh(i).y*db42*ePrev(3*j-1)
                      + (_dSh(i).z*db33+_dSh(i).y*db43+_dSh(i).x*db53)*ePrev(3*j);
      }
   }
}


void Elas3DT4::Dilatation(real_t coef)
{
   real_t c = _lambda*_det*coef;
   for (size_t j=1; j<=4; j++) {
      for (size_t i=1; i<=4; i++) {
         eMat(3*i-2,3*j-2) += c*_dSh(i).x*_dSh(j).x;
         eMat(3*i-2,3*j-1) += c*_dSh(i).x*_dSh(j).y;
         eMat(3*i-2,3*j  ) += c*_dSh(i).x*_dSh(j).z;
         eMat(3*i-1,3*j-2) += c*_dSh(i).y*_dSh(j).x;
         eMat(3*i-1,3*j-1) += c*_dSh(i).y*_dSh(j).y;
         eMat(3*i-1,3*j  ) += c*_dSh(i).y*_dSh(j).z;
         eMat(3*i  ,3*j-2) += c*_dSh(i).z*_dSh(j).x;
         eMat(3*i  ,3*j-1) += c*_dSh(i).z*_dSh(j).y;
         eMat(3*i  ,3*j  ) += c*_dSh(i).z*_dSh(j).z;
      }
   }
}


void Elas3DT4::DilatationToRHS(real_t coef)
{
   real_t c = _lambda*_det*coef;
   for (size_t j=1; j<=4; j++) {
      for (size_t i=1; i<=4; i++) {
         eRHS(3*j-2) -= c*(_dSh(i).x*_dSh(j).x*ePrev(3*j-1) +
                           _dSh(i).x*_dSh(j).y*ePrev(3*j-1) +
                           _dSh(i).x*_dSh(j).z*ePrev(3*j));
         eRHS(3*j-1) -= c*(_dSh(i).y*_dSh(j).x*ePrev(3*j-1) +
                           _dSh(i).y*_dSh(j).y*ePrev(3*j-1) +
                           _dSh(i).y*_dSh(j).z*ePrev(3*j));
         eRHS(3*j  ) -= c*(_dSh(i).z*_dSh(j).x*ePrev(3*j-1) +
                           _dSh(i).z*_dSh(j).y*ePrev(3*j-1) +
                           _dSh(i).z*_dSh(j).z*ePrev(3*j));
      }
   }
}


void Elas3DT4::BodyRHS(UserData<real_t>& ud)
{
   real_t fx = ud.BodyForce(_center,_time,1);
   real_t fy = ud.BodyForce(_center,_time,2);
   real_t fz = ud.BodyForce(_center,_time,3);
   for (size_t i=1; i<=4; i++) {
      eRHS(3*i-2) += fx*0.5*OFELI_TWELVETH*_det;
      eRHS(3*i-1) += fy*0.5*OFELI_TWELVETH*_det;
      eRHS(3*i  ) += fz*0.5*OFELI_TWELVETH*_det;
   }
}


void Elas3DT4::BodyRHS(const Vect<real_t>& f,
                             int           opt)
{
   if (opt==LOCAL_ARRAY) {
      for (size_t i=1; i<=4; i++) {
         eRHS(3*i-2) += 0.5*OFELI_TWELVETH*_det*f(i,1);
         eRHS(3*i-1) += 0.5*OFELI_TWELVETH*_det*f(i,2);
         eRHS(3*i  ) += 0.5*OFELI_TWELVETH*_det*f(i,3);
      }
   }
   else {
      for (size_t i=1; i<=4; i++) {
         eRHS(3*i-2) += 0.5*OFELI_TWELVETH*_det*f((*_theElement)(i)->n(),1);
         eRHS(3*i-1) += 0.5*OFELI_TWELVETH*_det*f((*_theElement)(i)->n(),2);
         eRHS(3*i  ) += 0.5*OFELI_TWELVETH*_det*f((*_theElement)(i)->n(),3);
      }
   }
}


void Elas3DT4::BoundaryRHS(const Vect<real_t>& f)
{
   for (size_t k=1; k<=3; k++) {
      if (_theSide->getCode(k) != CONTACT) {
         real_t c = _area*f(_theSide->n(),k);
         sRHS(k  ) += c;
         sRHS(k+3) += c;
         sRHS(k+6) += c;
      }
   }
}


void Elas3DT4::buildEigen(SkSMatrix<real_t>& K,
                          Vect<real_t>&      M)
{
   size_t i, j;
   MESH_EL {
      set(theElement);
      real_t c1=_G*_det, c2=2*c1;
      for (j=1; j<=4; j++) {
         real_t db11 = c2*_dSh(j).x, db22 = c2*_dSh(j).y;
         real_t db33 = c2*_dSh(j).z, db42 = c1*_dSh(j).z;
         real_t db43 = c1*_dSh(j).y, db53 = c1*_dSh(j).x;
         for (i=1; i<=4; i++) {
            eMat(3*i-2,3*j-2) += _dSh(i).x*db11 + _dSh(i).z*db42 + _dSh(i).y*db43;
            eMat(3*i-2,3*j-1) += _dSh(i).y*db53;
            eMat(3*i-2,3*j  ) += _dSh(i).z*db53;
            eMat(3*i-1,3*j-2) += _dSh(i).x*db43;
            eMat(3*i-1,3*j-1) += _dSh(i).y*db22 + _dSh(i).z*db42 + _dSh(i).x*db53;
            eMat(3*i-1,3*j  ) += _dSh(i).z*db43;
            eMat(3*i  ,3*j-2) += _dSh(i).x*db42;
            eMat(3*i  ,3*j-1) += _dSh(i).y*db42;
            eMat(3*i  ,3*j  ) += _dSh(i).z*db33 + _dSh(i).y*db43 + _dSh(i).x*db53;
         }
      }
      real_t c = _lambda*_det;
      for (j=1; j<=4; j++) {
         for (i=1; i<=4; i++) {
            eMat(3*i-2,3*j-2) += c*_dSh(i).x*_dSh(j).x;
            eMat(3*i-2,3*j-1) += c*_dSh(i).x*_dSh(j).y;
            eMat(3*i-2,3*j  ) += c*_dSh(i).x*_dSh(j).z;
            eMat(3*i-1,3*j-2) += c*_dSh(i).y*_dSh(j).x;
            eMat(3*i-1,3*j-1) += c*_dSh(i).y*_dSh(j).y;
            eMat(3*i-1,3*j  ) += c*_dSh(i).y*_dSh(j).z;
            eMat(3*i  ,3*j-2) += c*_dSh(i).z*_dSh(j).x;
            eMat(3*i  ,3*j-1) += c*_dSh(i).z*_dSh(j).y;
            eMat(3*i  ,3*j  ) += c*_dSh(i).z*_dSh(j).z;
         }
      }
      c = 0.5*OFELI_TWELVETH*_det*_rho;
      for (i=1; i<=4; i++) {
         eRHS(3*i-2) += c;
         eRHS(3*i-1) += c;
         eRHS(3*i  ) += c;
      }
      ElementAssembly(K);
      M.Assembly(TheElement,eRHS.get());
   }
}

} /* namespace OFELI */
