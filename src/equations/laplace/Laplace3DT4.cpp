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

                      Class Laplace3DT4 : Laplace Equation
          using 4-Node Tetrahedral finite element in three dimensions

  ==============================================================================*/


#include "equations/laplace/Laplace3DT4.h"

namespace OFELI {

Laplace3DT4::Laplace3DT4(Mesh&             ms,
                         SpMatrix<real_t>& A,
                         Vect<real_t>&     b)
            : Equation<real_t,4,4,3,3>(ms)
{
   _u = &b;
}


Laplace3DT4::Laplace3DT4(Mesh&         ms,
                         Vect<real_t>& b)
            : Equation<real_t,4,4,3,3>(ms)
{
   _u = &b;
}


Laplace3DT4::Laplace3DT4(Mesh& ms)
            : Equation<real_t,4,4,3,3>(ms)
{  }


Laplace3DT4::Laplace3DT4(Element* el)
{
   set(el);
}


Laplace3DT4::Laplace3DT4(Side* sd)
{
   set(sd);
}


void Laplace3DT4::set(const Element* el)
{
   _nb_dof = 1;
   Init(el);
   Tetra4 t(_theElement);
   _volume = t.getVolume();
   _det = t.getDet();
   for (size_t i=1; i<=4; i++)
      _dSh(i) = t.DSh(i);
   eMat = 0;
   eRHS = 0;
   eA0 = 0;
   eA1 = 0;
}


void Laplace3DT4::set(const Side* sd)
{
   _nb_dof = 1;
   Init(sd);
   Triang3 t(_theSide);
   _area = t.getArea();
   sMat = 0;
   sRHS = 0;
}


void Laplace3DT4::build()
{
   MESH_EL {
      set(theElement);
      LHS();
      BodyRHS(*_bf);
      if (_bc_given)
         updateBC(*_bc);
      element_assembly(*theElement,eMat,_A);
      ElementAssembly(*_u);
   }
   MESH_SD {
      set(theSide);
      BoundaryRHS(*_sf);
      SideAssembly(*_u);
   }
}


void Laplace3DT4::buildEigen(int opt)
{
   set(theElement);
   for (size_t i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eA0(i,j) += _area*_dSh(i)*_dSh(j);
   if (opt==0) {
      real_t c = 0.1*_volume;
      real_t d = 0.5*c;
      eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c; eA1(4,4) += c;
      eA1(1,2) += d; eA1(2,1) += d; eA1(1,3) += d; eA1(1,4) += d;
      eA1(3,1) += d; eA1(2,3) += d; eA1(3,2) += d; eA1(2,4) += d;
      eA1(4,1) += d; eA1(4,2) += d; eA1(4,3) += d; eA1(3,4) += d;
   }
   else {
      real_t c = 0.25*_volume;
      for (size_t i=1; i<=4; i++)
         eA1(i,i) += c;
   }
}


int Laplace3DT4::solve(Vect<real_t>& u)
{
   Vect<real_t> x(_A.size());
   LinearSolver<real_t> ls(_A,*_u,x);
   ls.setTolerance(1.e-7);
   int nb_it = ls.solve(CG_SOLVER,ILU_PREC);
   u.insertBC(*_theMesh,x,*_bc);
   return nb_it;
}


void Laplace3DT4::LHS(real_t coef)
{
   for (size_t i=1; i<=4; i++)
      for (size_t j=1; j<=4; j++)
         eMat(i,j) += coef*_volume*_dSh(i)*_dSh(j);
}


void Laplace3DT4::BodyRHS(const Vect<real_t>& f)
{
   real_t c = 0.25*_volume;
   for (size_t i=1; i<=4; i++)
      eRHS(i) += f((*_theElement)(i)->n())*c;
}


void Laplace3DT4::BoundaryRHS(const Vect<real_t>& h)
{
   if (_theSide->getCode(1)>0) {
      real_t c = OFELI_THIRD*_area;
      sRHS(1) += c*h((*_theSide)(1)->n());
      sRHS(2) += c*h((*_theSide)(2)->n());
      sRHS(3) += c*h((*_theSide)(3)->n());
   }
}


void Laplace3DT4::Axb(const Vect<real_t>& x,
                      Vect<real_t>&       b)
{
   MESH_EL {
      set(theElement);
      LHS();
      BodyRHS(*_bf);
      ElementAssembly(*_u);
      AxbAssembly(TheElement,x,b);
   }
   MESH_SD {
      set(theSide);
      BoundaryRHS(*_sf);
      SideAssembly(*_u);
      AxbAssembly(TheElement,x,b);
   }
}

} /* namespace OFELI */
