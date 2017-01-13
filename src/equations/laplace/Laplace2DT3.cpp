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

                      Class Laplace2DT3: Laplace Equation
              using 3-Node Triangular finite element in two dimensions

  ==============================================================================*/


#include "equations/laplace/Laplace2DT3.h"

namespace OFELI {

Laplace2DT3::Laplace2DT3(Mesh&             ms,
                         SpMatrix<real_t>& A,
                         Vect<real_t>&     b)
            : Equation<real_t,3,3,2,2>(ms)
{
   _u = &b;
   _A = &A;
}


Laplace2DT3::Laplace2DT3(Mesh&         ms,
                         Vect<real_t>& b)
            : Equation<real_t,3,3,2,2>(ms)
{
   _u = &b;
}


Laplace2DT3::Laplace2DT3(Mesh& ms)
            : Equation<real_t,3,3,2,2>(ms)
{
}


Laplace2DT3::Laplace2DT3(Mesh&         ms,
                         Vect<real_t>& b,
                         Vect<real_t>& Dbc,
                         Vect<real_t>& Nbc,
                         Vect<real_t>& u)
            : Equation<real_t,3,3,2,2>(ms)
{
   _bf = &b;
   _u = &u;
   _bc_given = true;
   _bc = &Dbc;
   _sf = &Nbc;
   _A = new SpMatrix<real_t>(*_theMesh);
}


Laplace2DT3::Laplace2DT3(Element* el)
{
   set(el);
}


Laplace2DT3::Laplace2DT3(Side* sd)
{
   set(sd);
}


void Laplace2DT3::set(const Element* el)
{
   _nb_dof = 1;
   Init(el);
   Triang3 tr(_theElement);
   _area = tr.getArea();
   _dSh(1) = tr.DSh(1); _dSh(2) = tr.DSh(2); _dSh(3) = tr.DSh(3);
   eMat = 0;
   eRHS = 0;
   eA0 = 0;
   eA1 = 0;
}


void Laplace2DT3::set(const Side* sd)
{
   _nb_dof = 1;
   Init(sd);
   Line2 ln(_theSide);
   _length = ln.getLength();
   sMat = 0;
   sRHS = 0;
}


void Laplace2DT3::build()
{
   *_A = 0;
   _b = new Vect<real_t>(_A->size());
   MESH_EL {
      set(theElement);
      LHS();
      BodyRHS(*_bf);
      if (_bc_given)
         updateBC(*_bc);
      ElementAssembly(*_A);
      ElementAssembly(*_b);
   }
   MESH_SD {
      set(theSide);
      BoundaryRHS(*_sf);
      SideAssembly(*_b);
   }
}


void Laplace2DT3::buildEigen(int opt)
{
   set(theElement);
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eA0(i,j) += _area*(_dSh(i),_dSh(j));
   if (opt==0) {
      real_t c = 0.5*OFELI_SIXTH*_area;
      eA1(1,1) += 2*c; eA1(2,2) += 2*c; eA1(3,3) += 2*c;
      eA1(1,2) +=   c; eA1(2,1) +=   c; eA1(1,3) +=   c;
      eA1(3,1) +=   c; eA1(2,3) +=   c; eA1(3,2) +=   c;
   }
   else {
      real_t c = OFELI_THIRD*_area;
      eA1(1,1) += c; eA1(2,2) += c; eA1(3,3) += c;
   }
}


int Laplace2DT3::solve(Vect<real_t>& u)
{
   Vect<real_t> x(_b->size());
   LinearSolver<real_t> ls(*_A,*_b,x);
   ls.setTolerance(1.e-7);
   int nb_it = ls.solve(CG_SOLVER,DILU_PREC);
   u.insertBC(*_theMesh,x,*_bc);
   delete _b;
   return nb_it;
}


int Laplace2DT3::run()
{
   build();
   int n = solve(*_u);
   delete _A;
   return n;
}


void Laplace2DT3::LHS(real_t coef)
{
   for (size_t i=1; i<=3; i++)
      for (size_t j=1; j<=3; j++)
         eMat(i,j) += coef*_area*(_dSh(i),_dSh(j));
}


void Laplace2DT3::BodyRHS(const Vect<real_t>& f)
{
   eRHS(1) += f((*_theElement)(1)->n())*_area*OFELI_THIRD;
   eRHS(2) += f((*_theElement)(2)->n())*_area*OFELI_THIRD;
   eRHS(3) += f((*_theElement)(3)->n())*_area*OFELI_THIRD;
}


void Laplace2DT3::BoundaryRHS(const Vect<real_t>& h)
{
   if (_theSide->getCode(1)>0) {
      sRHS(1) += 0.5*_length*h((*_theSide)(1)->n());
      sRHS(2) += 0.5*_length*h((*_theSide)(2)->n());
   }
}


void Laplace2DT3::Axb(const Vect<real_t>& x,
                            Vect<real_t>& b)
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
