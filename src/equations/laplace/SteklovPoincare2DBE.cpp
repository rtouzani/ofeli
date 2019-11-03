/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   This file is part of the OFELI Library

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

                    Implementation of class SteklovPoincare2DBE

  ==============================================================================*/


#include "equations/laplace/SteklovPoincare2DBE.h"
#include "solvers/GMRes.h"
#include "equations/AbsEqua_impl.h"
#include "shape_functions/Triang3.h"
#include "equations/AbsEqua_impl.h"
#include "linear_algebra/Vect_impl.h"


namespace OFELI {

SteklovPoincare2DBE::SteklovPoincare2DBE()
                    : AbsEqua<real_t>(), _ext(1)
{
   _A = nullptr;
   _b = nullptr;
}


SteklovPoincare2DBE::SteklovPoincare2DBE(Mesh& ms)
                    : AbsEqua<real_t>(), _ext(1)
{
   setMesh(ms);
}


SteklovPoincare2DBE::SteklovPoincare2DBE(Mesh&         ms,
                                         Vect<real_t>& u)
                    : AbsEqua<real_t>(), _ext(1)
{
   setMesh(ms);
   _u = &u;
}


SteklovPoincare2DBE::~SteklovPoincare2DBE()
{
}


void SteklovPoincare2DBE::setExterior()
{
   _ext = -1;
}


real_t SteklovPoincare2DBE::single_layer(size_t               j,
                                         const Point<real_t>& z) const
{
   return -0.125/OFELI_PI*_h*I2(0.25*_h*_h,z*_ttg[j],z.NNorm());
}


real_t SteklovPoincare2DBE::double_layer(size_t               j,
                                         const Point<real_t>& z) const
{
   real_t t = _nn[j]*z;
   return -0.25/OFELI_PI*t*I3(0.25*_h*_h,z*_ttg[j],z.NNorm()); 
}


void SteklovPoincare2DBE::setMesh(Mesh& ms)
{
   _theMesh = &ms;
   _nb_eq = _nb_boundary_sides = _theMesh->getNbBoundarySides();
   if (_nb_sides==0)
      throw OFELIException("SteklovPoincare2DBE::setMesh(ms): No boundary sides extracted.");
   _nn.setSize(_nb_eq);
   _length.setSize(_nb_eq);
   _center.setSize(_nb_eq);
   _ttg.setSize(_nb_eq);
   util();
   if (_A!=nullptr)
      delete _A;
   _A = new SpMatrix<real_t>(_nb_eq,_nb_eq);
   _ls.setMatrix(_A);
   _set_matrix = true;
   setSolver(GMRES_SOLVER,DIAG_PREC);
   if (_b!=nullptr)
      delete _b;
   _b = new Vect<real_t>(_nb_eq);
}


int SteklovPoincare2DBE::run()
{
   Side *sd1, *sd2;
   if (_bc==nullptr)
      throw OFELIException("SteklovPoincare2DBE::run(): No Dirichlet boundary condition given.");
   if (_u==nullptr)
      throw OFELIException("SteklovPoincare2DBE::run(): No solution vector given.");
   _u->clear();
   _A->clear();
   for (size_t s=0; s<_nb_boundary_sides; ++s) {
      sd1 = _theMesh->theBoundarySides[s];
      (*_b)[s] = -0.25*((*_bc)((*sd1)(1)->n())+(*_bc)((*sd1)(2)->n()));
      for (size_t t=0; t<_nb_boundary_sides; ++t) {
         sd2 = _theMesh->theBoundarySides[t];
         _h = _length[t];
         real_t g = 0.5*((*_bc)((*sd2)(1)->n())+(*_bc)((*sd2)(2)->n()));
         Point<real_t> z = _center[t] - _center[s];
         if (s==t)
            _A->add(s+1,s+1,0.25/OFELI_PI*_ext*_h*I1(0.5*_h));
         else {
            _A->add(s+1,t+1,-_ext*single_layer(t,z));
            _b->add(s+1,-_ext*g*double_layer(t,z));
         }
      }
   }
   int ret = solveLinearSystem(_A,*_b,*_u);
   return ret;
}


void SteklovPoincare2DBE::util()
{
   for (size_t s=0; s<_nb_boundary_sides; ++s) {
      the_side = _theMesh->theBoundarySides[s];
      the_element = the_side->getNeighborElement(1);
      if (the_element->getShape()!=TRIANGLE)
         throw OFELIException("SteklovPoincare2DBE::util(): This class is valid for triangles only.");
      Triang3 tr(the_element);
      Point<real_t> x1 = The_side(1)->getCoord(), x2 = The_side(2)->getCoord();
      _ttg[s] = x2 - x1;
      _center[s] = 0.5*(x1+x2);
      Point<real_t> c = _center[s] - tr.getCenter();
      Point<real_t> N(_ttg[s].y,-_ttg[s].x);
      if (c*N<0)
         N.x = -N.x, N.y = -N.y;
      _nn[s] = N;
      _length[s] = _ttg[s].Norm();
   }
}

} /* namespace OFELI */
