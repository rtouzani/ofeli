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

                          Implementation of class LCL2DT

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/


#include "equations/cl/LCL2DT.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/LocalVect.h"
#include "linear_algebra/LocalMatrix.h"
#include "linear_algebra/Point.h"
#include <algorithm>

using std::min;
using std::max;

namespace OFELI {


LCL2DT::LCL2DT(Mesh& m)
       : Muscl2DT(m)
{
   init();
   _init_alloc = true;
   _U = new Vect<real_t>(*_theMesh,1,ELEMENT_DOF);
}


LCL2DT::LCL2DT(Mesh&         m,
               Vect<real_t>& U)
       : Muscl2DT(m)
{
   _init_alloc = false;
   init();
   _U = &U;
}


void LCL2DT::init()
{
   _LU.setSize(_nb_sides);
   _RU.setSize(_nb_sides);
   _FU.setSize(_nb_sides);
   _v.setSize(_nb_sides,2);
   setReferenceLength(getMinimumEdgeLength());
   setCFL(0.2);
   setTimeStep(0.);
}


LCL2DT::~LCL2DT()
{
   if (_init_alloc)
      delete _U;
}


void LCL2DT::setInitialCondition(Vect<real_t>& u)
{
   *_U = u;
}


void LCL2DT::setInitialCondition(real_t u)
{
   *_U = u;
}


void LCL2DT::setReconstruction()
{
   if (Muscl2DT::setReconstruction(*_U,_LU,_RU,1)) {
      cerr << "ERROR: Reconstruction of u failed" << endl;
      exit(3);
   }
};


void LCL2DT::setBC(const Side&  sd,
                         real_t u)
{
   _RU.set(sd.n(),u);
}


void LCL2DT::setBC(int    code,
                   real_t u)
{
   mesh_sides(*_theMesh) {
      if (The_side.isOnBoundary() && The_side.getCode()==code)
         _RU.set(side_label,u);
   }
}


void LCL2DT::setBC(real_t u)
{
   mesh_sides(*_theMesh) {
      if (The_side.isOnBoundary())
         _RU.set(side_label,u);
   }
}


void LCL2DT::setVelocity(const LocalVect<real_t,2>& v)
{
   _v.setSize(_nb_sides,2);
   mesh_sides(*_theMesh) {
      size_t n = side_label;
      _v.set(n,1,v(1));
      _v.set(n,2,v(2));
   }
}


void LCL2DT::setVelocity(const Vect<real_t>& v)
{
   _v.setSize(_nb_sides,2);
   mesh_sides(*_theMesh) {
      size_t n = side_label;
      _v.set(n,1,v(n,1));
      _v.set(n,2,v(n,2));
   }
}


real_t LCL2DT::getRiemannUpwind(size_t s)
{
   real_t zeta = _n(s).x*_v(s,1) + _n(s).y*_v(s,2);
   (zeta>0.) ? _FU.set(s,zeta*_LU(s)) : _FU.set(s,zeta*_RU(s));
   return sqrt(_v(s,1)*_v(s,1) + _v(s,2)*_v(s,2));
}


real_t LCL2DT::getFluxUpwind()
{
   real_t lambda=1.e-10;
   mesh_sides(*_theMesh)
      lambda = max(getRiemannUpwind(side_label),lambda);
   return lambda;
}


void LCL2DT::forward()
{
   mesh_sides(*_theMesh) {
      size_t n = side_label;
      _U->add(The_side.getNeighborElement(1)->n(),-_FU(n)*_Lrate(n)*_TimeStep);
      if (The_side.isOnBoundary()==false)
         _U->add(The_side.getNeighborElement(2)->n(),_FU(n)*_Rrate(n)*_TimeStep);
   }
}


void LCL2DT::Forward(const Vect<real_t>& Flux,
                           Vect<real_t>& Field)
{
   mesh_sides(*_theMesh) {
      size_t n = side_label;
      Field.add(The_side.getNeighborElement(1)->n(),-Flux(n)*_Lrate(n)*_TimeStep);
      if (The_side.isOnBoundary()==false)
         Field.add(The_side.getNeighborElement(2)->n(),Flux(n)*_Rrate(n)*_TimeStep);
   }
}


real_t LCL2DT::runOneTimeStep()
{
   real_t V = getFluxUpwind();
   _TimeStep = _ReferenceLength/V*_CFL;
   forward();
   return _TimeStep;
}

} /* namespace OFELI */
