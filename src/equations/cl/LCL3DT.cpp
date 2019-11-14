/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                        Implementation of class LCL3DT

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/


#include "equations/cl/LCL3DT.h"
#include "mesh/Mesh.h"
#include "mesh/Element.h"
#include "mesh/Side.h"
#include "linear_algebra/LocalVect_impl.h"
#include "linear_algebra/Vect_impl.h"
#include <algorithm>

using std::min;
using std::max;

namespace OFELI {

LCL3DT::LCL3DT(Mesh& m) : Muscl3DT(m)
{
   _init_alloc = true;
   init();
   _U = new Vect<real_t>(*_theMesh);
}


LCL3DT::LCL3DT(Mesh&         m,
               Vect<real_t>& U)
       : Muscl3DT(m)
{
   _init_alloc = false;
   init();
   _U = &U;
}


LCL3DT::~LCL3DT()
{
   if (_init_alloc)
      delete _U;
}


void LCL3DT::init()
{
   _LU.setMesh(*_theMesh);
   _RU.setMesh(*_theMesh);
   _FU.setMesh(*_theMesh);
   _v.setMesh(*_theMesh,3);
   setReferenceLength(getMinimumEdgeLength());
   setCFL(0.2);
}


void LCL3DT::setInitialCondition(Vect<real_t>& u)
{
   _U = &u;
}


void LCL3DT::setInitialCondition(real_t u)
{
    *_U = u;
}


void LCL3DT::setReconstruction()
{
   if (Muscl3DT::setReconstruction(*_U,_LU,_RU,1)) {
      cerr << "ERROR: reconstruction of u failed" << endl;
      exit(3);
   }
}


void LCL3DT::setBC(const Side&  sd,
                         real_t u)
{
   _RU(sd.n()) = u;
}


void LCL3DT::setBC(int    code,
                   real_t u)
{
   MeshBoundarySides(*_theMesh)
      if (theSide->getCode()==code)
         _RU(theSideLabel) = u;
}


void LCL3DT::setBC(real_t u)
{
   MeshBoundarySides(*_theMesh)
      _RU(theSideLabel) = u;
}


void LCL3DT::setVelocity(const LocalVect<real_t,3>& v)
{
   MESH_SD {
      size_t s = theSideLabel;
      _v(s,1) = v[0];
      _v(s,2) = v[1];
      _v(s,3) = v[2];
   }
}


void LCL3DT::setVelocity(const Vect<real_t>& v)
{
   _v = v;
}


real_t LCL3DT::getRiemannUpwind(size_t s)
{
   real_t flux, lambda;
   real_t zeta = _n(s).x*_v(s,1) + _n(s).y*_v(s,2) + _n(s).z*_v(s,3);
   (zeta>0.) ? flux = zeta*_LU(s,1) : flux = zeta*_RU(s,1);
   _FU(s) = flux;
   lambda = sqrt(_v(s,1)*_v(s,1)  + _v(s,2)*_v(s,2) + _v(s,3)*_v(s,3));
   return lambda;
}


real_t LCL3DT::getFluxUpwind()
{
   real_t lambda=1.e-10;
   MESH_SD
      lambda = max(lambda,getRiemannUpwind(theSideLabel));
   return lambda;
}


void LCL3DT::forward()
{
   MESH_SD {
      size_t n = theSideLabel;
      (*_U)(theSide->getNeighborElement(1)->n()) -= _FU(n)*_Lrate(n)*_TimeStep;
      if (!theSide->isOnBoundary())
         (*_U)(theSide->getNeighborElement(2)->n()) += _FU(n)*_Rrate(n)*_TimeStep;
   }
}


void LCL3DT::Forward(const Vect<real_t>& Flux,
                           Vect<real_t>& Field)
{
   MESH_SD {
      size_t n = theSideLabel;
      Field(theSide->getNeighborElement(1)->n()) -= Flux(n)*_Lrate(n)*_TimeStep;
      if (!theSide->isOnBoundary())
         Field(theSide->getNeighborElement(2)->n()) += Flux(n)*_Rrate(n)*_TimeStep;
   }
}


real_t LCL3DT::runOneTimeStep()
{
   real_t V = getFluxUpwind();
   _TimeStep = _ReferenceLength/V*_CFL;
   forward();
   return _TimeStep;
}

} /* namespace OFELI */
