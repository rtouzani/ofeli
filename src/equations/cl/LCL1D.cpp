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

                          Implementation of class LCL1D

       Author: S. Clain
               clain@mip.ups-tlse.fr

  ==============================================================================*/


#include "equations/cl/LCL1D.h"
#include <algorithm>

using std::min;
using std::max;

namespace OFELI {


LCL1D::LCL1D(Mesh &m) : Muscl1D(m)
{
   _Init();
   _init_alloc = true;
   _U = new Vect<double>(*_theMesh);
}


LCL1D::LCL1D(Mesh &m, Vect<double> &U) : Muscl1D(m)
{
   _Init();
   _U = &U;
   _init_alloc = false;
}


void LCL1D::_Init()
{
   _LU.setMesh(*_theMesh);
   _RU.setMesh(*_theMesh);
   _FU.setMesh(*_theMesh);
   _v.setMesh(*_theMesh);
   setReferenceLength(getMinimumLength());
   setCFL(0.2);
   setTimeStep(0.);
}


LCL1D::~LCL1D()
{
   if (_init_alloc)
      delete _U;
}


void LCL1D::setInitialCondition(Vect<double> &u)
{
   *_U = u;
}


void LCL1D::setInitialCondition(double u)
{
   *_U = u;
}


void LCL1D::setReconstruction()
{
   if (Muscl1D::setReconstruction(*_U,_LU,_RU,1)) {
      cerr << "ERROR: reconstruction of u failed" << endl;
      exit(3);
   }
};


void LCL1D::setBC(const Side &sd, double v)
{
   _RU(sd.n()) = v;
}


void LCL1D::setBC(int code, double u)
{
   MeshBoundarySides(*_theMesh) {
      if ((TheSide.getCode()==code))
         _RU(theSideLabel) = u;
   }
}


void LCL1D::setBC(double u)
{
   MeshBoundarySides(*_theMesh)
      _RU(theSideLabel) = u;
}


void LCL1D::setVelocity(double v)
{
   MESH_SD
      _v(theSideLabel) = v;
}


void LCL1D::setVelocity(Vect<double> &v)
{
   _v = v;
}


double LCL1D::computeRiemannUpwind(size_t s)
{
   double zeta, flux, lambda;
   zeta = _v(s);
   (zeta>0.) ? flux = zeta*_LU(s) : flux = zeta*_RU(s);
   _FU(s) = flux;
   lambda = fabs(_v(s));
   return lambda;
}


double LCL1D::getFluxUpwind()
{
   double lambda_max=1.e-9;
   MESH_SD
      lambda_max = max(computeRiemannUpwind(theSideLabel),lambda_max);
   return lambda_max;
}


void LCL1D::forward()
{
   MESH_SD {
      size_t n = theSideLabel;
      (*_U)(theSide->getNeighborElement(1)->n()) -= _FU(n)*_Lrate(n)*_TimeStep;
      if (!theSide->isOnBoundary())
         (*_U)(theSide->getNeighborElement(2)->n()) += _FU(n)*_Rrate(n)*_TimeStep;
   }
}


double LCL1D::runOneTimeStep()
{
   double V = getFluxUpwind();
   _TimeStep = _ReferenceLength/V*_CFL;
   forward();
   return _TimeStep;
}


void LCL1D::Forward(const Vect<double> &Flux, Vect<double> &Field)
{
   MESH_SD {
      size_t n = theSideLabel;
      Field(theSide->getNeighborElement(1)->n()) -= Flux(n)*_Lrate(n)*_TimeStep;
      if (!theSide->isOnBoundary())
         Field(theSide->getNeighborElement(2)->n()) += Flux(n)*_Rrate(n)*_TimeStep;
   }
}

} /* namespace OFELI */
