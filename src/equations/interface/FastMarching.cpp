/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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

                     Implementation of class 'FastMarching'

  ==============================================================================*/


#include "equations/interface/FastMarching.h"
#include "linear_algebra/Vect_impl.h"


namespace OFELI {

FastMarching::FastMarching()
             : _theFMM1D(nullptr), _theFMM2D(nullptr), _theFMM3D(nullptr)
{
}


FastMarching::FastMarching(const Grid&   g,
                           Vect<real_t>& T)
             : _theFMM1D(nullptr), _theFMM2D(nullptr), _theFMM3D(nullptr)
{
   _dim = g.getDim();
   if (_dim==1)
      _theFMM1D = new FastMarching1DG(g,T);
   else if (_dim==2)
      _theFMM2D = new FastMarching2DG(g,T);
   else if (_dim==3)
      _theFMM3D = new FastMarching3DG(g,T);
}


FastMarching::FastMarching(const Grid&   g,
                           Vect<real_t>& T,
                           Vect<real_t>& F)
             : _theFMM1D(nullptr), _theFMM2D(nullptr), _theFMM3D(nullptr)
{
   _dim = g.getDim();
   if (_dim==1)
      _theFMM1D = new FastMarching1DG(g,T,F);
   else if (_dim==2)
      _theFMM2D = new FastMarching2DG(g,T,F);
   else if (_dim==3)
      _theFMM3D = new FastMarching3DG(g,T,F);
}


FastMarching::~FastMarching()
{
   if (_theFMM1D!=nullptr)
      delete _theFMM1D;
   else if (_theFMM2D!=nullptr)
      delete _theFMM2D;
   else if (_theFMM3D!=nullptr)
      delete _theFMM3D;
}


void FastMarching::set(const Grid&   g,
                       Vect<real_t>& T)
{
   if (_dim==1 && _theFMM1D==nullptr) {
      _theFMM1D = new FastMarching1DG;
      _theFMM1D->set(g,T);
   }
   else if (_dim==2 && _theFMM2D==nullptr) {
      _theFMM2D = new FastMarching2DG;
      _theFMM2D->set(g,T);
   }
   else if (_dim==3 && _theFMM3D==nullptr) {
      _theFMM3D = new FastMarching3DG;
      _theFMM3D->set(g,T);
   }
}


void FastMarching::set(const Grid&   g,
                       Vect<real_t>& T,
                       Vect<real_t>& F)
{
   if (_dim==1 && _theFMM1D==nullptr) {
      _theFMM1D = new FastMarching1DG;
      _theFMM1D->set(g,T,F);
   }
   else if (_dim==2 && _theFMM2D==nullptr) {
      _theFMM2D = new FastMarching2DG;
      _theFMM2D->set(g,T,F);
   }
   else if (_dim==3 && _theFMM3D==nullptr) {
      _theFMM3D = new FastMarching3DG;
      _theFMM3D->set(g,T,F);
   }
}


int FastMarching::run()
{
   int ret = 0;
   if (_dim==1)
      ret = _theFMM1D->run();
   else if (_dim==2)
      ret = _theFMM2D->run();
   else if (_dim==3)
      ret = _theFMM3D->run();
   return ret;
}


real_t FastMarching::getResidual()
{
   real_t res = 0.;
   if (_dim==1)
      res = _theFMM1D->getResidual();
   else if (_dim==2)
      res = _theFMM2D->getResidual();
   else if (_dim==3)
      res = _theFMM3D->getResidual();
   return res;
}

} /* namespace OFELI */
