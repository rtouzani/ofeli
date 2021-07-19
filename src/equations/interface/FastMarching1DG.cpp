/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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


#include "equations/interface/FastMarching1DG.h"
#include "linear_algebra/Vect_impl.h"


namespace OFELI {

FastMarching1DG::FastMarching1DG()
                : _nx(0)
{
   _b = nullptr;
}


FastMarching1DG::FastMarching1DG(const Grid&   g,
                                 Vect<real_t>& T)
{
   set(g,T);
}


FastMarching1DG::FastMarching1DG(const Grid&   g,
                                 Vect<real_t>& T,
                                 Vect<real_t>& F)
{
   set(g,T);
   _b = &F;
}


FastMarching1DG::~FastMarching1DG()
{
   if (_b!=nullptr)
      delete _b;
}


void FastMarching1DG::set(const Grid&   g,
                          Vect<real_t>& T)
{
   _theGrid = &g;
   if (_theGrid->getDim()!=1)
      throw OFELIException("In FastMarching1DG::set(...): Wrong choice of space dimension");
   _nx = _theGrid->getNx();
   _hx = _theGrid->getHx();
   _u = &T;
   _b = new Vect<real_t>(_nx+1);
   _U.resize(_nx+1);
   *_b = 1.;
}


void FastMarching1DG::set(const Grid&   g,
                          Vect<real_t>& T,
                          Vect<real_t>& F)
{
   set(g,T);
   *_b = F;
}


size_t FastMarching1DG::Neigs()
{
   size_t i=_p->i;
   _neigs.clear();
   if (i>1 && _U[i-1].state!=Pt::FROZEN)
      _neigs.push_back(&_U[i-1]);
   if (i<=_nx && _U[i+1].state!=Pt::FROZEN)
      _neigs.push_back(&_U[i+1]);
   return _neigs.size();
}


void FastMarching1DG::init()
{
   for (size_t i=1; i<=_nx+1; ++i) {
      _p = &_U[i-1];
      _p->i = i;
      _p->x = _theGrid->getCoord(i).x;
      _p->sgn = Sgn((*_u)(i));
      _p->v = fabs((*_u)(i));
      _p->state = Pt::FAR;
      if (_p->v<INFINITY)
         _p->state = Pt::FROZEN;
   }
   for (size_t i=1; i<=_nx+1; ++i) {
      _p = &_U[i-1];
      if (_p->state==Pt::FROZEN) {
         int ng = Neigs();
         for (int n=0; n<ng; ++n) {
            _np = &_U[_neigs[n]->i-1];
            if (_np->state != Pt::FROZEN) {
               if (eval()==0) { 
                  if (_Narrow.find(_np)<0) {
                     _np->state = Pt::ALIVE;
                     _Narrow.insert(_np);
                  }
                  else {
                     _np->state = Pt::FROZEN;
                     _Narrow.remove();
                  }
               }
            }
         }
      }
   }
}


int FastMarching1DG::run()
{
   init();
   int nb = _Narrow.size();
   for (int n=0; n<_Narrow.size(); ++n) {
      _p = _Narrow[n];
      if (_p->state==Pt::FROZEN)
         nb--;
      _p->state = Pt::FROZEN;
      int ng = Neigs();
      for (int k=0; k<ng; ++k) {
         _np = &_U[_neigs[k]->i-1];
         if (eval()==0) {
            if (_Narrow.find(_np)<0) {
               _np->state = Pt::ALIVE;
               _Narrow.insert(_np);
               nb++;
            }
            else {
               if (_np->state!=Pt::FROZEN)
                  nb--;
               _np->state = Pt::FROZEN;
            }
         }
      }
   }
   for (size_t i=1; i<=_nx+1; ++i)
      (*_u)(i) = _U[i-1].v*_U[i-1].sgn;
   return 0;
}


real_t FastMarching1DG::eval()
{
   size_t i = _np->i;
   real_t vx=0., x=0.;
   real_t a = 1.0/(_hx*_hx);
   bool c1 = (i>1   ) && _U[i-2].state!=Pt::FAR;
   bool c2 = (i<=_nx) && _U[i  ].state!=Pt::FAR;
   if (c1 && !c2)
      vx = _U[i-1].v;
   else if (!c1 && c2)
      vx = _U[i+1].v;
   else if (!c1 && !c2)
      a = 0., vx = 0.;
   else if (i>1 && i<=_nx)
      vx = fmin(_U[i-2].v,_U[i].v);
   int ret = MaxQuad(a,-a*vx,a*vx*vx-1./((*_b)(i)*(*_b)(i)),x);
   if (ret==0)
      _np->v = x;
   return ret;
}


real_t FastMarching1DG::getResidual()
{
   real_t err=0., ux=0.;
   for (size_t i=1; i<=_nx; ++i) {
      ux = 1.0/_hx*((*_u)(i+1)-(*_u)(i));
      err += fabs(ux*ux - 1./((*_b)(i)*(*_b)(i)));
   }
   return err/_nx;
}

} /* namespace OFELI */
