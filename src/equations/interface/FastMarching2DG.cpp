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

                    Implementation of class 'FastMarching2DG'

  ==============================================================================*/


#include "equations/interface/FastMarching2DG.h"
#include "linear_algebra/Vect_impl.h"


namespace OFELI {

FastMarching2DG::FastMarching2DG()
                : _nx(0), _ny(0)
{
}


FastMarching2DG::FastMarching2DG(const Grid&   g,
                                 Vect<real_t>& T)
{
   set(g,T);
}


FastMarching2DG::FastMarching2DG(const Grid&   g,
                                 Vect<real_t>& T,
                                 Vect<real_t>& F)
{
   set(g,T,F);
   _b = F;
}


FastMarching2DG::~FastMarching2DG()
{
}


void FastMarching2DG::set(const Grid&   g,
                          Vect<real_t>& T)
{
   _theGrid = &g;
   if (_theGrid->getDim()!=2)
      throw OFELIException("In FastMarching2DG::set(...): Wrong choice of space dimension");
   _nx = _theGrid->getNx(), _ny = _theGrid->getNy();
   _hx = _theGrid->getHx(), _hy = _theGrid->getHy();
   _u = &T;
   _b.setSize(_nx+1,_ny+1);
   _U.resize((_nx+1)*(_ny+1));
   _b = 1.;
}


void FastMarching2DG::set(const Grid&   g,
                          Vect<real_t>& T,
                          Vect<real_t>& F)
{
   set(g,T);
   _b = F;
}


size_t FastMarching2DG::Neigs()
{
   size_t i=_p->i, j=_p->j;
   _neigs.clear();
   if (i>1 && _U[IJ(i-1,j)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJ(i-1,j)]);
   if (i<=_nx && _U[IJ(i+1,j)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJ(i+1,j)]);
   if (j>1 && _U[IJ(i,j-1)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJ(i,j-1)]);
   if (j<=_ny && _U[IJ(i,j+1)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJ(i,j+1)]);
   return _neigs.size();
}


void FastMarching2DG::init()
{
   int ret = 0;
   for (size_t i=1; i<=_nx+1; ++i) {
      for (size_t j=1; j<=_ny+1; ++j) {
         _p = &_U[IJ(i,j)];
         _p->i = i, _p->j = j;
         _p->x = _theGrid->getCoord(i,j).x;
         _p->y = _theGrid->getCoord(i,j).y;
         _p->sgn = Sgn((*_u)(i,j));
         _p->v = fabs((*_u)(i,j));
         _p->state = Pt::FAR;
         if (_p->v<INFINITY)
            _p->state = Pt::FROZEN;
      }
   }
   for (size_t i=1; i<=_nx+1; ++i) {
      for (size_t j=1; j<=_ny+1; ++j) {
         _p = &_U[IJ(i,j)];
         if (_p->state==Pt::FROZEN) {
            Neigs();
            for (auto const& N: _neigs) {
               _np = &_U[IJ(N->i,N->j)];
               if (_np->state != Pt::FROZEN) {
                  ret = eval();
                  if (ret==0) { 
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
}


int FastMarching2DG::run()
{
   init();
   int ret = 0;
   int nb = _Narrow.size();
   for (int n=0; n<_Narrow.size(); ++n) {
      _p = _Narrow[n];
      if (_p->state==Pt::FROZEN)
         nb--;
      _p->state = Pt::FROZEN;
      Neigs();
      for (auto const& N: _neigs) {
         _np = &_U[IJ(N->i,N->j)];
         ret = eval();
         if (ret==0) {
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
      for (size_t j=1; j<=_ny+1; ++j)
         (*_u)(i,j) = _U[IJ(i,j)].v*_U[IJ(i,j)].sgn;
   return 0;
}


real_t FastMarching2DG::eval()
{
   size_t i = _np->i, j = _np->j;
   real_t vx=0., vy=0., x=0.;
   real_t a = 1.0/(_hx*_hx), b = 1.0/(_hy*_hy);
   bool c1 = (i>1   ) && _U[IJ(i-1,j)].state!=Pt::FAR;
   bool c2 = (i<=_nx) && _U[IJ(i+1,j)].state!=Pt::FAR;
   bool d1 = (j>1   ) && _U[IJ(i,j-1)].state!=Pt::FAR;
   bool d2 = (j<=_ny) && _U[IJ(i,j+1)].state!=Pt::FAR;
   if (c1 && !c2)
      vx = _U[IJ(i-1,j)].v;
   else if (!c1 && c2)
      vx = _U[IJ(i+1,j)].v;
   else if (!c1 && !c2)
      a = 0., vx = 0.;
   else if (i>1 && i<=_nx)
      vx = fmin(_U[IJ(i-1,j)].v,_U[IJ(i+1,j)].v);
   if (d1 && !d2)
      vy = _U[IJ(i,j-1)].v;
   else if (!d1 && d2)
      vy = _U[IJ(i,j+1)].v;
   else if (!d1 && !d2)
      b = 0., vy = 0.;
   else if (j>1 && j<=_ny)
      vy = fmin(_U[IJ(i,j-1)].v,_U[IJ(i,j+1)].v);
   int ret = MaxQuad(a+b,-a*vx-b*vy,a*vx*vx+b*vy*vy-1./(_b(i,j)*_b(i,j)),x);
   if (ret==0)
      _np->v = x;
   return ret;
}


real_t FastMarching2DG::getResidual()
{
   real_t err=0., ux=0., uy=0.;
   for (size_t i=1; i<=_nx; ++i) {
      for (size_t j=1; j<=_ny; ++j) {
         ux = 0.5/_hx*((*_u)(i+1,j)-(*_u)(i,j) + (*_u)(i+1,j+1)-(*_u)(i,j+1));
         uy = 0.5/_hy*((*_u)(i,j+1)-(*_u)(i,j) + (*_u)(i+1,j+1)-(*_u)(i+1,j));
         err += fabs(ux*ux + uy*uy - 1./(_b(i,j)*_b(i,j)));
      }
   }
   return err/(_nx*_ny);
}

} /* namespace OFELI */
