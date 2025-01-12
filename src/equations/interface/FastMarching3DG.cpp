/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani

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

                     Implementation of class 'FastMarching3DG'

  ==============================================================================*/


#include "equations/interface/FastMarching3DG.h"
#include "linear_algebra/Vect_impl.h"


namespace OFELI {

FastMarching3DG::FastMarching3DG()
                : _nx(0), _ny(0), _nz(0)
{
}


FastMarching3DG::FastMarching3DG(const Grid&   g,
                                 Vect<real_t>& T)
{
   set(g,T);
}


FastMarching3DG::FastMarching3DG(const Grid&   g,
                                 Vect<real_t>& T,
                                 Vect<real_t>& F)
{
   set(g,T);
   _b = F;
}


FastMarching3DG::~FastMarching3DG()
{
}


void FastMarching3DG::set(const Grid&   g,
                          Vect<real_t>& T)
{
   _theGrid = &g;
   if (_theGrid->getDim()!=3)
      throw OFELIException("In FastMarching3DG::set(...): Wrong choice of space dimension");
   _nx = _theGrid->getNx(), _ny = _theGrid->getNy(), _nz = _theGrid->getNz();
   _hx = _theGrid->getHx(), _hy = _theGrid->getHy(), _hz = _theGrid->getHz();
   _u = &T;
   _b.setSize(_nx+1,_ny+1,_nz+1);
   _U.resize((_nx+1)*(_ny+1)*(_nz+1));
   _b = 1.;
}


void FastMarching3DG::set(const Grid&   g,
                          Vect<real_t>& T,
                          Vect<real_t>& F)
{
   set(g,T);
   _b = F;
}


size_t FastMarching3DG::Neigs()
{
   size_t i=_p->i, j=_p->j, k=_p->k;
   _neigs.clear();
   if (i>1 && _U[IJK(i-1,j,k)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJK(i-1,j,k)]);
   if (i<=_nx && _U[IJK(i+1,j,k)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJK(i+1,j,k)]);
   if (j>1 && _U[IJK(i,j-1,k)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJK(i,j-1,k)]);
   if (j<=_ny && _U[IJK(i,j+1,k)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJK(i,j+1,k)]);
   if (k>1 && _U[IJK(i,j,k-1)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJK(i,j,k-1)]);
   if (k<=_nz && _U[IJK(i,j,k+1)].state!=Pt::FROZEN)
      _neigs.push_back(&_U[IJK(i,j,k+1)]);
   return _neigs.size();
}


void FastMarching3DG::init()
{
   for (size_t i=1; i<=_nx+1; ++i) {
      for (size_t j=1; j<=_ny+1; ++j) {
         for (size_t k=1; k<=_nz+1; ++k) {
            _p = &_U[IJK(i,j,k)];
            _p->i = i, _p->j = j, _p->k = k;
            _p->x = _theGrid->getCoord(i,j,k).x;
            _p->y = _theGrid->getCoord(i,j,k).y;
            _p->z = _theGrid->getCoord(i,j,k).z;
            _p->sgn = Sgn((*_u)(i,j,k));
            _p->v = fabs((*_u)(i,j,k));
            _p->state = Pt::FAR;
            if (_p->v<INFINITY)
               _p->state = Pt::FROZEN;
         }
      }
   }
   for (size_t i=1; i<=_nx+1; ++i) {
      for (size_t j=1; j<=_ny+1; ++j) {
         for (size_t k=1; k<=_nz+1; ++k) {
            _p = &_U[IJK(i,j,k)];
            if (_p->state==Pt::FROZEN) {
               Neigs();
               for (auto const& N: _neigs) {
                  _np = &_U[IJK(N->i,N->j,N->k)];
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
   }
}


int FastMarching3DG::run()
{
   init();
   for (int n=0; n<_Narrow.size(); ++n) {
      _p = _Narrow[n];
      _p->state = Pt::FROZEN;
      Neigs();
      for (auto const& N: _neigs) {
         _np = &_U[IJK(N->i,N->j,N->k)];
         if (eval()==0) {
            if (_Narrow.find(_np)<0) {
               _np->state = Pt::ALIVE;
               _Narrow.insert(_np);
            }
            else
               _np->state = Pt::FROZEN;
         }
      }
   }
   for (size_t i=1; i<=_nx+1; ++i)
      for (size_t j=1; j<=_ny+1; ++j)
         for (size_t k=1; k<=_nz+1; ++k)
            (*_u)(i,j,k) = _U[IJK(i,j,k)].v*_U[IJK(i,j,k)].sgn;
   return 0;
}


real_t FastMarching3DG::eval()
{
   size_t i = _np->i, j = _np->j, k = _np->k;
   real_t vx=0., vy=0., vz=0., x=0.;
   real_t a = 1.0/(_hx*_hx), b = 1.0/(_hy*_hy), c = 1.0/(_hz*_hz);
   bool c1 = (i>1   ) && _U[IJK(i-1,j,k)].state!=Pt::FAR;
   bool c2 = (i<=_nx) && _U[IJK(i+1,j,k)].state!=Pt::FAR;
   bool d1 = (j>1   ) && _U[IJK(i,j-1,k)].state!=Pt::FAR;
   bool d2 = (j<=_ny) && _U[IJK(i,j+1,k)].state!=Pt::FAR;
   bool e1 = (k>1   ) && _U[IJK(i,j,k-1)].state!=Pt::FAR;
   bool e2 = (k<=_nz) && _U[IJK(i,j,k+1)].state!=Pt::FAR;
   if (c1 && !c2)
      vx = _U[IJK(i-1,j,k)].v;
   else if (!c1 && c2)
      vx = _U[IJK(i+1,j,k)].v;
   else if (!c1 && !c2)
      a = 0., vx = 0.;
   else if (i>1 && i<=_nx)
      vx = fmin(_U[IJK(i-1,j,k)].v,_U[IJK(i+1,j,k)].v);
   if (d1 && !d2)
      vy = _U[IJK(i,j-1,k)].v;
   else if (!d1 && d2)
      vy = _U[IJK(i,j+1,k)].v;
   else if (!d1 && !d2)
      b = 0., vy = 0.;
   else if (j>1 && j<=_ny)
      vy = fmin(_U[IJK(i,j-1,k)].v,_U[IJK(i,j+1,k)].v);
   if (e1 && !e2)
      vz = _U[IJK(i,j,k-1)].v;
   else if (!e1 && e2)
      vz = _U[IJK(i,j,k+1)].v;
   else if (!e1 && !e2)
      c = 0., vz = 0.;
   else if (k>1 && k<=_nz)
      vz = fmin(_U[IJK(i,j,k-1)].v,_U[IJK(i,j,k+1)].v);
   int ret = MaxQuad(a+b+c,-a*vx-b*vy-c*vz,a*vx*vx+b*vy*vy+c*vz*vz-1./(_b(i,j,k)*_b(i,j,k)),x);
   if (ret==0)
      _np->v = x;
   return ret;
}


real_t FastMarching3DG::getResidual()
{
   real_t ux=0., uy=0., uz=0., err=0.;
   for (size_t i=1; i<=_nx; ++i) {
      for (size_t j=1; j<=_ny; ++j) {
         for (size_t k=1; k<=_nz; ++k) {
            ux = 0.25/_hx*((*_u)(i+1,j,k)-(*_u)(i,j,k) + (*_u)(i+1,j+1,k)-(*_u)(i,j+1,k) +
                           (*_u)(i+1,j,k+1)-(*_u)(i,j,k+1) + (*_u)(i+1,j+1,k+1)-(*_u)(i,j+1,k+1));
            uy = 0.25/_hy*((*_u)(i,j+1,k)-(*_u)(i,j,k) + (*_u)(i+1,j+1,k)-(*_u)(i+1,j,k) +
                           (*_u)(i,j+1,k+1)-(*_u)(i,j,k+1) + (*_u)(i+1,j+1,k+1)-(*_u)(i+1,j,k+1));
            uz = 0.25/_hy*((*_u)(i,j,k+1)-(*_u)(i,j,k) + (*_u)(i+1,j,k+1)-(*_u)(i+1,j,k) +
                           (*_u)(i,j+1,k+1)-(*_u)(i,j+1,k) + (*_u)(i+1,j+1,k+1)-(*_u)(i+1,j+1,k));
            err += fabs(ux*ux + uy*uy + uz*uz - 1./(_b(i,j,k)*_b(i,j,k)));
         }
      }
   }
   return err/(_nx*_ny*_nz);
}

} /* namespace OFELI */
