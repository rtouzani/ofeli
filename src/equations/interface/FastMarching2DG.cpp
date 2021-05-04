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


FastMarching2DG::FastMarching2DG(const Grid&         g,
                                 Vect<real_t>&       T,
                                 const Vect<real_t>& F)
{
   set(g,T,F);
}


void FastMarching2DG::set(const Grid&   g,
                          Vect<real_t>& T)
{
   _theGrid = &g;
   if (_theGrid->getDim()!=2)
      throw OFELIException("In FastMarching2DG::set(...): Wrong choice of space dimension");
   _nx = _theGrid->getNx(), _ny = _theGrid->getNy();
   _hx = _theGrid->getHx(), _hy = _theGrid->getHy();
   _T = &T;
   _F.setSize(_nx+1,_ny+1);
   _u.resize((_nx+1)*(_ny+1));
   _F = 1.;
}


void FastMarching2DG::set(const Grid&         g,
                          Vect<real_t>&       T,
                          const Vect<real_t>& F)
{
   set(g,T);
   _F = F;
}


int FastMarching2DG::Neigs()
{
   int i=_p->i, j=_p->j;
   _neigs.clear();
   if (i>1 && _u[IJ(i-1,j)].state!=Pt::FROZEN)
      _neigs.push_back(&_u[IJ(i-1,j)]);
   if (i<=_nx && _u[IJ(i+1,j)].state!=Pt::FROZEN)
      _neigs.push_back(&_u[IJ(i+1,j)]);
   if (j>1 && _u[IJ(i,j-1)].state!=Pt::FROZEN)
      _neigs.push_back(&_u[IJ(i,j-1)]);
   if (j<=_ny && _u[IJ(i,j+1)].state!=Pt::FROZEN)
      _neigs.push_back(&_u[IJ(i,j+1)]);
   return _neigs.size();
}


void FastMarching2DG::init()
{
   for (int i=1; i<=_nx+1; ++i) {
      for (int j=1; j<=_ny+1; ++j) {
         _p = &_u[IJ(i,j)];
         _p->i = i, _p->j = j;
         _p->x = _theGrid->getCoord(i,j).x;
         _p->y = _theGrid->getCoord(i,j).y;
         _p->sgn = Sgn((*_T)(i,j));
         _p->v = fabs((*_T)(i,j));
         _p->state = Pt::FAR;
         if (_p->v<INFINITY)
            _p->state = Pt::FROZEN;
      }
   }
   for (int i=1; i<=_nx+1; ++i) {
      for (int j=1; j<=_ny+1; ++j) {
         _p = &_u[IJ(i,j)];
         if (_p->state==Pt::FROZEN) {
            Neigs();
            for (auto const& N: _neigs) {
               _np = &_u[IJ(N->i,N->j)];
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


int FastMarching2DG::run()
{
   init();
   int nb = _Narrow.size();
   for (size_t n=0; n<_Narrow.size(); ++n) {
      _p = _Narrow[n];
      if (_p->state==Pt::FROZEN)
         nb--;
      _p->state = Pt::FROZEN;
      Neigs();
      for (auto const& N: _neigs) {
         _np = &_u[IJ(N->i,N->j)];
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
   for (int i=1; i<=_nx+1; ++i)
      for (int j=1; j<=_ny+1; ++j)
         (*_T)(i,j) = _u[IJ(i,j)].v*_u[IJ(i,j)].sgn;
   return 0;
}


real_t FastMarching2DG::eval()
{
   int i = _np->i, j = _np->j;
   real_t vx=0., vy=0., x=0.;
   real_t a = 1.0/(_hx*_hx), b = 1.0/(_hy*_hy);
   bool c1 = (i>1   ) && _u[IJ(i-1,j)].state!=Pt::FAR;
   bool c2 = (i<=_nx) && _u[IJ(i+1,j)].state!=Pt::FAR;
   bool d1 = (j>1   ) && _u[IJ(i,j-1)].state!=Pt::FAR;
   bool d2 = (j<=_ny) && _u[IJ(i,j+1)].state!=Pt::FAR;
   if (c1 && !c2)
      vx = _u[IJ(i-1,j)].v;
   else if (!c1 && c2)
      vx = _u[IJ(i+1,j)].v;
   else if (!c1 && !c2)
      a = 0., vx = 0.;
   else if (i>1 && i<=_nx)
      vx = fmin(_u[IJ(i-1,j)].v,_u[IJ(i+1,j)].v);
   if (d1 && !d2)
      vy = _u[IJ(i,j-1)].v;
   else if (!d1 && d2)
      vy = _u[IJ(i,j+1)].v;
   else if (!d1 && !d2)
      b = 0., vy = 0.;
   else if (j>1 && j<=_ny)
      vy = fmin(_u[IJ(i,j-1)].v,_u[IJ(i,j+1)].v);
   int ret = MaxQuad(a+b,-a*vx-b*vy,a*vx*vx+b*vy*vy-1./(_F(i,j)*_F(i,j)),x);
   if (ret==0)
      _np->v = x;
   return ret;
}


real_t FastMarching2DG::getResidual()
{
   real_t err=0., ux=0., uy=0.;
   for (int i=1; i<=_nx; ++i) {
      for (int j=1; j<=_ny; ++j) {
         ux = 0.5/_hx*((*_T)(i+1,j)-(*_T)(i,j) + (*_T)(i+1,j+1)-(*_T)(i,j+1));
         uy = 0.5/_hy*((*_T)(i,j+1)-(*_T)(i,j) + (*_T)(i+1,j+1)-(*_T)(i+1,j));
         err += fabs(ux*ux + uy*uy - 1./(_F(i,j)*_F(i,j)));
      }
   }
   return err/(_nx*_ny);
}


void FastMarching2DG::ExtendSpeed(Vect<real_t>& v)
{
   for (int i=1; i<=_nx; i++) {
      for (int j=1; j<=_ny; j++)
         UpdateExt(i,j,v);
      for (int j=_ny; j>=1; j--)
         UpdateExt(i,j,v);
   }
   for (int i=_nx; i>=1; i--) {
      for (int j=1; j<=_ny; j++)
         UpdateExt(i,j,v);
      for (int j=_ny; j>=1; j--)
         UpdateExt(i,j,v);
   }
}


void FastMarching2DG::UpdateExt(int           i,
                                int           j,
                                Vect<real_t>& v)
{
   real_t f = 0.5/_hx*((*_T)(i+1,j)+(*_T)(i+1,j+1)-(*_T)(i,j)-(*_T)(i,j+1));
   real_t g = 0.5/_hy*((*_T)(i,j+1)+(*_T)(i+1,j+1)-(*_T)(i,j)-(*_T)(i+1,j));
   real_t c=0.5*(f/_hx+g/_hy), d=0.5*(f/_hx-g/_hy);
   if (c==0.)
      v(i+1,j  ) = v(i,j+1) + c/d*(v(i,j)-v(i+1,j+1));
   else
      v(i+1,j+1) = v(i,j  ) + d/c*(v(i,j+1)-v(i+1,j));
}


} /* namespace OFELI */
