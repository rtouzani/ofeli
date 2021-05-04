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
}


FastMarching1DG::FastMarching1DG(const Grid&   g,
                                 Vect<real_t>& T)
{
   set(g,T);
}


FastMarching1DG::FastMarching1DG(const Grid&         g,
                                 Vect<real_t>&       T,
                                 const Vect<real_t>& F)
{
   set(g,T);
   _F = F;
}


void FastMarching1DG::set(const Grid&   g,
                          Vect<real_t>& T)
{
   _theGrid = &g;
   if (_theGrid->getDim()!=1)
      throw OFELIException("In FastMarching1DG::set(...): Wrong choice of space dimension");
   _nx = _theGrid->getNx();
   _hx = _theGrid->getHx();
   _T = &T;
   _F.setSize(_nx+1);
   _u.resize(_nx+1);
   _F = 1.;
}


void FastMarching1DG::set(const Grid&         g,
                          Vect<real_t>&       T,
                          const Vect<real_t>& F)
{
   set(g,T);
   _F = F;
}


int FastMarching1DG::Neigs()
{
   int i=_p->i;
   _neigs.clear();
   if (i>1 && _u[i-1].state!=Pt::FROZEN)
      _neigs.push_back(&_u[i-1]);
   if (i<=_nx && _u[i+1].state!=Pt::FROZEN)
      _neigs.push_back(&_u[i+1]);
   return _neigs.size();
}


void FastMarching1DG::init()
{
   for (int i=1; i<=_nx+1; ++i) {
      _p = &_u[i-1];
      _p->i = i;
      _p->x = _theGrid->getCoord(i).x;
      _p->sgn = Sgn((*_T)(i));
      _p->v = fabs((*_T)(i));
      _p->state = Pt::FAR;
      if (_p->v<INFINITY)
         _p->state = Pt::FROZEN;
   }
   for (int i=1; i<=_nx+1; ++i) {
      _p = &_u[i-1];
      if (_p->state==Pt::FROZEN) {
         int ng = Neigs();
         for (int n=0; n<ng; ++n) {
            _np = &_u[_neigs[n]->i-1];
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
   for (size_t n=0; n<_Narrow.size(); ++n) {
      _p = _Narrow[n];
      if (_p->state==Pt::FROZEN)
         nb--;
      _p->state = Pt::FROZEN;
      int ng = Neigs();
      for (int k=0; k<ng; ++k) {
         _np = &_u[_neigs[k]->i-1];
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
      (*_T)(i) = _u[i-1].v*_u[i-1].sgn;
   return 0;
}


real_t FastMarching1DG::eval()
{
   int i = _np->i;
   real_t vx=0., x=0.;
   real_t a = 1.0/(_hx*_hx);
   bool c1 = (i>1   ) && _u[i-2].state!=Pt::FAR;
   bool c2 = (i<=_nx) && _u[i  ].state!=Pt::FAR;
   if (c1 && !c2)
      vx = _u[i-1].v;
   else if (!c1 && c2)
      vx = _u[i+1].v;
   else if (!c1 && !c2)
      a = 0., vx = 0.;
   else if (i>1 && i<=_nx)
      vx = fmin(_u[i-2].v,_u[i].v);
   int ret = MaxQuad(a,-a*vx,a*vx*vx-1./(_F(i)*_F(i)),x);
   if (ret==0)
      _np->v = x;
   return ret;
}


real_t FastMarching1DG::getResidual()
{
   real_t err=0., ux=0.;
   for (int i=1; i<=_nx; ++i) {
      ux = 1.0/_hx*((*_T)(i+1)-(*_T)(i));
      err += fabs(ux*ux - 1./(_F(i)*_F(i)));
   }
   return err/_nx;
}


void FastMarching1DG::ExtendSpeed(Vect<real_t>& v)
{
   for (int i=1; i<=_nx; i++)
         UpdateExt(i,v);
   for (int i=_nx; i>=1; i--)
         UpdateExt(i,v);
}


void FastMarching1DG::UpdateExt(int           i,
                                Vect<real_t>& v)
{
}

} /* namespace OFELI */
