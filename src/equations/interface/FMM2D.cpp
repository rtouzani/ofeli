/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2018 Rachid Touzani

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

                        Implementation of class 'FMM2D'

  ==============================================================================*/


#include "equations/interface/FMM2D.h"
#include "mesh/Grid.h"
#include <cmath>


namespace OFELI {

FMM2D::FMM2D(const Grid&   g,
             Vect<real_t>& phi,
             bool          HA)
      : FMM(g,phi,HA)
{
   set();
}


FMM2D::FMM2D(const Grid&         g,
             Vect<real_t>&       phi,
             const Vect<real_t>& F,
             bool                HA)
      : FMM(g,phi,HA)
{
   set();
   _F = F;
}


void FMM2D::set()
{
// Initialization 
   for (int i=1; i<=_nx+1; ++i) {
      for (int j=1; j<=_ny+1; ++j) {
         real_t i1=1, i2=1, j1=1, j2=1;
         if (i>=2)
            i1 = (*_u)(i,j)*(*_u)(i-1,j);
         if (i<=_nx)
            i2 = (*_u)(i,j)*(*_u)(i+1,j);
         if (j>=2)
            j1 = (*_u)(i,j)*(*_u)(i,j-1);
         if (j<=_ny)
            j2 = (*_u)(i,j)*(*_u)(i,j+1);
         if (i1<=0. || i2<=0. || j1<=0. || j2<=0.)
            _AlivePt(i,j) = (*_u)(i,j);
      }
   }
cout<<"**\n"<<_AlivePt;

// Appropriate initial solution is necessary for high accuracy
   if (_high_accuracy) {
      for (int i=1; i<=_nx+1; ++i) {
         for (int j=1; j<=_ny+1; ++j) {
            real_t i1=1., i2=1., j1=1., j2=1.;
            if (i>3 && _AlivePt(i-1,j)!=_inf)
               i1 = (*_u)(i,j)*(*_u)(i-2,j);
            if (i<=_nx-1 && _AlivePt(i+1,j)<_inf)
               i2 = (*_u)(i,j)*(*_u)(i+2,j);
            if (j>3 && _AlivePt(i,j-1)<_inf)
               j1 = (*_u)(i,j)*(*_u)(i,j-2);
            if (j<=_ny-1 && _AlivePt(i,j+1)<_inf)
               j2 = (*_u)(i,j)*(*_u)(i,j+2);
            if (i1<0. || i2<0. || j1<0. || j2<0.)
               _AlivePt(i,j) = (*_u)(i,j);
         }
      }
   }
   _F.setSize(_nx+1,_ny+1);
   _F = 1.;
}


void FMM2D::run()
{
   IPoint p;
   size_t ind;
   vector<IPoint> vs(4);

// This vector is an interface between the vector _AlivePt and the Heap
// Vector TAlive contains 3 types of points: Alive, Narrow, and faraway
   _TAlive.setSize(_nx+1,_ny+1);
   _TAlive = _AlivePt;

// Initialization of the heap
   init();

   while (_Narrow.getSize()!=0) {
      IPoint pt = _Narrow.Current(); // head of the heap
      int ix=pt.i, iy=pt.j;
//    Add point to Alive
      _AlivePt(ix,iy) = pt.sgn*pt.val;
//    Set point as frozen
      _TAlive(ix,iy) = pt.sgn*pt.val;

//    This point and its neighbours have the same sign
      int s = Sgn(_AlivePt(ix,iy));

//    Generating neighbour point
      pt.getNeighbour(vs);
      for (size_t k=0; k<4; ++k) {
         int m=vs[k].i, n=vs[k].j;
         if ((m>0 && m<=_nx+1) && (n>0 && n<=_ny+1)) {
//          faraway point
            if (_TAlive(m,n)==_inf) {
               eval(vs[k],s);
               vs[k].sgn = s;
//             add to Narrow points
               _Narrow.Add(vs[k]);
            }
//          Narrow point to update
            else if ((_TAlive(m,n)<_inf) && (_AlivePt(m,n)==_inf)) {
               eval(vs[k],s);
//             Point is already in the heap
               if (_Narrow.Find(vs[k],ind))
//                updating ...
                  _Narrow.Update(vs[k].val,ind);
            }
         }
      }
   }
   *_u = _AlivePt;
}


void FMM2D::eval(IPoint& pt,
                 int     sign)
{
   vector<IPoint> neig(4);

// Initialization of quadratic equation coefficients
   real_t c=-_hx*_hx, b=0., a=0., aa=2.25;
   getDisplacement(neig);
   for (size_t k=0; k<2; ++k) {
      real_t val1=_inf, val2=_inf;
      for (size_t l=0; l<2; ++l) {
         IPoint pni = pt;
         pni += neig[2*k+l];
         int ix=pni.i, iy=pni.j;
         real_t v1=fabs(_TAlive(ix,iy));
         if (ix>0 && ix<=_nx+1 && iy>0 && iy<=_ny+1 && v1<_inf) {
            if (v1<val1) {
               val1 = v1;
               IPoint pni2 = pt;
               pni2 += 2*neig[2*k+l];
               int jx=pni2.i, jy=pni2.j;
               if (jx>0 && jx<=_nx+1 && jy>0 && jy<=_ny+1) {
                  real_t v2=fabs(_TAlive(jx,jy));
                  val2 = _inf;
                  if (v2<_inf && v2<=v1)
                     val2 = v2;
               }
            }
         }
      }

//    Computation with high accuracy
      if (_high_accuracy && val2<_inf) {
         real_t temp = OFELI_THIRD*(4.0*val1 - val2);
         a += aa;
         b -= aa*temp;
         c += aa*temp*temp;
      }
      else {
         if (val1<_inf) {
            a += 1.0;
            b -= val1;
            c += val1*val1;
         }
      }
   }

// Solve the quadratic equation
   real_t xmin=0, ymin=0, max_sol=_inf;
   int ix=pt.i, iy=pt.j;
   if (MaxQuad(a,b,c,max_sol)==0) {
      if (ix==1)
         xmin = fabs(_AlivePt(2,iy));
      else if (ix==_nx+1)
         xmin = fabs(_AlivePt(_nx,iy));
      else
         xmin = fmin(fabs(_AlivePt(ix-1,iy)),fabs(_AlivePt(ix+1,iy)));
      if (iy==1)
         ymin = fabs(_AlivePt(ix,2));
      else if (iy==_ny+1)
         ymin = fabs(_AlivePt(ix,_ny));
      else
         ymin = fmin(fabs(_AlivePt(ix,iy-1)),fabs(_AlivePt(ix,iy+1)));
      max_sol = fmin(xmin,ymin) + _hx;
   }

// if a better distance is found
   if (fabs(_TAlive(ix,iy)) > max_sol)
      _TAlive(ix,iy) = sign*max_sol;

// updating
   pt.val = fabs(_TAlive(ix,iy));
   pt.sgn = sign;
}


void FMM2D::init()
{
   _added.setSize(_nx+1,_ny+1);
   for (int i=1; i<=_nx+1; ++i) {
      for (int j=1; j<=_ny+1; ++j) {
         if (_AlivePt(i,j)<_inf) {
            int s = Sgn(_AlivePt(i,j));
            if (i>1 && _TAlive(i-1,j)==_inf)
               add(i-1,j,s);
            if (i<=_nx && _TAlive(i+1,j)==_inf)
               add(i+1,j,s);
            if (j>1 && _TAlive(i,j-1)==_inf)
               add(i,j-1,s);
            if (j<=_ny && _TAlive(i,j+1)==_inf)
               add(i,j+1,s);
         }
      }
   }
}


void FMM2D::add(int i,
                int j,
                int s)
{
   if (_added(i,j)==0) {
      IPoint pt(i,j);
      eval(pt,s);
      _Narrow.Add(pt), _added(i,j) = 1;
   }
}


real_t FMM2D::check_error()
{
   real_t err=0;
   for (int i=1; i<=_nx; i++) {
      for (int j=1; j<=_ny; j++) {
         real_t dx = (*_u)(i+1,j)-(*_u)(i,j)+(*_u)(i+1,j+1)-(*_u)(i,j+1),
                dy = (*_u)(i,j+1)-(*_u)(i,j)+(*_u)(i+1,j+1)-(*_u)(i+1,j);
         err += dx*dx/(_hx*_hx) + dy*dy/(_hy*_hy);
      }
   }
   return fabs(0.5*sqrt(err/(_nx*_ny))-1);
}


void FMM2D::ExtendSpeed(Vect<real_t>& F)
{
   for (int i=1; i<=_nx; i++) {
      for (int j=1; j<=_ny; j++)
         UpdateExt(i,j,F);
      for (int j=_ny; j>=1; j--)
         UpdateExt(i,j,F);
   }
   for (int i=_nx; i>=1; i--) {
      for (int j=1; j<=_ny; j++)
         UpdateExt(i,j,F);
      for (int j=_ny; j>=1; j--)
         UpdateExt(i,j,F);
   }
}


void FMM2D::UpdateExt(int           i,
                      int           j,
                      Vect<real_t>& F)
{
   real_t f = 0.5/_hx*((*_u)(i+1,j)+(*_u)(i+1,j+1)-(*_u)(i,j)-(*_u)(i,j+1)),
          g = 0.5/_hy*((*_u)(i,j+1)+(*_u)(i+1,j+1)-(*_u)(i,j)-(*_u)(i+1,j));
   real_t c=0.5*(f/_hx+g/_hy), d=0.5*(f/_hx-g/_hy);
   if (c==0.)
      F(i+1,j  ) = F(i,j+1) + c/d*(F(i,j)-F(i+1,j+1));
   else
      F(i+1,j+1) = F(i,j  ) + d/c*(F(i,j+1)-F(i+1,j));
}

} /* namespace OFELI */
