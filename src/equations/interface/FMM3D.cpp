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

                          Implementation of class 'FMM3D'

  ==============================================================================*/

#include "equations/interface/FMM3D.h"
#include "mesh/Grid.h"
#include <cmath>

namespace OFELI {

FMM3D::FMM3D(const Grid&   g,
             Vect<real_t>& phi,
             bool          HA)
      : FMM(g,phi,HA)
{
   set();
   _F = 1.;
}


FMM3D::FMM3D(const Grid&         g,
             Vect<real_t>&       phi,
             const Vect<real_t>& F,
             bool                HA)
      : FMM(g,phi,HA)
{
   set();
   _F = F;
}


void FMM3D::set()
{
   real_t i1, i2, j1, j2, k1, k2;
   for (int k=1; k<=_nz+1; ++k) {
      for (int i=1; i<=_nx+1; ++i) {
         for (int j=1; j<=_ny+1; ++j) {
            i1 = i2 = j1 = j2 = k1 = k2 = 1.;
            if (i>=2)
               i1 = (*_u)(i,j,k)*(*_u)(i-1,j,k);
            if (i<=_nx)
               i2 = (*_u)(i,j,k)*(*_u)(i+1,j,k);
            if (j>=2)
               j1 = (*_u)(i,j,k)*(*_u)(i,j-1,k);
            if (j<=_ny)
               j2 = (*_u)(i,j,k)*(*_u)(i,j+1,k);
            if (k>=2)
               k1 = (*_u)(i,j,k)*(*_u)(i,j,k-1);
            if (k<=_nz)
               k2 = (*_u)(i,j,k)*(*_u)(i,j,k+1);
            if ((i1<0.) || (i2<0.) || (j1<0.) || (j2<0.) || (k1<0.) || (k2<0.))
               _AlivePt(i,j,k) = (*_u)(i,j,k);
         }
      }
   }

// For a high accuracy version an appropriate initial solution is necessary
   if (_high_accuracy) {
      for (int k=1; k<=_nz+1; ++k) {
         for (int i=1; i<=_nx+1; ++i) {
            for (int j=1; j<=_ny+1; ++j) {
               i1 = i2 = j1 = j2 = k1 = k2 = 1.;
               if (i>3 && _AlivePt(i-1,j,k)<_inf)
                  i1 = (*_u)(i,j,k)*(*_u)(i-2,j,k);
               if (i+2<=_nx+1 && _AlivePt(i+1,j,k)<_inf)
                  i2 = (*_u)(i,j,k)*(*_u)(i+2,j,k);
               if (j>3 && _AlivePt(i,j-1,k)<_inf)
                  j1 = (*_u)(i,j,k)*(*_u)(i,j-2,k);
               if (j+2<=_ny+1 && _AlivePt(i,j+1,k)<_inf)
                  j2 = (*_u)(i,j,k)*(*_u)(i,j+2,k);
               if (k>3 && _AlivePt(i,j,k-1)<_inf)
                  k1 = (*_u)(i,j,k)*(*_u)(i,j,k-2);
               if (k+2<=_ny+1 && _AlivePt(i,j,k+1)<_inf)
                  k2 = (*_u)(i,j,k)*(*_u)(i,j,k+2);
               if ((i1<0.) || (i2<0.) || (j1<0.) || (j2<0.) || (k1<0.) || (k2<0.))
                  _AlivePt(i,j,k) = (*_u)(i,j,k);
            }
         }
      }
   }

   _F.setSize(_nx+1,_ny+1,_nz+1);
}


void FMM3D::run()
{
   IPoint pt, p;
   vector<IPoint> vs(6);
   size_t ind;

// This vector is an interface between the vector _AlivePt and the Heap
   _TAlive.setSize(_nx+1,_ny+1,_nz+1);

// The vector TAlive contains 3 types of points
// Alive points, Narrow points, and faraway points
   _TAlive = _AlivePt;

// Initialization of the heap
   init();

// Solution
   while (_Narrow.getSize()!=0) {
      pt = _Narrow.Current();
//    add point to Alive
      _AlivePt(pt.i,pt.j,pt.k) = pt.sgn*pt.val;
//    update the point as frozen
      _TAlive(pt.i,pt.j,pt.k) = pt.sgn*pt.val;

//    This point and its neighbour have the same sign
      int sign = Sgn(_AlivePt(pt.i,pt.j,pt.k));
      pt.getNeighbour(vs);

//    Treating the point's neighbours
      for (size_t l=0; l<6; ++l) {
         int ii=vs[l].i, jj=vs[l].j, kk=vs[l].k;
         if ((ii>0 && ii<=_nx+1) && (jj>0 && jj<=_ny+1) && (kk>0 && kk<=_nz+1)) {
//          faraway point
            if (_TAlive(ii,jj,kk)==_inf) {
               eval(vs[l],sign);
//             add to narrow
               _Narrow.Add(vs[l]);
            }

//          a Narrow point
            else if ((_TAlive(ii,jj,kk)<_inf) && (_AlivePt(ii,jj,kk)<_inf)) {
               eval(vs[l],sign);
               if (_Narrow.Find(vs[l],ind)) // the point is already in the heap
                  _Narrow.Update(vs[l].val,ind);
            }
         }
      }
   }
   *_u = _AlivePt;
}


void FMM3D::eval(IPoint& pt,
                 int     sign)
{
   vector<IPoint> neig(6);
   IPoint pni, pni2;
   real_t v1, v2;

// Initialization of quadratic equation coefficients
   real_t c=-_hx*_hx, b=0, a=0;
   getDisplacement(neig,true);
   int ii=pt.i, jj=pt.j, kk=pt.k;

   for (int j=0; j<3; ++j) {
      real_t val1=_inf, val2=_inf;
      for (int i=0; i<2; ++i) {
         IPoint pni = pt;
         pni += neig[2*j+i]; // Next neighbour
         if ( (pni.i>0) && (pni.i<=_nx+1) &&
              (pni.j>0) && (pni.j<=_ny+1) &&
              (pni.k>0) && (pni.k<=_nz+1) &&
              (fabs(_TAlive(pni.i,pni.j,pni.k))<_inf)) {
            v1 = fabs(_TAlive(pni.i,pni.j,pni.k));

            if (v1 < val1) {
               val1 = v1;
               IPoint pni2 = pt;
               pni2 += 2*neig[2*j+i];
               if ((pni2.i>0) && (pni2.i<=_nx+1) &&
                   (pni2.j>0) && (pni2.j<=_ny+1) &&
                   (pni2.k>0) && (pni2.k<=_nz+1)) {
                  v2 = fabs(_TAlive(pni2.i,pni2.j,pni2.k));
                  val2 = _inf;
                  if (v2<_inf && v2<=v1)
                     val2 = v2;
               }
            }
         }
      }

//    High Accuracy
      if ((_high_accuracy==true) && (val2<_inf)) {
         real_t temp = OFELI_THIRD*(4.0*val1 - val2);
         real_t aa = 2.25;
         a += aa;
         b -= aa*temp;
         c += aa*temp*temp;
      }
      else {
         a += 1.0;
         if (val1<_inf) {
            b -= val1;
            c += val1*val1;
         }
      }
   }

// Solving the quadratic equation
   real_t xmin, ymin, zmin;
   real_t max_sol = _inf;
   int ret = FMM::MaxQuad(a,b,c,max_sol);
   if (ret==0) {
      if (ii==1)
         xmin = fabs(_AlivePt(ii+1,jj,kk));
      else if (ii==_nx+1)
         xmin = fabs(_AlivePt(ii-1,jj,kk));
      else
         xmin = fmin(fabs(_AlivePt(ii-1,jj,kk)),fabs(_AlivePt(ii+1,jj,kk)));

      if (jj==1)
         ymin = fabs(_AlivePt(ii,jj+1,kk));
      else if (jj==_ny+1)
         ymin = fabs(_AlivePt(ii,jj-1,kk));
      else
         ymin = fmin(fabs(_AlivePt(ii,jj-1,kk)),fabs(_AlivePt(ii,jj+1,kk)));

      if (kk==1)
         zmin = fabs(_AlivePt(ii,jj,kk+1));
      else if (jj==_nx+1)
         zmin = fabs(_AlivePt(ii,jj,kk-1) );
      else
         zmin = fmin(fabs(_AlivePt(ii,jj,kk-1)),fabs(_AlivePt(ii,jj,kk+1)));
      max_sol = fmin(fmin(xmin,ymin),zmin) + _hx;
   }

// if we find a better distance
   if (Abs(_TAlive(ii,jj,kk)) > max_sol)
      _TAlive(ii,jj,kk) = sign*max_sol;

// updating
   pt.val = fabs(_TAlive(ii,jj,kk));
   pt.sgn = sign;
}


void FMM3D::init()
{
   int sign=0;
   for (int i=1; i<_nx+1; ++i) {
      for (int j=1; j<_ny+1; ++j) {
         for (int k=1; k<_nz+1; ++k) {
            if (_AlivePt(i,j,k)<_inf) {
               IPoint pt;
               sign = Sgn(_AlivePt(i,j,k));
               if (i>=2 && _TAlive(i-1,j,k)==_inf) {
                  pt.set(i-1,j,k);
                  eval(pt,sign);
                  _Narrow.Add(pt);
               }
               if (i<=_nx && _TAlive(i+1,j,k)==_inf) {
                  pt.set(i+1,j,k);
                  eval(pt,sign);
                  _Narrow.Add(pt);
               }
               if (j>=2 && _TAlive(i,j-1,k)==_inf) {
                  pt.set(i,j-1,k);
                  eval(pt,sign);
                  _Narrow.Add(pt);
               }
               if (j<=_ny && _TAlive(i,j+1,k)==_inf) {
                  pt.set(i,j+1,k);
                  eval(pt,sign);
                  _Narrow.Add(pt);
               }
               if (k>=2 && _TAlive(i,j,k-1)==_inf) {
                  pt.set(i,j,k-1);
                  eval(pt,sign);
                  _Narrow.Add(pt);
               }
               if (k<=_nz && _TAlive(i,j,k+1)==_inf) {
                  pt.set(i,j,k+1);
                  eval(pt,sign);
                  _Narrow.Add(pt);
               }
            }
         }
      }
   }
}


void FMM3D::ExtendSpeed(Vect<real_t>& F)
{
  /*   for (int i=1; i<=_nx; i++) {
      for (int j=1; j<=_ny; j++)
         for (int k=1; k<=_nz; k++)
            UpdateExt(i,j,k,F);
      for (int j=_ny; j>=1; j--)
         for (int k=_nz; k>=1; k--)
            UpdateExt(i,j,k,F);
   }
   for (int i=_nx; i>=1; i--) {
      for (int j=1; j<=_ny; j++)
         UpdateExt(i,j,k,F);
      for (int j=_ny; j>=1; j--)
         UpdateExt(i,j,k,F);
         }*/
}


void FMM3D::UpdateExt(int           i,
                      int           j,
                      int           k,
                      Vect<real_t>& F)
{
   Point<real_t> f;
   f.x = 0.25/_hx*((*_u)(i+1,j  ,k  ) + (*_u)(i  ,j  ,k+1) + (*_u)(i+1,j+1,k  ) +
                   (*_u)(i+1,j+1,k+1) - (*_u)(i  ,j  ,k  ) - (*_u)(i  ,j  ,k+1) -
                   (*_u)(i  ,j+1,k  ) - (*_u)(i  ,j+1,k+1));
   f.y = 0.25/_hy*((*_u)(i  ,j+1,k  ) + (*_u)(i  ,j+1,k+1) + (*_u)(i+1,j+1,k  ) +
                   (*_u)(i+1,j+1,k+1) - (*_u)(i  ,j  ,k  ) - (*_u)(i  ,j  ,k+1) -
                   (*_u)(i+1,j  ,k  ) - (*_u)(i+1,j  ,k+1));
   f.x = 0.25/_hz*((*_u)(i+1,j+1,k+1) + (*_u)(i+1,j  ,k+1) + (*_u)(i  ,j,k+1) +
                   (*_u)(i  ,j+1,k+1) - (*_u)(i  ,j  ,k  ) - (*_u)(i  ,j+1,k  ) -
                   (*_u)(i+1,j  ,k  ) - (*_u)(i+1,j+1,k  ));
   real_t a = f.x + f.y + f.z, b = f.x - f.y - f.z,
          c = f.x - f.y + f.z, d = f.x + f.y - f.z;
   if (a != 0.)
      F(i+1,j+1,k+1) = F(i,j,k) + b/a*(F(i,j+1,k+1)-F(i+1,j,k))
                                + c/a*(F(i,j+1,k)-F(i+1,j,k+1))
                                + d/a*(F(i,j,k+1)-F(i+1,j+1,k));
   else if (b != 0.)
      F(i+1,j,k) = F(i,j,k) + a/b*(F(i,j,k)-F(i+1,j+1,k+1))
                            + c/b*(F(i,j+1,k)-F(i+1,j,k+1))
                            + d/b*(F(i,j,k+1)-F(i+1,j+1,k));
   else if (c != 0.)
      F(i,j+1,k) = F(i,j,k) + a/c*(F(i+1,j+1,k+1)-F(i,j,k))
                            + b/c*(F(i+1,j,k)-F(i,j+1,k+1))
                            + d/c*(F(i+1,j+1,k)-F(i,j,k+1));
   else if (d != 0.)
      F(i,j,k+1) = F(i+1,j+1,k) + b/a*(F(i+1,j,k)-F(i+1,j+1,k+1))
                                + c/a*(F(i+1,j+1,k+1)-F(i,j,k))
                                + d/a*(F(i+1,j,k+1)-F(i,j+1,k));
   else
      ;
}


real_t FMM3D::check_error()
{
   real_t err=0;
   for (int i=1; i<=_nx; i++) {
      for (int j=1; j<=_ny; j++) {
         for (int k=1; k<=_nz; k++) {
            real_t dx = (*_u)(i+1,j+1,k+1) + (*_u)(i+1,j+1,k  ) + (*_u)(i+1,j  ,k+1) +
                        (*_u)(i+1,j  ,k  ) - (*_u)(i  ,j+1,k+1) - (*_u)(i  ,j+1,k  ) -
                        (*_u)(i  ,j  ,k+1) - (*_u)(i  ,j  ,k  );
            real_t dy = (*_u)(i  ,j+1,k+1) + (*_u)(i  ,j+1,k  ) + (*_u)(i+1,j+1,k+1) +
                        (*_u)(i+1,j+1,k  ) - (*_u)(i  ,j  ,k+1) - (*_u)(i  ,j  ,k  ) -
                        (*_u)(i+1,j  ,k+1) - (*_u)(i+1,j  ,k  );
            real_t dz = (*_u)(i  ,j  ,k+1) + (*_u)(i  ,j+1,k+1) + (*_u)(i+1,j  ,k+1) +
                        (*_u)(i+1,j+1,k+1) - (*_u)(i  ,j  ,k  ) - (*_u)(i  ,j+1,k  ) -
                        (*_u)(i+1,j  ,k  ) - (*_u)(i+1,j+1,k  );
            err += dx*dx/(_hx*_hx) + dy*dy/(_hy*_hy) + dz*dz/(_hz*_hz);
         }
      }
   }
   return fabs(0.25*sqrt(err/(_nx*_ny*_nz))-1);
}

} /* namespace OFELI */
