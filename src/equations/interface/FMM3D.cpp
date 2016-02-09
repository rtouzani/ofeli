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

                          Implementation of class 'FMM3D'

  ==============================================================================*/

#include "equations/interface/FMM3D.h"
#include "linear_algebra/LocalVect.h"
#include "mesh/Grid.h"

namespace OFELI {

FMM3D::FMM3D(const Grid&         g,
                   Vect<real_t>* phi,
                   bool HA)
      : FMM(g,phi,HA)
{
   real_t i1, i2, j1, j2, k1, k2;
   for (size_t k=1; k<=_nz+1; ++k) {
      for (size_t i=1; i<=_nx+1; ++i) {
         for (size_t j=1; j<=_ny+1; ++j) {
            i1 = i2 = j1 = j2 = k1 = k2 = 1.;
            if (i>=2)
	      i1 = (*_phi)(i,j,k)*(*_phi)(i-1,j,k);
            if (i<=_nx)
               i2 = (*_phi)(i,j,k)*(*_phi)(i+1,j,k);
            if (j>=2)
               j1 = (*_phi)(i,j,k)*(*_phi)(i,j-1,k);
            if (j<=_ny)
               j2 = (*_phi)(i,j,k)*(*_phi)(i,j+1,k);
            if (k>=2)
               k1 = (*_phi)(i,j,k)*(*_phi)(i,j,k-1);
            if (k<=_nz)
	      k2 = (*_phi)(i,j,k)*(*_phi)(i,j,k+1);
            if ((i1<0.) || (i2<0.) || (j1<0.) || (j2<0.) || (k1<0.) || (k2<0.))
               _AlivePt(i,j,k) = (*_phi)(i,j,k);
            }
         }
      }

//    if one wants to use the method with high accuracy
//    an appropriate initial solution is necessary
      if (_high_accuracy) {
         for (size_t k=1; k<=_nz+1; ++k) {
            for (size_t i=1; i<=_nx+1; ++i) {
               for (size_t j=1; j<=_ny+1; ++j) {
                  i1 = i2 = j1 = j2 = k1 = k2 = 1.;
                  if (i>3 && _AlivePt(i-1,j,k)!=_inf)
                     i1 = (*_phi)(i,j,k)*(*_phi)(i-2,j,k);
                  if (i+2<=_nx+1 && _AlivePt(i+1,j,k)!=_inf)
                     i2 = (*_phi)(i,j,k)*(*_phi)(i+2,j,k);
                  if (j>3 && _AlivePt(i,j-1,k)!=_inf)
                     j1 = (*_phi)(i,j,k)*(*_phi)(i,j-2,k);
                  if (j+2<=_ny+1 && _AlivePt(i,j+1,k)!=_inf)
                     j2 = (*_phi)(i,j,k)*(*_phi)(i,j+2,k);
                  if (k>3 && _AlivePt(i,j,k-1)!=_inf)
                     k1 = (*_phi)(i,j,k)*(*_phi)(i,j,k-2);
                  if (k+2<=_ny+1 && _AlivePt(i,j,k+1)!=_inf)
                     k2 = (*_phi)(i,j,k)*(*_phi)(i,j,k+2);
                  if ((i1<0.) || (i2<0.) || (j1<0.) || (j2<0.) || (k1<0.) || (k2<0.))
                     _AlivePt(i,j,k) = (*_phi)(i,j,k);
               }
            }
         }
      }
}


void FMM3D::solve()
{
   IPoint pt, p;
   LocalVect<IPoint,6> Vs;
   size_t ii, jj, kk, ind;
   int sign;

// Allocate the heap
   Heap NarrowPt((_nx+1)*(_ny+1)*(_nz+1));

// This vector is an interface between the vector _AlivePt and the Heap
   _TAlive.setSize(_nx+1,_ny+1,_nz+1);

// The vector TAlive contains 3 types of points
// Alive points, Narrow points, and faraway points
   _TAlive = _AlivePt;

// Initialization of the heap
   InitHeap(NarrowPt);

// Solution
   while (NarrowPt.getSize() > 0 ) {
      pt = NarrowPt.Service();
      ii = pt.getX(); jj = pt.getY(), kk = pt.getZ();
//    add point to Alive
      _AlivePt(pt.getX(),pt.getY(),pt.getZ()) = pt.getSgn()*pt.getValue();
//    update the point as frozen
      _TAlive(pt.getX(),pt.getY(),pt.getZ()) = pt.getSgn()*pt.getValue();

//    This point and its neighbour have the same sign
      sign = Sgn(_AlivePt(pt.getX(),pt.getY(),pt.getZ()));
      pt.GenerateNeighbour(Vs);

//    Treating the point's neighbours
      for (size_t l=0; l<6; ++l) {
         ii = Vs[l].getX(); jj = Vs[l].getY(); kk = Vs[l].getZ();
         if ((ii>=1 && ii<=_nx+1) && (jj>=1 && jj<=_ny+1) && (kk>=1 && kk<=_nz+1)) {
//          faraway point
            if (_TAlive(ii,jj,kk) == _inf) {
               Evaluate(Vs[l],sign);
//             add to narrow
               NarrowPt.Add(Vs[l]);
            }

//          a Narrow point
            else if ((_TAlive(ii,jj,kk)<_inf) && (_AlivePt(ii,jj,kk)==_inf)) {
               Evaluate(Vs[l],sign);
               if (NarrowPt.Find(Vs[l],ind)) // the point is already in the heap
                  NarrowPt.Update(Vs[l].getValue(),ind);
            }
         }
      }
   }
   *_phi = _AlivePt;
}


void FMM3D::Evaluate(IPoint& pt,
                     int     sign)
{
   LocalVect<IPoint,6> Neigbour;
   IPoint pni, pni2;
   real_t val1, val2, v1, v2;
   real_t temp, aa; // constants defined for fast marching with high accuracy

// Initialization of quadratic equation coefficients
   real_t c=-_gd->getHx()*_gd->getHx(), b=0.0, a=0.0;
   GenerateDisplacement(Neigbour,true);
   size_t ii=size_t(pt.getX()), jj=size_t(pt.getY()), kk=size_t(pt.getZ());

   for (int j=0; j<3; ++j) {
      val1 = val2 = _inf;
      for (int i=0; i<2; ++i) {
         pni = pt + Neigbour[2*j+i]; // Next neighbour
         if ( (pni.getX()>=1) && (size_t(pni.getX())<=_nx+1) &&
              (pni.getY()>=1) && (size_t(pni.getY())<=_ny+1) &&
              (pni.getZ()>=1) && (size_t(pni.getZ())<=_nz+1) &&
              (std::abs(_TAlive(pni.getX(),pni.getY(),pni.getZ())) < _inf)) {
            v1 = std::abs(_TAlive(pni.getX(),pni.getY(),pni.getZ()));

            if (v1 < val1) {
               val1 = v1;
               pni2 = pt;
               pni2 = pni2 + 2*Neigbour[2*j+i];
               if ((pni2.getX()>=1) && (size_t(pni2.getX())<=_nx+1) &&
                   (pni2.getY()>=1) && (size_t(pni2.getY())<=_ny+1) &&
                   (pni2.getZ()>=1) && (size_t(pni2.getZ())<=_nz+1)) {
                  v2 = std::abs(_TAlive(pni2.getX(),pni2.getY(),pni2.getZ()));
                  if ((v2<_inf) && (v2<=v1))
                     val2 = v2;
                  else
                     val2 = _inf;
               }
            }
         }
      }

//    High Accuracy
      if ((_high_accuracy==true) && (val2!=_inf)) {
         temp = OFELI_THIRD*(4.0*val1 - val2);
         aa = 9.0/4.0;
         a += aa;
         b -= 2.0*aa*temp;
         c += aa*temp*temp;
      }
      else {
         if (val1 != _inf) {
            a += 1.0;
            b -= 2.0*val1;
            c += val1*val1;
         }
      }
   }

// Solving the quadratic equation
   real_t xmin, ymin, zmin;
   real_t max_sol = _inf;
   int ret = MaxQuadratic(a,b,c,max_sol);
   if (ret==0) {
      if (ii==1)
         xmin = std::abs(_AlivePt(ii+1,jj,kk) );
      else if (ii==_nx+1)
         xmin = std::abs(_AlivePt(ii-1,jj,kk));
      else
         xmin = std::min(std::abs(_AlivePt(ii-1,jj,kk)),std::abs(_AlivePt(ii+1,jj,kk)));

      if (jj==1)
         ymin = std::abs(_AlivePt(ii,jj+1,kk));
      else if (jj==_ny+1)
         ymin = std::abs(_AlivePt(ii,jj-1,kk));
      else
         ymin = std::min(std::abs(_AlivePt(ii,jj-1,kk)),std::abs(_AlivePt(ii,jj+1,kk)));

      if (kk==1)
         zmin = Abs(_AlivePt(ii,jj,kk+1));
      else if (jj==_nx+1)
         zmin = Abs(_AlivePt(ii,jj,kk-1) );
      else
         zmin = std::min(std::abs(_AlivePt(ii,jj,kk-1)),std::abs(_AlivePt(ii,jj,kk+1)));
      max_sol = Min(xmin,ymin,zmin) + _gd->getHx();
   }

// if we find a better distance
   if (Abs(_TAlive(ii,jj,kk)) > max_sol)
      _TAlive(ii,jj,kk) = sign*max_sol;

// updating
   pt.setValue(Abs(_TAlive(ii,jj,kk)));
   pt.setSgn(sign);
}


void FMM3D::InitHeap(Heap& NarrowPt)
{
   int sign;
   for (size_t i=1; i<_nx+1; ++i) {
      for (size_t j=1; j<_ny+1; ++j) {
         for (size_t k=1; k<_nz+1; ++k) {
            if (_AlivePt(i,j,k) < _inf) {
               IPoint pt;
               sign = Sgn(_AlivePt(i,j,k));
               if ((i>=2) && _TAlive(i-1,j,k)==_inf) {
                  pt.setX(i-1); pt.setY(j); pt.setZ(k);
                  Evaluate(pt,sign);
                  NarrowPt.Add(pt);
               }
               if ((i<=_nx) && _TAlive(i+1,j,k)==_inf) {
                  pt.setX(i+1); pt.setY(j); pt.setZ(k);
                  Evaluate(pt,sign);
                  NarrowPt.Add(pt);
               }
               if ((j>=2) && _TAlive(i,j-1,k)==_inf) {
                  pt.setX(i); pt.setY(j-1); pt.setZ(k);
                  Evaluate(pt,sign);
                  NarrowPt.Add(pt);
               }
               if ((j<=_ny) && _TAlive(i,j+1,k)==_inf) {
                  pt.setX(i); pt.setY(j+1); pt.setZ(k);
                  Evaluate(pt,sign);
                  NarrowPt.Add(pt);
               }
               if ((k>=2) && _TAlive(i,j,k-1)==_inf) {
                  pt.setX(i); pt.setY(j); pt.setZ(k-1);
                  Evaluate(pt,sign);
                  NarrowPt.Add(pt);
               }
               if ((k<=_nz) && _TAlive(i,j,k+1)==_inf) {
                  pt.setX(i); pt.setY(j); pt.setZ(k+1);
                  Evaluate(pt,sign);
                  NarrowPt.Add(pt);
               }
            }
         }
      }
   }
}


void FMM3D::ExtendSpeed(Vect<real_t>& F)
{
   real_t dx=_gd->getHx(), dy=_gd->getHy(), dz=_gd->getHz();
   size_t i, j, k=0;
   for (i=1; i<=_nx; i++) {
      for (j=1; j<=_ny; j++)
         for (k=1; k<=_nz; k++)
            UpdateExt(i,j,k,dx,dy,dz,F);
      for (j=_ny; j>=1; j--)
         UpdateExt(i,j,k,dx,dy,dz,F);
   }
   for (i=_nx; i>=1; i--) {
      for (j=1; j<=_ny; j++)
         UpdateExt(i,j,k,dx,dy,dz,F);
      for (j=_ny; j>=1; j--)
         UpdateExt(i,j,k,dx,dy,dz,F);
   }
}


void FMM3D::UpdateExt(size_t        i,
                      size_t        j,
                      size_t        k,
                      real_t        dx,
                      real_t        dy,
                      real_t        dz,
                      Vect<real_t>& F)
{
   Point<real_t> f;
   f.x = 0.25/dx*((*_phi)(i+1,j  ,k  ) + (*_phi)(i  ,j  ,k+1) + (*_phi)(i+1,j+1,k  ) +
                  (*_phi)(i+1,j+1,k+1) - (*_phi)(i  ,j  ,k  ) - (*_phi)(i  ,j  ,k+1) -
                  (*_phi)(i  ,j+1,k  ) - (*_phi)(i  ,j+1,k+1));
   f.y = 0.25/dy*((*_phi)(i  ,j+1,k  ) + (*_phi)(i  ,j+1,k+1) + (*_phi)(i+1,j+1,k  ) +
                  (*_phi)(i+1,j+1,k+1) - (*_phi)(i  ,j  ,k  ) - (*_phi)(i  ,j  ,k+1) -
                  (*_phi)(i+1,j  ,k  ) - (*_phi)(i+1,j  ,k+1));
   f.x = 0.25/dz*((*_phi)(i+1,j+1,k+1) + (*_phi)(i+1,j  ,k+1) + (*_phi)(i  ,j,k+1) +
                  (*_phi)(i  ,j+1,k+1) - (*_phi)(i  ,j  ,k  ) - (*_phi)(i  ,j+1,k  ) -
                  (*_phi)(i+1,j  ,k  ) - (*_phi)(i+1,j+1,k  ));
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
   real_t dx, dy, dz, err=0;
   real_t hx=_gd->getHx(), hy=_gd->getHy(), hz=_gd->getHz();
   for (size_t i=1; i<=_nx; i++) {
      for (size_t j=1; j<=_ny; j++) {
         for (size_t k=1; k<=_nz; k++) {
            dx = (*_phi)(i+1,j+1,k+1) + (*_phi)(i+1,j+1,k  ) + (*_phi)(i+1,j  ,k+1) +
                 (*_phi)(i+1,j  ,k  ) - (*_phi)(i  ,j+1,k+1) - (*_phi)(i  ,j+1,k  ) -
                 (*_phi)(i  ,j  ,k+1) - (*_phi)(i  ,j  ,k  );
            dy = (*_phi)(i  ,j+1,k+1) + (*_phi)(i  ,j+1,k  ) + (*_phi)(i+1,j+1,k+1) +
                 (*_phi)(i+1,j+1,k  ) - (*_phi)(i  ,j  ,k+1) - (*_phi)(i  ,j  ,k  ) -
                 (*_phi)(i+1,j  ,k+1) - (*_phi)(i+1,j  ,k  );
            dz = (*_phi)(i  ,j  ,k+1) + (*_phi)(i  ,j+1,k+1) + (*_phi)(i+1,j  ,k+1) +
                 (*_phi)(i+1,j+1,k+1) - (*_phi)(i  ,j  ,k  ) - (*_phi)(i  ,j+1,k  ) -
                 (*_phi)(i+1,j  ,k  ) - (*_phi)(i+1,j+1,k  );
            err += dx*dx/(hx*hx) + dy*dy/(hy*hy) + dz*dz/(hz*hz);
         }
      }
   }
   return fabs(0.25*sqrt(err/(_nx*_ny*_nz)) - 1);
}

} /* namespace OFELI */
