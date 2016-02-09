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

                        Implementation of class 'FMM2D'

  ==============================================================================*/


#include "equations/interface/FMM2D.h"
#include "linear_algebra/LocalVect.h"
#include "mesh/Grid.h"
#include <algorithm>

namespace OFELI {

FMM2D::FMM2D(const Grid&         g,
                   Vect<real_t>* phi,
                   bool          HA)
      : FMM(g,phi,HA)
{
   real_t i1, i2, j1, j2, i11, i21, j11, j21;

// Initialization 
   for (size_t i=1; i<=_nx+1; ++i) {
      for (size_t j=1; j<=_ny+1; ++j) {
         i1 = i2 = j1 = j2 = i11 = i21 = j11 = j21 = 1.;
         if (i>=2)
            i1 = (*_phi)(i,j)*(*_phi)(i-1,j);
         if (i+1<=_nx+1)
            i2 = (*_phi)(i,j)*(*_phi)(i+1,j);
            if (j>=2)
               j1 = (*_phi)(i,j)*(*_phi)(i,j-1);
            if (j<=_ny)
               j2 = (*_phi)(i,j)*(*_phi)(i,j+1);
            if ((i1<0.) || (i2<0.) || (j1<0.) || (j2<0.))
               _AlivePt(i,j) = (*_phi)(i,j);
      }
   }

// if one wants to use the method with high accuracy
// an appropriate initial solution is necessary
   if (_high_accuracy) {
     for (size_t i=1; i<=_nx+1; ++i) {
        for (size_t j=1; j<=_ny+1; ++j) {
            i1 = i2 = j1 = j2 = 1.;
            if (i>3 && _AlivePt(i-1,j) != _inf)
               i1 = (*_phi)(i,j)*(*_phi)(i-2,j);
            if (i<=_nx-1 && _AlivePt(i+1,j) != _inf)
               i2 = (*_phi)(i,j)*(*_phi)(i+2,j);
            if (j>3 && _AlivePt(i,j-1) != _inf)
               j1 = (*_phi)(i,j)*(*_phi)(i,j-2);
            if (j<=_ny-1 && _AlivePt(i,j+1)!=_inf)
               j2 = (*_phi)(i,j)*(*_phi)(i,j+2);
            if ((i1<0.) || (i2<0.) || (j1<0.) || (j2<0.))
               _AlivePt(i,j) = (*_phi)(i,j);
         }
      }
   }
}


void FMM2D::solve()
{
   IPoint pt, p;
   LocalVect<IPoint,6> Vs;
   size_t m, n, ind;
   int sign;

// Allocate the heap
   Heap NarrowPt((_nx+1)*(_ny+1));

// This vector is an interface between the vector _AlivePt and the Heap
   _TAlive.setSize((_nx+1)*(_ny+1));

// The vector TAlive contains 3 type of points
// Alive points, Narrow points,and faraway points
   _TAlive = _AlivePt;

// Initialization of the heap
   InitHeap(NarrowPt);

   while (NarrowPt.getSize() > 0) {
      pt = NarrowPt.Service(); // head of the heap
//    Add point to Alive
      _AlivePt(pt.getX(),pt.getY()) = pt.getSgn()*pt.getValue();
//    Update point as frozen
      _TAlive(pt.getX(),pt.getY()) = pt.getSgn()*pt.getValue();
	
//    This point and its neighbour have the same sign
      sign = Sgn(_AlivePt(pt.getX(),pt.getY()));

//    generating point neighbour
      pt.GenerateNeighbour(Vs);
      for (size_t k=0; k<4; ++k) {
         m = Vs[k].getX(); n = Vs[k].getY();
         if ((m>=1 && m<=_nx+1) && (n>=1 && n<=_ny+1)) {
//          faraway point
            if (_TAlive(m,n) == _inf) {
               Evaluate(Vs[k],sign);
               Vs[k].setSgn(sign);
//             add to Narrow points
               NarrowPt.Add(Vs[k]);
            }
//          a Narrow point to update
            else if ((_TAlive(m,n)<_inf) && (_AlivePt(m,n)==_inf)) {
               Evaluate(Vs[k],sign);
//             The point is already in the heap
               if (NarrowPt.Find(Vs[k],ind))
//                updating ...
                  NarrowPt.Update(Vs[k].getValue(),ind);
            }
         }
      }
   }
   *_phi = _AlivePt;
}


void FMM2D::Evaluate(IPoint& pt,
                     int     sign)
{
   LocalVect<IPoint,6> Neighbour;
   IPoint pi, pni, pni2;
   real_t val1, v1, val2, v2, temp, aa;

// Initialization of quadratic equation coefficients
   real_t c=-_gd->getHx()*_gd->getHx(), b=0., a=0.;
   GenerateDisplacement(Neighbour);
   for (size_t j=0; j<2; ++j) {
      val1 = val2 = _inf;
      for (size_t i=0; i<2; ++i) {
         pni = pt + Neighbour[2*j+i];
         if ( (pni.getX()>=1) && (size_t(pni.getX())<=_nx+1) &&  
              (pni.getY()>=1) && (size_t(pni.getY())<=_ny+1) && 
              (std::abs(_TAlive(pni.getX(),pni.getY())) < _inf)) {
            v1 = std::abs(_TAlive(pni.getX(),pni.getY()));
            if (v1 < val1) {
               val1 = v1;
               pni2 = pt;
               pni2 = pni2 + 2*Neighbour[2*j+i];
               if ((pni2.getX()>=1) && (size_t(pni2.getX())<=_nx+1) &&
                   (pni2.getY()>=1) && (size_t(pni2.getY())<=_ny+1)) {
                  v2 = std::abs(_TAlive(pni2.getX(),pni2.getY()));
                  if ((v2<_inf) && (v2<=v1))
                     val2 = v2;
                  else
                     val2 = _inf;
               }
            }
         }
      }

//    computation with high accuracy
      if (_high_accuracy && (val2!=_inf)) {
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
   real_t xmin, ymin, max_sol=_inf;
   int ret = MaxQuadratic(a,b,c,max_sol);
   if (ret==0) {
      if (pt.getX()==1)
         xmin = Abs(_AlivePt(pt.getX()+1,pt.getY()));
      else if (pt.getX() == int(_nx+1))
         xmin = Abs(_AlivePt(pt.getX()-1,pt.getY()));
      else
         xmin = std::min(Abs(_AlivePt(pt.getX()-1,pt.getY())),Abs(_AlivePt(pt.getX()+1,pt.getY())));
      if (pt.getY() == 1)
         ymin = Abs(_AlivePt(pt.getX(),pt.getY()+1));
      else if (pt.getX() == int(_nx+1))
         ymin = Abs(_AlivePt(pt.getX(),pt.getY()-1));
      else
         ymin = std::min(Abs(_AlivePt(pt.getX(),pt.getY()-1)),Abs(_AlivePt(pt.getX(),pt.getY()+1)));
      max_sol = std::min(xmin,ymin) + _gd->getHx();
   }

// if a better distance is found
   if (Abs(_TAlive(pt.getX(),pt.getY())) > max_sol)
      _TAlive(pt.getX(),pt.getY()) = sign*max_sol;

// updating
   pt.setValue(Abs(_TAlive(pt.getX(),pt.getY())));
   pt.setSgn(sign);
}


void FMM2D::InitHeap(Heap& NarrowPt)
{
   int sign;
   IPoint pt;
   for (size_t i=1; i<=_nx+1; ++i) {
      for (size_t j=1; j<=_ny+1; ++j ){
         if (_AlivePt(i,j) < _inf) {
            sign = Sgn(_AlivePt(i,j));
            if ((i-1>=1) && _TAlive(i-1,j)==_inf) {
               pt.setX(i-1); pt.setY(j);
               Evaluate(pt,sign);
               NarrowPt.Add(pt);
            }
            if ((i+1<=_nx+1) && _TAlive(i+1,j)==_inf) {
               pt.setX(i+1); pt.setY(j);
               Evaluate(pt,sign);
               NarrowPt.Add(pt);
            }
            if ((j-1>=1) && _TAlive(i,j-1)==_inf) {
               pt.setX(i); pt.setY(j-1);
               Evaluate(pt,sign);
               NarrowPt.Add(pt);
            }
            if ((j+1<=_ny+1) && _TAlive(i,j+1)==_inf) {
               pt.setX(i); pt.setY(j+1);
               Evaluate(pt,sign);
               NarrowPt.Add(pt);
            }
         }
      }
   }
}


real_t FMM2D::check_error()
{
   real_t dx, dy, err=0;
   real_t hx=_gd->getHx(), hy=_gd->getHy();
   for (size_t i=1; i<=_nx; i++) {
      for (size_t j=1; j<=_ny; j++) {
         dx = (*_phi)(i+1,j)-(*_phi)(i,j)+(*_phi)(i+1,j+1)-(*_phi)(i,j+1);
         dy = (*_phi)(i,j+1)-(*_phi)(i,j)+(*_phi)(i+1,j+1)-(*_phi)(i+1,j);
         err += dx*dx/(hx*hx) + dy*dy/(hy*hy);
      }
   }
   return fabs(0.5*sqrt(err/(_nx*_ny))-1);
}


void FMM2D::ExtendSpeed(Vect<real_t>& F)
{
   real_t dx=_gd->getHx(), dy=_gd->getHy();
   size_t i, j;
   for (i=1; i<=_nx; i++) {
      for (j=1; j<=_ny; j++)
         UpdateExt(i,j,dx,dy,F);
      for (j=_ny; j>=1; j--)
         UpdateExt(i,j,dx,dy,F);
   }
   for (i=_nx; i>=1; i--) {
      for (j=1; j<=_ny; j++)
         UpdateExt(i,j,dx,dy,F);
      for (j=_ny; j>=1; j--)
         UpdateExt(i,j,dx,dy,F);
   }
}


void FMM2D::UpdateExt(size_t        i,
                      size_t        j,
                      real_t        dx,
                      real_t        dy,
                      Vect<real_t>& F)
{
   Point<real_t> f;
   f.x = 0.5/dx*((*_phi)(i+1,j)+(*_phi)(i+1,j+1)-(*_phi)(i,j)-(*_phi)(i,j+1));
   f.y = 0.5/dy*((*_phi)(i,j+1)+(*_phi)(i+1,j+1)-(*_phi)(i,j)-(*_phi)(i+1,j));
   real_t c=0.5*(f.x/dx+f.y/dy), d=0.5*(f.x/dx-f.y/dy);
   if (c == 0.)
      F(i+1,j) = F(i,j+1) + c/d*(F(i,j)-F(i+1,j+1));
   else
      F(i+1,j+1) = F(i,j) + d/c*(F(i,j+1)-F(i+1,j));
}

} /* namespace OFELI */
