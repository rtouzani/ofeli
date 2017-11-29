/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

FMM2D::FMM2D(const Grid&   g,
             Vect<real_t>& phi,
             bool          HA)
      : FMM(g,phi,HA)
{
   real_t i1, i2, j1, j2, i11, i21, j11, j21;

// Initialization 
cout<<"PHI\n"<<*_phi<<endl;
   for (size_t i=1; i<=_nx+1; ++i) {
      for (size_t j=1; j<=_ny+1; ++j) {
         i1 = i2 = j1 = j2 = i11 = i21 = j11 = j21 = 1.;
         if (i>=2)
            i1 = (*_phi)(i,j)*(*_phi)(i-1,j);
         if (i<=_nx)
            i2 = (*_phi)(i,j)*(*_phi)(i+1,j);
         if (j>=2)
            j1 = (*_phi)(i,j)*(*_phi)(i,j-1);
         if (j<=_ny)
            j2 = (*_phi)(i,j)*(*_phi)(i,j+1);
         if (i1<=0. || i2<=0. || j1<=0. || j2<=0.)
            _AlivePt(i,j) = (*_phi)(i,j);
      }
   }
cout<<"ALIVE_PT\n"<<_AlivePt<<endl;

// If one wants to use the method with high accuracy
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
            if (i1<0. || i2<0. || j1<0. || j2<0.)
               _AlivePt(i,j) = (*_phi)(i,j);
         }
      }
   }
}


void FMM2D::solve()
{
   IPoint p;
   LocalVect<IPoint,6> Vs;
   size_t ind;

// Allocate the heap
   Heap NarrowPt((_nx+1)*(_ny+1));

// This vector is an interface between the vector _AlivePt and the Heap
   _TAlive.setSize(_nx+1,_ny+1);

// Vctor TAlive contains 3 types of points: Alive, Narrow, and faraway
   _TAlive = _AlivePt;

// Initialization of the heap
cout<<"TAlive 1\n"<<_TAlive;
   InitHeap(NarrowPt);
cout<<"TAlive 2\n"<<_TAlive;

   while (NarrowPt.getSize() > 0) {
cout<<"\n\nNb. of narrow points: "<<NarrowPt.getSize()<<endl;
      IPoint pt = NarrowPt.Current(); // head of the heap
      size_t ix=pt.getX(), iy=pt.getY();
cout<<"POINT  ("<<ix<<","<<iy<<")"<<endl;
//    Add point to Alive
      _AlivePt(ix,iy) = pt.getSgn()*pt.getValue();
//    Update point as frozen
      _TAlive(ix,iy) = pt.getSgn()*pt.getValue();
cout<<"Point is frozen"<<endl;

//    This point and its neighbour have the same sign
      int s = Sgn(_AlivePt(ix,iy));

//    generating neighbour point
      pt.GenerateNeighbour(Vs);
      for (size_t k=0; k<4; ++k) {
         size_t m=Vs[k].getX(), n=Vs[k].getY();
cout<<"Neighbour: "<<k<<"  mn = "<<m<<" "<<n<<endl;
         if ((m>=1 && m<=_nx+1) && (n>=1 && n<=_ny+1)) {
//          faraway point
            if (_TAlive(m,n) == _inf) {
cout<<"Faraway point: "<<m<<"  "<<n<<endl;
               Evaluate(Vs[k],s);
               Vs[k].setSgn(s);
//             add to Narrow points
               NarrowPt.Add(Vs[k]);
            }
//          Narrow point to update
            else if ((_TAlive(m,n)<_inf) && (_AlivePt(m,n)==_inf)) {
cout<<"Narrow point to update: "<<m<<"  "<<n<<endl;
               Evaluate(Vs[k],s);
//             Point is already in the heap
               if (NarrowPt.Find(Vs[k],ind)) {
//                updating ...
cout<<"Already in heap: "<<m<<"  "<<n<<" ... update"<<endl;
                  NarrowPt.Update(Vs[k].getValue(),ind);
               }
            }
         }
      }
   }
   *_phi = _AlivePt;
cout<<*_phi;
}


void FMM2D::Evaluate(IPoint& pt,
                     int     sign)
{
   LocalVect<IPoint,6> neig;
   real_t val1, v1, val2, v2;

// Initialization of quadratic equation coefficients
   real_t c=-_gd->getHx()*_gd->getHx(), b=0., a=0., aa=9.0/4.0;
   GenerateDisplacement(neig);
   for (size_t j=0; j<2; ++j) {
      val1 = val2 = _inf;
      for (size_t i=0; i<2; ++i) {
         IPoint pni = pt + neig[2*j+i];
         size_t ix=pni.getX(), iy=pni.getY();
         if (ix>=1 && ix<=_nx+1 && iy>=1 && iy<=_ny+1 && 
             std::abs(_TAlive(ix,iy)) < _inf) {
            v1 = std::abs(_TAlive(ix,iy));
            if (v1 < val1) {
               val1 = v1;
               IPoint pni2 = pt + 2*neig[2*j+i];
               size_t jx=pni2.getX(), jy=pni2.getY();
               if (jx>=1 && jx<=_nx+1 && jy>=1 && jy<=_ny+1) {
                  v2 = std::abs(_TAlive(jx,jy));
                  if (v2<_inf && v2<=v1)
                     val2 = v2;
                  else
                     val2 = _inf;
               }
            }
         }
      }

//    Computation with high accuracy
      if (_high_accuracy && val2!=_inf) {
         real_t temp = (4.0*val1 - val2)/3.0;
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
   real_t xmin=0, ymin=0, max_sol=_inf;
cout<<"** max_sol: "<<max_sol<<" *** "<<a<<" "<<b<<" "<<c<<endl;
   size_t ix=pt.getX(), iy=pt.getY();
   if (FMM::MaxQuadratic(a,b,c,max_sol)==0) {
      if (ix==1)
         xmin = std::abs(_AlivePt(2,iy));
      else if (ix==_nx+1)
         xmin = std::abs(_AlivePt(_nx,iy));
      else
         xmin = std::min(std::abs(_AlivePt(ix-1,iy)),std::abs(_AlivePt(ix+1,iy)));
      if (iy==1)
         ymin = std::abs(_AlivePt(ix,2));
      //      else if (ix==_nx+1)
      else if (iy==_ny+1)
         ymin = std::abs(_AlivePt(ix,_ny));
      else
         ymin = std::min(std::abs(_AlivePt(ix,iy-1)),std::abs(_AlivePt(ix,iy+1)));
      max_sol = std::min(xmin,ymin) + _gd->getHx();
cout<<"***"<<endl;
   }
cout<<"++ max_sol: "<<max_sol<<endl;

// if a better distance is found
   if (Abs(_TAlive(ix,iy)) > max_sol)
      _TAlive(ix,iy) = sign*max_sol;

// updating
   pt.setValue(Abs(_TAlive(ix,iy)));
cout<<"UPDATED: "<<ix<<" "<<iy<<": "<<pt.getValue()<<endl;
   pt.setSgn(sign);
}


void FMM2D::InitHeap(Heap& NarrowPt)
{
   IPoint pt;
   Vect<int> added(_nx+1,_ny+1);
cout<<"in InitHeap\n"<<_AlivePt<<endl;

   for (int i=1; i<=int(_nx)+1; ++i) {
      for (int j=1; j<=int(_ny)+1; ++j ) {
cout<<"INIT: "<<i<<" "<<j<<endl;
         if (_AlivePt(i,j)<_inf) {
cout<<"OK"<<endl;
            int s = Sgn(_AlivePt(i,j));
            if (i>=2 && _TAlive(i-1,j)==_inf) {
cout<<"Processing ("<<i-1<<","<<j<<")"<<endl;
               pt.setXY(i-1,j);
               Evaluate(pt,s);
               if (added(i-1,j)==0)
                  NarrowPt.Add(pt), added(i-1,j) = 1;
               NarrowPt.Add(pt);
            }
            if (i<=int(_nx) && _TAlive(i+1,j)==_inf) {
cout<<"Processing ("<<i+1<<","<<j<<")"<<endl;
               pt.setXY(i+1,j);
               Evaluate(pt,s);
               if (added(i+1,j)==0)
                  NarrowPt.Add(pt), added(i+1,j) = 1;
               NarrowPt.Add(pt), added(i-1,j) = 1;
            }
            if (j>=2 && _TAlive(i,j-1)==_inf) {
cout<<"Processing ("<<i<<","<<j-1<<")"<<endl;
               pt.setXY(i,j-1);
               Evaluate(pt,s);
               if (added(i,j-1)==0)
                  NarrowPt.Add(pt), added(i,j-1) = 1;
                NarrowPt.Add(pt), added(i-1,j) = 1;
            }
            if (j<=int(_ny) && _TAlive(i,j+1)==_inf) {
cout<<"Processing ("<<i<<","<<j+1<<")"<<endl;
               pt.setXY(i,j+1);
               Evaluate(pt,s);
               if (added(i,j+1)==0)
                  NarrowPt.Add(pt), added(i,j+1) = 1;
               NarrowPt.Add(pt), added(i-1,j) = 1;
            }
         }
      }
   }
cout<<"NARROW\n"<<NarrowPt<<endl;
}


real_t FMM2D::check_error()
{
   real_t err=0, hx=_gd->getHx(), hy=_gd->getHy();
   for (int i=1; i<=_nx; i++) {
      for (int j=1; j<=_ny; j++) {
         real_t dx = (*_phi)(i+1,j)-(*_phi)(i,j)+(*_phi)(i+1,j+1)-(*_phi)(i,j+1),
                dy = (*_phi)(i,j+1)-(*_phi)(i,j)+(*_phi)(i+1,j+1)-(*_phi)(i+1,j);
         err += dx*dx/(hx*hx) + dy*dy/(hy*hy);
      }
   }
   return fabs(0.5*sqrt(err/(_nx*_ny))-1);
}


void FMM2D::ExtendSpeed(Vect<real_t>& F)
{
   real_t dx=_gd->getHx(), dy=_gd->getHy();
   for (int i=1; i<=_nx; i++) {
      for (int j=1; j<=_ny; j++)
         UpdateExt(i,j,dx,dy,F);
      for (int j=_ny; j>=1; j--)
         UpdateExt(i,j,dx,dy,F);
   }
   for (int i=_nx; i>=1; i--) {
      for (int j=1; j<=_ny; j++)
         UpdateExt(i,j,dx,dy,F);
      for (int j=_ny; j>=1; j--)
         UpdateExt(i,j,dx,dy,F);
   }
}


void FMM2D::UpdateExt(int           i,
                      int           j,
                      real_t        dx,
                      real_t        dy,
                      Vect<real_t>& F)
{
   real_t f = 0.5/dx*((*_phi)(i+1,j)+(*_phi)(i+1,j+1)-(*_phi)(i,j)-(*_phi)(i,j+1)),
          g = 0.5/dy*((*_phi)(i,j+1)+(*_phi)(i+1,j+1)-(*_phi)(i,j)-(*_phi)(i+1,j));
   real_t c=0.5*(f/dx+g/dy), d=0.5*(f/dx-g/dy);
   if (c == 0.)
      F(i+1,j) = F(i,j+1) + c/d*(F(i,j)-F(i+1,j+1));
   else
      F(i+1,j+1) = F(i,j) + d/c*(F(i,j+1)-F(i+1,j));
}

} /* namespace OFELI */
