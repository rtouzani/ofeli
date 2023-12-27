/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2024 Rachid Touzani

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

                        Implementation of class 'GeoModel'

  ==============================================================================*/

#include "solvers/GeoModel.h"
#include "linear_algebra/Vect_impl.h"
#include <cmath>

namespace OFELI {

GeoModel::GeoModel()
         : _nu(0), _cu(2), _cw(2), _npu(1), _npw(1), _cp(nullptr), _pv(nullptr)
{
}


GeoModel::GeoModel(const Vect<real_t>& b,
	           Vect<real_t>&       p)
         : _cp(&p), _pv(&b)
{
}


GeoModel::GeoModel(const Vect<real_t>& b,
                   const Vect<real_t>& h,
	           Vect<real_t>&       p)
         : _cp(&p), _h(&h), _pv(&b)
{
}


void GeoModel::setData(const Vect<real_t>& b,
                       Vect<real_t>&       p)
{
   _cp = &p;
   _pv = &b;
}


void GeoModel::setData(const Vect<real_t>& b,
                       const Vect<real_t>& h,
                       Vect<real_t>&       p)
{
   _cp = &p;
   _pv = &b;
   _h = &h;
}


void GeoModel::setBSplinePar(size_t n,
                             size_t c,
                             size_t np)
{
   _nu = n;
   _cu = c;
   _npu = np;
}


void GeoModel::setBSplineSurfacePar(size_t m,
                                    size_t n,
                                    size_t c,
                                    size_t d,
                                    size_t npu,
                                    size_t npw)
{
   _nu = n;
   _nw = m;
   _cu = c;
   _cw = d;
   _npu = npu;
   _npw = npw;
}


void GeoModel::setBezierPar(size_t n,
                            size_t nc)
{
   _cu = n;
   _npu = nc;
}


void GeoModel::setBezierSurfacePar(size_t m,
                                   size_t n,
                                   size_t npu,
                                   size_t npw)
{
   _nu = n;
   _nw = m;
   _npu = npu;
   _npw = npw;
}


void GeoModel::setNurbsPar(size_t n,
                           size_t c,
                           size_t np)
{
   _nu = n;
   _cu = c;
   _npu = np;
}


real_t GeoModel::factrl(int n)
{
   static size_t ntop=6;
   static real_t a[33]={1.0,1.0,2.0,6.0,24.0,120.0,720.0}; 
   if (n < 0)
     throw OFELIException("GeoModel::factrl(n): Negative factorial");
   if (n > 32)
     throw OFELIException("GeoModel::factrl(n): Factorial value too large");
   while (ntop < n) {
      size_t j = ntop++;
      a[ntop] = a[j]*ntop;
   }
   return a[n];
}


void GeoModel::knot(size_t       n,
                    size_t       c,
                    vector<int>& x)
{
   x[0] = 0;
   for (size_t i=1; i<n+c; ++i) {
      x[i] = x[i-1];
      if (i+1>c && i<n+1)
         x[i]++;
   }
}


void GeoModel::BSplineBasis(size_t             n,
                            size_t             c,
                            real_t             t,
                            const vector<int>& x,
                            vector<real_t>&    N)
{
// Calculate the first order basis functions
   vector<real_t> temp(35);
   for (size_t i=0; i<n+c-1; ++i) {
      temp[i] = 0.;
      if (t>=x[i] && t<x[i+1])
         temp[i] = 1.;
   }

// Calculate the higher order basis functions
   for (size_t k=2; k<=c; ++k) {
      for (size_t i=0; i<n+c-k; ++i) {
         real_t d=0., e=0.;
         if (temp[i] != 0.)
            d = (t-x[i])*temp[i]/(x[i+k-1]-x[i]);
         if (temp[i+1] != 0.)
            e = (x[i+k]-t)*temp[i+1]/(x[i+k]-x[i+1]);
         temp[i] = d + e;
      }
   }
   if (t==real_t(x[n+c-1]))
      temp[n-1] = 1.;
   for (size_t i=0; i<n; ++i)
      N[i] = temp[i];
}


real_t GeoModel::BernsteinBasis(int    n,
                                int    i,
                                real_t t)
{
   return factrl(n)/(factrl(i)*factrl(n-i))*pow(t,i)*pow(1-t,n-i);
}


void GeoModel::RationalBasis(size_t             c,
                             real_t             t,
                             size_t             n,
                             const vector<int>& x,
                             vector<real_t>&    N)
{
// Calculate the first order nonrational basis functions n[i]
   vector<real_t> temp(36);
   for (size_t i=0; i<n+c-1; ++i) {
      temp[i] = 0.;
      if (t>=x[i] && t<x[i+1])
         temp[i] = 1.;
   }

// Calculate the higher order nonrational basis functions
   for (size_t k=2; k<=c; ++k) {
      for (size_t i=0; i<n+c-k; ++i) {
//       If the lower order basis function is zero skip the calculation
         real_t d=0., e=0.;
         if (temp[i] != 0.)    
            d = (t-x[i])*temp[i]/(x[i+k-1]-x[i]);
         if (temp[i+1] != 0.)
            e = (x[i+k]-t)*temp[i+1]/(x[i+k]-x[i+1]);
         temp[i] = d + e;
      }
   }

   if (t==real_t(x[n+c-1]))
      temp[n-1] = 1.;
// Calculate sum for denominator of rational basis functions
   real_t s = 0.;
   for (size_t i=0; i<n; ++i)
      s += temp[i]*(*_h)[i];

// Form rational basis functions and put in r vector
   for (size_t i=0; i<n; ++i) {
      N[i] = 0.;
      if (s != 0.)
         N[i] = temp[i]*(*_h)[i]/s;
   }
}


void GeoModel::BSpline()
{
   vector<real_t> nbasis(20);
   for (size_t i=0; i<=_nu; ++i)
      nbasis[i] = 0.;
   vector<int> x(30);  // allows for 20 data points with basis function of order 5
   for (size_t i=0; i<=_nu+_cu; ++i)
      x[i] = 0;

// Generate the uniform open knot vector
   knot(_nu,_cu,x);

// Calculate the points on the bspline curve
   size_t icount = 0;
   real_t t = 0.;
   real_t step = real_t(x[_nu+_cu-1])/real_t(_npu-1);

   for (size_t k=0; k<_npu; ++k) {
      if (real_t(x[_nu+_cu-1])-t < NURBS_THRESHOLD)
         t = x[_nu+_cu-1];
 //   Generate the basis function for this value of t
      BSplineBasis(_nu,_cu,t,x,nbasis);
      for (size_t j=0; j<3; ++j) {    // Generate a point on the curve
         size_t jcount = j;
         (*_cp)[icount+j] = 0.;
         for (size_t i=0; i<_nu; ++i) {
            (*_cp)[icount+j] += nbasis[i]*(*_pv)[jcount];
            jcount += 3;
	 }
      }
      icount += 3, t += step;
   }
}


void GeoModel::Bezier()
{
   size_t icount=0;
   real_t t=0., step=1.0/real_t(_npu-1);
   for (size_t k=0; k<_npu; ++k) {
      if ((1.0-t) < NURBS_THRESHOLD)
         t = 1.0;
      for (size_t j=0; j<3; ++j) {
         size_t jcount = j;
         (*_cp)[icount+j] = 0.;
         for (size_t i=0; i<_nu; ++i) {
	   (*_cp)[icount+j] += BernsteinBasis(_nu-1,int(i),t)*(*_pv)[jcount];
            jcount += 3;
         }
      }
      icount += 3, t += step;
   }
}


void GeoModel::BSplineSurface()
{
   vector<int> x(30), y(30);
   for (size_t i=0; i<_nu+_cu; ++i)
      x[i] = 0;
   for (size_t i=0; i<_nw+_cw; ++i)
      y[i] = 0;
   vector<real_t> N(30), M(30);
   for (size_t i=0; i<_nu; ++i)
      N[i] = 0.;
   for (size_t i=0; i<_nw; ++i)
      M[i] = 0.;
   for (size_t i=0; i<3*_npu*_npw; ++i)
      (*_cp)[i] = 0.;

// Generate the open uniform knot vectors
   knot(_nu,_cu,x);
   knot(_nw,_cw,y);
   size_t icount = 0;

// Calculate the points on the \bsp surface
   real_t stepu = real_t(x[_nu+_cu-1])/real_t(_npu-1);
   real_t stepw = real_t(y[_nw+_cw-1])/real_t(_npw-1);
   real_t u = 0.;
   for (size_t iu=0; iu<_npu; ++iu) {
      if (real_t(x[_nu+_cu-1])-u < NURBS_THRESHOLD)
         u = x[_nu+_cu-1];
      BSplineBasis(_nu,_cu,u,x,N);
      real_t w = 0.;
      for (size_t iw=0; iw<_npw; ++iw) {
         if (real_t(y[_nw+_cw-1])-w < NURBS_THRESHOLD)
            w = y[_nw+_cw-1];
         BSplineBasis(_nw,_cw,w,y,M);
         for (size_t i=0; i<_nu; ++i) {
            if (N[i] != 0.) {
               size_t jbas = 3*_nw*i;
               for (size_t j=0; j<_nw; ++j) {
                  if (M[j] != 0.) {
                     size_t k = jbas + 3*j;
                     real_t MN = N[i]*M[j];
                     (*_cp)[icount  ] += (*_pv)[k  ]*MN;
                     (*_cp)[icount+1] += (*_pv)[k+1]*MN;
                     (*_cp)[icount+2] += (*_pv)[k+2]*MN;
                  }
               }
            }
         }
         icount += 3, w += stepw;
      }
      u += stepu;
   }
}


void GeoModel::BezierSurface()
{
   size_t icount = 0;
   real_t stepu = 1.0/real_t(_npu-1), stepw = 1.0/real_t(_npw-1);
   real_t u=0.0;

   for (size_t uinc=1; uinc<=_npu; ++uinc) {  /* for fixed u calculate various w's */
      if (1.0-u < NURBS_THRESHOLD)
//       fix up the u = 1 value because of float
         u = 1.0;
      real_t w=0.0;
      for (size_t winc=1; winc<=_npw; ++winc) {
         if (1.0-w < NURBS_THRESHOLD)
            w = 1.0;
         for (size_t i=0; i<=_nu; ++i) {
//          Bernstein basis function in the u direction (see Eq.(5.2))
	   real_t jin = BernsteinBasis(_nu,int(i),u);
            if (jin != 0.) {
               size_t jbas = 3*(_npw+1)*i;
               for (size_t j=0; j<=_npw; ++j) {
//                Bernstein basis function in the w direction (see Eq.(5.2)) 
		 real_t kjm = BernsteinBasis(_nw,int(j),w); 
                  if (kjm != 0.) {
//                   Calculate the surface points
                     size_t j1 = jbas + 3*j;	   
                     (*_cp)[icount  ] += (*_pv)[j1]*jin*kjm;
                     (*_cp)[icount+1] += (*_pv)[j1+1]*jin*kjm;
                     (*_cp)[icount+2] += (*_pv)[j1+2]*jin*kjm;
                  }
               }
            }
         }
         icount += 3, w += stepw;
      }
      u += stepu;
   }
}


void GeoModel::Nurbs()
{
   vector<real_t> N(20);
   for (size_t i=0; i<=_nu; ++i)
      N[i] = 0.;
   vector<int> x(30);
   for (size_t i=0; i<_nu+_cu; ++i)
      x[i] = 0;

   knot(_nu,_cu,x);
   size_t icount = 0;

// Calculate the points on the rational B-spline curve
   real_t t = 0;
   real_t step = real_t(x[_nu+_cu-1])/real_t(_npu-1);
   for (size_t i=0; i<_npu; ++i) {
      if (real_t(x[_nu+_cu-1])-t < NURBS_THRESHOLD)
         t = real_t(x[_nu+_cu-1]);
      RationalBasis(_cu,t,_nu,x,N);
      for (size_t j=0; j<3; ++j) {
         size_t jcount = j;
         (*_cp)[icount+j] = 0.;
         for (size_t k=0; k<_nu; ++k) {
            (*_cp)[icount+j] += N[k]*(*_pv)[jcount];
            jcount += 3;
         }
      }
      icount += 3;
      t += step;
   }
}

} /* namespace OFELI */
