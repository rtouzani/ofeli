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

                        Implementation of class 'FuncApprox'

  ==============================================================================*/

#include "solvers/FuncApprox.h"
#include "linear_algebra/Vect_impl.h"
#include "linear_algebra/DMatrix_impl.h"
#include "util/util.h"
#include <cmath>

namespace OFELI {

FuncApprox::FuncApprox()
           : _alloc(false), _n1(0), _n2(1), _c1(2), _c2(2), _np1(1), _np2(1),
             _p(nullptr), _q(nullptr)
{
   _ls_opt = -1;
}


FuncApprox::~FuncApprox()
{
   if (_alloc)
      delete _B;
}


void FuncApprox::setLagrange(int                 n,
                             const Vect<real_t>& x,
                             const Vect<real_t>& y,
                             Fct&                f)
{
   _ap = LAGRANGE;
   if (n<=0)
      throw OFELIException("FuncApprox::setLagrange(..): Illegal degree of interpolation polynomial.");
   _x = &x;
   _y = &y;
   _ffct = &f;
   _ffct->setVar("x");
   _degree = n;
}


void FuncApprox::getLeastSquare(Fct& f)
{
   if (_ls_opt<0)
      throw OFELIException("FuncApprox::getLeastSquare(Fct&): No least square approximation defined.");
   string p = "", v="x";
   f.setVar("x");
   if (_ls_opt==0) {
      for (size_t i=0; i<_n2; ++i) {
         Fct &ff = *(*_fct)[i];
         ff.setVar("x");
         if (ff.getVar(1) != v)
            throw OFELIException("FuncApprox::getLeastSquare(Fct&): Inconsistency in variable name.");
         p = p + MonomialExpression((*_p)[i],(*_fct)[i]->getExpression(),i);
      }
      f.set(p);
   }
   else if (_ls_opt==2) {
      p = toString((*_p)[0]) + sstring((*_p)[1]) + "*" + v;
      for (size_t i=2; i<_n2; ++i)
         p = p + sstring((*_p)[i]) + "*" + v + "^" + toString(i);
   }
   else if (_ls_opt==3)
      p = toString((*_p)[0]) + sstring((*_p)[1]) + "*" + v;
   else
      throw OFELIException("FuncApprox::getLeastSquare(Fct&): Function expression unavailable for this choice.");
   f.set(p);
}


void FuncApprox::setLeastSquare(const vector<Fct *>& f,
                                const Vect<real_t>&  x,
                                const Vect<real_t>&  y,
                                Vect<real_t>&        a)
{
   _ap = LEAST_SQUARE;
   _ls_opt = 0;
   _fct = &f;
   _x = &x;
   _y = &y;
   _p = &a;
   _n1 = _x->size();
   _n2 = _p->size();
   _A.setSize(_n2);
   _b.setSize(_n2);
   _B = new DMatrix<real_t>(_n1,_n2);
   _alloc = 1;
   for (size_t i=0; i<_n2; ++i) {
      Fct &ff = *(*_fct)[i];
      ff.setVar("x");
      if (ff.getNbVar()>1)
         throw OFELIException("FuncApprox::setLeastSquare(..): Functions must have one variable only.");
      string x = ff.getVar(1);
      for (size_t j=0; j<_n1; ++j)
         (*_B)(j+1,i+1) = ff((*_x)[j]);
   }
   for (size_t i=1; i<=_n2; ++i) {
      real_t s = 0.;
      for (size_t k=1; k<=_n1; ++k)
         s += (*_B)(i,k)*(*_y)(i);
      _b(i) = s;
      for (size_t j=1; j<=_n2; ++j) {
         s = 0.;
         for (size_t k=1; k<=_n1; ++k)
            s += (*_B)(k,i)*(*_B)(k,j);
         _A(i,j) = s;
      }
   }
}


void FuncApprox::setLeastSquare(DMatrix<real_t>&    B,
                                const Vect<real_t>& y,
                                Vect<real_t>&       a)
{
   _ap = LEAST_SQUARE;
   _ls_opt = 1;
   _B = &B;
   _y = &y;
   _p = &a;
   _n1 = _x->size();
   _n2 = _p->size();
   _A.setSize(_n2);
   _b.setSize(_n2);
   _alloc = false;
   for (size_t i=1; i<=_n2; ++i) {
      real_t s = 0.;
      for (size_t k=1; k<=_n1; ++k)
         s += (*_B)(k,i)*(*_y)(i);
      _b(i) = s;
      for (size_t j=1; j<=_n2; ++j) {
         s = 0.;
         for (size_t k=1; k<=_n1; ++k)
            s += (*_B)(k,i)*(*_B)(k,j);
         _A(i,j) = s;
      }
   }
}


void FuncApprox::setLeastSquare(const Vect<real_t>& x,
                                const Vect<real_t>& y,
                                size_t              N,
                                Vect<real_t>&       a)
{
   _ap = LEAST_SQUARE;
   _ls_opt = 2;
   _x = &x;
   _y = &y;
   _p = &a;
   _n1 = _x->size();
   _n2 = N + 1;
   _A.setSize(_n2);
   _b.setSize(_n2);
   _B = new DMatrix<real_t>(_n1,_n2);
   _alloc = true;
   if (_p->size()<_n2)
      throw OFELIException("FuncApprox::setLeastSquare(..): Incorrect size of vector a.");
   for (size_t i=1; i<=_n1; ++i) {
      (*_B)(i,1) = 1.;
      for (size_t j=2; j<=_n2; ++j)
         (*_B)(i,j) = (*_B)(i,j-1)*(*_x)(i);
   }
   for (size_t i=1; i<=_n2; ++i) {
      real_t s = 0.;
      for (size_t k=1; k<=_n1; ++k)
         s += (*_B)(k,i)*(*_y)(i);
      _b(i) = s;
      for (size_t j=1; j<=_n2; ++j) {
         s = 0.;
         for (size_t k=1; k<=_n1; ++k)
            s += (*_B)(k,i)*(*_B)(k,j);
         _A(i,j) = s;
      }
   }
}


void FuncApprox::setLeastSquare(const Vect<real_t>& x,
                                const Vect<real_t>& y,
                                real_t&             a0,
                                real_t&             a1)
{
   _ap = LEAST_SQUARE;
   _ls_opt = 3;
   _x = &x;
   _y = &y;
   _a0 = &a0;
   _a1 = &a1;
   _n1 = _x->size();
   _n2 = 2;
   _A.setSize(2);
   _b.setSize(2);
   _alloc = false;
}

  
void FuncApprox::setBSpline(size_t              n,
                            size_t              c,
                            size_t              np,
                            const Vect<real_t>& b,
                            Vect<real_t>&       p)
{
   _ap = BSPLINE;
   _n1 = n;
   _c1 = c;
   _np1 = np;
   _p = &p;
   _q = &b;
}


void FuncApprox::setBSplineSurface(size_t              m,
                                   size_t              n,
                                   size_t              c,
                                   size_t              d,
                                   size_t              npu,
                                   size_t              npw,
                                   const Vect<real_t>& b,
                                   Vect<real_t>&       p)
{
   _ap = BSPLINE_SURFACE;
   _n1 = n, _n2 = m;
   _c1 = c, _c2 = d;
   _np1 = npu, _np2 = npw;
   _p = &p;
   _q = &b;
}


void FuncApprox::setBezier(size_t              n,
                           size_t              nc,
                           const Vect<real_t>& b,
                           Vect<real_t>&       p)
{
   _ap = BEZIER;
   _n1 = n;
   _np1 = nc;
   _p = &p;
   _q = &b;
}


void FuncApprox::setBezierSurface(size_t              m,
                                  size_t              n,
                                  size_t              npu,
                                  size_t              npw,
                                  const Vect<real_t>& b,
                                  Vect<real_t>&       p)
{
   _ap = BEZIER_SURFACE;
   _n1 = n, _n2 = m;
   _np1 = npu, _np2 = npw;
   _p = &p;
   _q = &b;
}


void FuncApprox::setNurbs(size_t              n,
                          size_t              c,
                          size_t              np,
                          const Vect<real_t>& b,
                          const Vect<real_t>& h,
                          Vect<real_t>&       p)
{
   _ap = NURBS;
   _n1 = n;
   _c1 = c;
   _np1 = np;
   _p = &p;
   _h = &h;
   _q = &b;
}


void FuncApprox::setNurbsSurface(size_t              m,
                                 size_t              n,
                                 size_t              c,
                                 size_t              d,
                                 size_t              npu,
                                 size_t              npw,
                                 const Vect<real_t>& b,
                                 Vect<real_t>&       p)
{
   _ap = NURBS_SURFACE;
   _n1 = n, _n2 = m;
   _c1 = c, _c2 = d;
   _np1 = npu, _np2 = npw;
   _p = &p;
   _q = &b;
   _qq.resize(b.size());
   _qq = b;
}


real_t FuncApprox::factrl(int n)
{
   static size_t ntop=6;
   static real_t a[33]={1.0,1.0,2.0,6.0,24.0,120.0,720.0}; 
   if (n < 0)
     throw OFELIException("FuncApprox::factrl(n): Negative factorial");
   if (n > 32)
     throw OFELIException("FuncApprox::factrl(n): Factorial value too large");
   while (int(ntop) < n) {
      size_t j = ntop++;
      a[ntop] = a[j]*ntop;
   }
   return a[n];
}


void FuncApprox::knot(size_t       n,
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


real_t FuncApprox::sumrbas(const vector<real_t>& N,
                           const vector<real_t>& M)
{
   real_t s = 0.0;
   for (size_t i=0; i<_n1; ++i) {
      for (size_t j=0; j<_n2; ++j)
         s += _qq[4*_n2*i+4*j+3]*N[i]*M[j];
   }
   return s;
}


void FuncApprox::BSplineBasis(size_t             n,
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


real_t FuncApprox::BernsteinBasis(int    n,
                                  int    i,
                                  real_t t)
{
   return factrl(n)/(factrl(i)*factrl(n-i))*pow(t,i)*pow(1-t,n-i);
}


void FuncApprox::RationalBasis(size_t             c,
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
      for (size_t i=1; i<=n+c-k; ++i) {
//       If the lower order basis function is zero skip the calculation
         real_t d=0., e=0.;
         if (temp[i-1] != 0.)    
            d = (t-x[i-1])*temp[i-1]/(x[i+k-2]-x[i-1]);
         if (temp[i] != 0.)
            e = (x[i+k-1]-t)*temp[i]/(x[i+k-1]-x[i]);
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


void FuncApprox::SRationalBasis(size_t             c,
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


int FuncApprox::runBSpline()
{
   vector<real_t> N(20);
   clear(N);
   vector<int> x(30);  // allows for 20 data points with basis function of order 5
   clear(x);

// Generate the uniform open knot vector
   knot(_n1,_c1,x);

// Calculate the points on the bspline curve
   size_t icount = 0;
   real_t t = 0.;
   real_t step = real_t(x[_n1+_c1-1])/real_t(_np1-1);

   for (size_t k=0; k<_np1; ++k) {
      if (real_t(x[_n1+_c1-1])-t < FUNC_APPROX_THRESHOLD)
         t = x[_n1+_c1-1];
 //   Generate the basis function for this value of t
      BSplineBasis(_n1,_c1,t,x,N);
      for (size_t j=0; j<3; ++j) {    // Generate a point on the curve
         size_t jcount = j;
         (*_p)[icount+j] = 0.;
         for (size_t i=0; i<_n1; ++i) {
            (*_p)[icount+j] += N[i]*(*_q)[jcount];
            jcount += 3;
	 }
      }
      icount += 3, t += step;
   }
   return 0;
}


int FuncApprox::runBezier()
{
   size_t icount=0;
   real_t t=0., step=1.0/real_t(_np1-1);
   for (size_t k=0; k<_np1; ++k) {
      if ((1.0-t) < FUNC_APPROX_THRESHOLD)
         t = 1.0;
      for (size_t j=0; j<3; ++j) {
         size_t jcount = j;
         (*_p)[icount+j] = 0.;
         for (size_t i=0; i<_n1; ++i) {
            (*_p)[icount+j] += BernsteinBasis(_c1-1,int(i),t)*(*_q)[jcount];
            jcount += 3;
         }
      }
      icount += 3, t += step;
   }
   return 0;
}


int FuncApprox::runBSplineSurface()
{
   vector<int> x(30), y(30);
   clear(x);
   clear(y);
   vector<real_t> N(30), M(30);
   clear(N);
   clear(M);
   clear(*_p);

// Generate the open uniform knot vectors
   knot(_n1,_c1,x);
   knot(_n2,_c2,y);
   size_t icount = 0;

// Calculate the points on the \bsp surface
   real_t stepu = real_t(x[_n1+_c1-1])/real_t(_np1-1);
   real_t stepw = real_t(y[_n2+_c2-1])/real_t(_np2-1);
   real_t u = 0.;
   for (size_t i1=0; i1<_np1; ++i1) {
      if (real_t(x[_n1+_c1-1])-u < FUNC_APPROX_THRESHOLD)
         u = x[_n1+_c1-1];
      BSplineBasis(_n1,_c1,u,x,N);
      real_t w = 0.;
      for (size_t i2=0; i2<_np2; ++i2) {
         if (real_t(y[_n2+_c2-1])-w < FUNC_APPROX_THRESHOLD)
            w = y[_n2+_c2-1];
         BSplineBasis(_n2,_c2,w,y,M);
         for (size_t i=0; i<_n1; ++i) {
            if (N[i] != 0.) {
               size_t jbas = 3*_n2*i;
               for (size_t j=0; j<_n2; ++j) {
                  if (M[j] != 0.) {
                     size_t k = jbas + 3*j;
                     real_t MN = N[i]*M[j];
                     (*_p)[icount  ] += (*_q)[k  ]*MN;
                     (*_p)[icount+1] += (*_q)[k+1]*MN;
                     (*_p)[icount+2] += (*_q)[k+2]*MN;
                  }
               }
            }
         }
         icount += 3, w += stepw;
      }
      u += stepu;
   }
   return 0;
}


int FuncApprox::runBezierSurface()
{
   size_t icount = 0;
   real_t stepu = 1.0/real_t(_np1-1), stepw = 1.0/real_t(_np2-1);
   real_t u=0.0;

   for (size_t i1=1; i1<=_np1; ++i1) {
      if (1.0-u < FUNC_APPROX_THRESHOLD)
//       fix up the u = 1 value because of float
         u = 1.0;
      real_t w=0.0;
      for (size_t i2=0; i2<_np2; ++i2) {
         if (1.0-w < FUNC_APPROX_THRESHOLD)
            w = 1.0;
         for (size_t i=0; i<=_n1; ++i) {
//          Bernstein basis function in the u direction (see Eq.(5.2))
            real_t jin = BernsteinBasis(_n1,int(i),u);
            if (jin != 0.) {
               size_t jbas = 3*(_np2+1)*i;
               for (size_t j=0; j<=_np2; ++j) {
//                Bernstein basis function in the w direction (see Eq.(5.2)) 
                  real_t kjm = BernsteinBasis(_n2,int(j),w); 
                  if (kjm != 0.) {
//                   Calculate the surface points
                     size_t j1 = jbas + 3*j;	   
                     (*_p)[icount  ] += (*_q)[j1  ]*jin*kjm;
                     (*_p)[icount+1] += (*_q)[j1+1]*jin*kjm;
                     (*_p)[icount+2] += (*_q)[j1+2]*jin*kjm;
                  }
               }
            }
         }
         icount += 3, w += stepw;
      }
      u += stepu;
   }
   return 0;
}


int FuncApprox::runNurbs()
{
   vector<real_t> N(20);
   clear(N);
   vector<int> x(30);
   clear(x);
   knot(_n1,_c1,x);
   size_t icount = 0;

// Calculate the points on the rational B-spline curve
   real_t t = 0;
   real_t step = real_t(x[_n1+_c1-1])/real_t(_np1-1);
   for (size_t i=0; i<_np1; ++i) {
      if (real_t(x[_n1+_c1-1])-t < FUNC_APPROX_THRESHOLD)
         t = real_t(x[_n1+_c1-1]);
      RationalBasis(_c1,t,_n1,x,N);
      for (size_t j=0; j<3; ++j) {
         size_t jcount = j;
         (*_p)[icount+j] = 0.;
         for (size_t k=0; k<_n1; ++k) {
            (*_p)[icount+j] += N[k]*(*_q)[jcount];
            jcount += 3;
         }
      }
      icount += 3;
      t += step;
   }
   return 0;
}


int FuncApprox::runNurbsSurface()
{
   vector<real_t> bold(4), niku(_n1*_np1), mjlw(_n2*_np2);
   vector<real_t> rsumij(_np1*_np2), savrsumij(_np1*_np2);
   clear(niku);
   clear(mjlw);
   for (size_t i=0; i<_np1*_np2; ++i)
      rsumij[i] = savrsumij[i] = 1.;
   int ret = 0;
   size_t ibnum = 0;
   int pm = 1, change = 4;
   for (size_t i=0; i<30000; ++i) {
      ret = rbsurf(ibnum,bold,niku,mjlw,rsumij,savrsumij);
      ibnum = 20;
      if (i==0) {
         bold[0] = _qq[ibnum  ];
         bold[1] = _qq[ibnum+1];
         bold[2] = _qq[ibnum+2];
         bold[3] = _qq[ibnum+3];
      }
      _qq[ibnum+3] += pm*change;
      pm = -pm;
   }
   return ret;
}


int FuncApprox::rbsurf(size_t          ibnum,
                       vector<real_t>& bold,
                       vector<real_t>& ni,
                       vector<real_t>& mj,
                       vector<real_t>& rsumij,
                       vector<real_t>& savrsumij)
{
   static int itest=0;
 
// Allows for 20 data points with basis function of order 5
   vector<int> x(30), y(30);
   vector<real_t> N(30), M(30);
   if (itest != int(_n1+_c1+_n2+_c2+_np1+_np2)) {
      clear(x);
      clear(y);
      clear(N);
      clear(M);
 
//    Generate the open uniform knot vectors
      knot(_n1,_c1,x);
      knot(_n2,_c2,y);

//    Calculate and store the basis functions
      real_t stepu = real_t(x[_n1+_c1-1])/real_t(_np1-1);
      real_t stepw = real_t(y[_n2+_c2-1])/real_t(_np2-1);
 
//    Calculate the Ni's at each u parametric value
      real_t u = 0., w = 0.;
      size_t j = 0;
      for (size_t i1=0; i1<_np1; ++i1) {
         if (real_t(x[_n1+_c1-1])-u < FUNC_APPROX_THRESHOLD)
            u = x[_n1+_c1-1];
         SRationalBasis(_c1,u,_n1,x,N);
         for (size_t i=0; i<_n1; ++i)
            ni[j++] = N[i];
         u += stepu;
      }

//    Calculate the Mj's at each w parametric value
      j = 0;
      for (size_t i2=0; i2<_np2; ++i2) {
         if (real_t(y[_n2+_c2-1])-w < FUNC_APPROX_THRESHOLD)
            w = y[_n2+_c2-1];
         SRationalBasis(_c2,w,_n2,y,M);
         for (size_t i=0; i<_n2; ++i)
            mj[j++] = M[i];
         w += stepw;
      }
 
//    Calculate the sum function at each parametric value
      j = 0;
      for (size_t i1=0; i1<_np1; ++i1) {
         for (size_t i=0; i<_n1; ++i)
            N[i] = ni[i1*_n1+i];
         for (size_t i2=0; i2<_np2; ++i2) {
            for (size_t i=0; i<_np2; ++i)
               M[i] = mj[i1*_n2+i];
            real_t s = sumrbas(N,M);
            if (s==0.)
               throw OFELIException("FuncAPprox::rbsurf(...): Sum of the basis functions = 0 This is most likely caused\n" 
                                    "by all zero homogeneous weighting factors."); 
            rsumij[j++] = 1./s;
         }
      }
      itest = 0;
   }

// Generate the complete rational B-spline surface
   if (itest==0) {
      clear(*_p);
      size_t icount = 0;
      for (size_t i1=0; i1<_np1; ++i1) {
         for (size_t i2=0; i2<_np2; ++i2) {
            size_t scount = i1*_np2 + i2;
            for (size_t i=0; i<_n1; ++i) {
               size_t jbas = 4*_n2*i, n = i1*_n1 + i;
               if (ni[n] != 0.) {
                  for (size_t j=0; j<_n2; ++j) {
                     size_t k = jbas + 4*j, m = i2*_n2 + j;
                     if (mj[m] !=0.) {
                        real_t pbasis = _qq[k+3]*ni[n]*mj[m]*rsumij[scount];
                        (*_p)[icount  ] += _qq[k  ]*pbasis;
                        (*_p)[icount+1] += _qq[k+1]*pbasis;
                        (*_p)[icount+2] += _qq[k+2]*pbasis;
                     }
                  }
               }
            }
            icount += 3;
         }
      }
      itest = _n1 + _c1 + _n2 + _c2 + _np1 + _np2;
      return 0;
   }
 
// Calculate the incremental change to the surface
   real_t bx = _qq[ibnum  ] - bold[0];
   real_t by = _qq[ibnum+1] - bold[1];
   real_t bz = _qq[ibnum+2] - bold[2];
   real_t bh = _qq[ibnum+3] - bold[3];

   if (bx != 0. || by != 0. || bz != 0. || bh != 0.) {
 
//    if true surface unchanged -- no calculation is needed
//    if not true surface has changed -- do incremental calculation 
 
//    Calculate the i,j index for Bij from ibnum, where ibnum is
//    assumed to be the lineal number of the x-component of Bij
      size_t iindex = (ibnum/4/_n2);     // depends on integer arithmetic to work
      size_t jindex = (ibnum-4*_n2*iindex)/4;

//    For the special case of the homogeneous weighting factor changing
//    the sum function must be recalculated for each value of u,w
      if (bh != 0.) {
         savrsumij = rsumij;
         size_t tcount = 0;
         for (size_t i1=0; i1<_np1; ++i1) {
            size_t n = i1*_n1 + iindex;
            for (size_t i2=0; i2<_np2; ++i2) {
               rsumij[tcount] = 1./(1./rsumij[tcount] + ni[n]*mj[i2*_n2+jindex]*(_qq[ibnum+3]-bold[3]));
               tcount++;
            }
         }
      }

//    Calculate the change in the surface for each u,w
      size_t icount = 0;
      for (size_t i1=0; i1<_np1; ++i1) {
         size_t n = i1*_n1 + iindex;
         if (ni[n] != 0.) {
            for (size_t i2=0; i2<_np2; ++i2) {
               size_t m = i2*_n2 + jindex;
               if (mj[m] != 0.) {
                  size_t scount = i1*_np2 + i2;
                  if (bh==0.) {
                     real_t pbasis = _qq[ibnum+3]*ni[n]*mj[m]*rsumij[scount];
                     (*_p)[icount  ] += bx*pbasis;
                     (*_p)[icount+1] += by*pbasis;
                     (*_p)[icount+2] += bz*pbasis;
                  }
                  else {
                     real_t pbasis = ni[n]*mj[m]*rsumij[scount];
                     real_t s = rsumij[scount]/savrsumij[scount];
                     (*_p)[icount  ] = (*_p)[icount  ]*s + bh*_qq[ibnum  ]*pbasis;
                     (*_p)[icount+1] = (*_p)[icount+1]*s + bh*_qq[ibnum+1]*pbasis;
                     (*_p)[icount+2] = (*_p)[icount+2]*s + bh*_qq[ibnum+2]*pbasis;
                  }
               }
               icount += 3;
            }
         }
         else
            icount += 3*_np2;
      }
   }
   bold[0] = _qq[ibnum  ];
   bold[1] = _qq[ibnum+1];
   bold[2] = _qq[ibnum+2];
   bold[3] = _qq[ibnum+3];
   return 0;
}


int FuncApprox::runLagrange()
{
   string p="", v=_ffct->getVar(1);
   for (int i=0; i<=_degree; ++i) {
      string q = "";
      real_t d = 1.;
      for (int j=0; j<=_degree; ++j) {
         if (i!=j)
            d *= (*_x)[i] - (*_x)[j];
      }
      d = (*_y)[i]/d;
      if (i>0) {
         string s = " + ";
         if (d<0)
            s = " - ", d = fabs(d);
         p = p + s;
      }
      for (int j=0; j<=_degree; ++j) {
         if (i!=j) {
            if ((*_x)[j]==0.0)
               q = q + v;
            else
               q = q + "(" + v + "-" + toString((*_x)[j]) + ")";
            if (j<_degree && (i<_degree || j<_degree-1))
               q = q + "*";
         }
      }
      p = p + toString(d) + "*" + q;
   }
   return _ffct->set(p);
}


int FuncApprox::runLeastSquare()
{
   real_t s=0., s1=0., s2=0.;
   if (_ls_opt<4) {
      for (size_t i=1; i<=_n2; ++i) {
         s = 0.;
         for (size_t k=1; k<=_np1; ++k)
            s += (*_B)(k,i)*(*_y)(i);
         _b(i) = s;
         for (size_t j=1; j<=_n2; ++j) {
            s = 0.;
            for (size_t k=1; k<=_n1; ++k)
               s += (*_B)(k,i)*(*_B)(k,j);
            _A(i,j) = s;
         }
      }
      _A.solveQR(_b,*_p);
   }
   else {
         s1 = s2 = 0.;
         for (size_t i=0; i<_n1; ++i) {
            s1 += (*_x)[i];
            s2 += (*_x)[i]*(*_x)[i];
         }
         s = 1./(_n1*s2 - s1*s1);
         *_a0 = s*(s2*(*_y)[0] - s1*(*_y)[1]);
         *_a1 = s*(_n1*(*_y)[1] - s1*(*_y)[0]);
   }
   return 0;
}


int FuncApprox::run()
{
   int ret = 0;
   switch (_ap) {

      case LAGRANGE:
         ret = runLagrange();
         break;

      case LEAST_SQUARE:
         ret = runLeastSquare();
         break;

      case BSPLINE:
         ret = runBSpline();
         break;

      case BSPLINE_SURFACE:
         ret = runBSplineSurface();
         break;

      case BEZIER:
         ret = runBezier();
         break;

      case BEZIER_SURFACE:
         ret = runBezierSurface();
         break;

      case NURBS:
         ret = runNurbs();
         break;

      case NURBS_SURFACE:
         ret = runNurbsSurface();
         break;

      default:
         ret = 1;
         break;
   }
   return ret;
}

} /* namespace OFELI */