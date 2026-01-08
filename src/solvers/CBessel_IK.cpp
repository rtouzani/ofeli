/*==============================================================================
 
                                   O  F  E  L  I
 
                            Object  Finite  Element  Library
 
  ==============================================================================

  Copyright (C) 1998 - 2026 Rachid Touzani
 
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
 
       Implementation of functions to compute modified Bessel functions
                        In, Kn and their derivatives

     Algorithms and coefficient values from "Computation of Special Functions", 
     Zhang and Jin, John Wiley and Sons, 1996.

    (C) 2003, C. Bond. All rights reserved.
 
 ==============================================================================*/

//  CBessel_IK.cpp -- complex modified Bessel functions.

#include "solvers/Bessel.h"

namespace OFELI {

int CBessel_IK_01(complex_t z, complex_t &ci0, complex_t &ci1, complex_t &ck0,
                  complex_t &ck1, complex_t &ci0p, complex_t &ci1p, complex_t &ck0p,
                  complex_t &ck1p)
{
   complex_t zr, zr2, cr, ca, cb, cs, ct, cw;
   real_t w0;
   int k, kz;
   static real_t a[] = {
                         0.125,
                         7.03125e-2,
                         7.32421875e-2,
                         1.1215209960938e-1,
                         2.2710800170898e-1,
                         5.7250142097473e-1,
                         1.7277275025845,
                         6.0740420012735,
                         2.4380529699556e1,
                         1.1001714026925e2,
                         5.5133589612202e2,
                         3.0380905109224e3
                        };
   static real_t b[] = {
                        -0.375,
                        -1.171875e-1,
                        -1.025390625e-1,
                        -1.4419555664063e-1,
                        -2.7757644653320e-1,
                        -6.7659258842468e-1,
                        -1.9935317337513,
                        -6.8839142681099,
                        -2.7248827311269e1,
                        -1.2159789187654e2,
                        -6.0384407670507e2,
                        -3.3022722944809e3
                        };
   static real_t a1[] = {
                         0.125,
                         0.2109375,
                         1.0986328125,
                         1.1775970458984e1,
                         2.1461706161499e2,
                         5.9511522710323e3,
                         2.3347645606175e5,
                         1.2312234987631e7,
                         8.401390346421e08,
                         7.2031420482627e10
                        };

   real_t a0 = abs(z);
   complex_t z2=z*z, z1=z;
   if (a0 == 0.0) {
      ci0 = complex_t(1,0);
      ci1 = complex_t(0,0);
      ck0 = ck1 = complex_t(OFELI_INFINITY,0);
      ci0p = complex_t(0,0);
      ci1p = complex_t(0.5,0.0);
      ck0p = ck1p = complex_t(-OFELI_INFINITY;,0);
      return 0;
   }
   if (z.real() < 0.0)
      z1 = -z;
   if (a0 <= 18.0) {
      ci0 = cr = complex_t(1,0);
      for (k=1; k<=50; k++) {
         cr *= 0.25*z2/real_t(k*k);
         ci0 += cr;
         if (abs(cr/ci0) < OFELI_EPSMCH)
            break;
      }
      ci1 = cr = complex_t(1,0);
      for (size_t k=1; k<=50; k++) {
         cr *= 0.25*z2/real_t(k*(k+1.0));
         ci1 += cr;
         if (abs(cr/ci1) < OFELI_EPSMCH)
            break;
      }
      ci1 *= 0.5*z1;
   }
   else {
      if (a0 >= 50.0)
         kz = 7;
      else if (a0 >= 35.0)
         kz = 9;
      else
         kz = 12;
      ca = exp(z1)/sqrt(2.0*M_PI*z1);
      ci0 = complex_t(1,0);
      zr = 1.0/z1;
      for (k=0; k<kz; k++)
         ci0 += a[k]*pow(zr,k+1.0);
      ci0 *= ca;
      ci1 = complex_t(1,0);
      for (k=0; k<kz; k++)
         ci1 += b[k]*pow(zr,k+1.0);
      ci1 *= ca;
   }
   if (a0 <= 9.0) {
      cs = complex_t(0,0);
      ct = -log(0.5*z1) - EULER_CONST;
      w0 = 0.0;
      cr = complex_t(1,0);
      for (size_t k=1; k<=50; k++) {
         w0 += 1.0/k;
         cr *= 0.25*z2/real_t(k*k);
         cs += cr*(w0+ct);
         if (abs((cs-cw)/cs) < OFELI_EPSMCH)
            break;
         cw = cs;
      }
      ck0 = ct + cs;
   }
   else {
      cb = 0.5/z1;
      zr2 = 1.0/z2;
      ck0 = complex_t(1,0);
      for (k=0; k<10; k++)
         ck0 += a1[k]*pow(zr2,k+1.0);
      ck0 *= cb/ci0;
   }
   ck1 = (1.0/z1 - ci1*ck0)/ci0;
   if (z.real() < 0.0) {
      if (z.imag() < 0.0) {
         ck0 += complex_t(0,1)*M_PI*ci0;
         ck1 = -ck1 + complex_t(0,1)*M_PI*ci1;
      }
      else if (z.imag() > 0.0) {
         ck0 -= complex_t(0,1)*M_PI*ci0;
         ck1 = -ck1 - complex_t(0,1)*M_PI*ci1;
      }
      ci1 = -ci1;
   }
   ci0p = ci1;
   ci1p = ci0 - 1.0*ci1/z;
   ck0p = -ck1;
   ck1p = -ck0 - 1.0*ck1/z;
   return 0;
}


int CBessel_IK_na(complex_t z, int &nm, CVector &ci, CVector &ck, CVector &cip, CVector &ckp)
{
   complex_t ci0, ci1, ck0, ck1, ckk, cf, cf1, cf2, cs;
   int ecode, n=ci.size();
   real_t a0 = abs(z);
   nm = n;
   if (a0 < 1.0e-100) {
      for (size_t k=0; k<=n; ++k) {
         ci[k]  = cip[k] = complex_t(0,0);
         ck[k]  = complex_t(-OFELI_INFINITY,0);
         ckp[k] = complex_t( OFELI_INFINITY,0);
      }
      ci[0]  = complex_t(1,0);
      cip[1] = complex_t(0.5,0.0);
      return 0;
   }
   ecode = CBessel_IK_01(z,ci[0],ci[1],ck[0],ck[1],cip[0],cip[1],ckp[0],ckp[1]);
   if (n < 2)
      return 0;
   ci0 = ci[0]; ci1 = ci[1];
   ck0 = ck[0]; ck1 = ck[1];
   int m = msta1(a0,200);
   if (m < n)
      nm = m;
   else
      m = msta2(a0,n,15);
   cf2 = complex_t(0,0);
   cf1 = complex_t(1.0e-100,0.0);
   for (int k=m; k>=0; k--) {
      cf = 2.0*(k+1.0)*cf1/z + cf2;
      if (k <= nm)
         ci[k] = cf;
      cf2 = cf1; cf1 = cf;
   }
   cs = ci0/cf;
   ci *= cs;
   for (size_t k=2; k<=nm; k++) {
      if (abs(ci[k-1]) > abs(ci[k-2]))
         ckk = (1.0/z-ci[k]*ck[k-1])/ci[k-1];
      else
         ckk = (ci[k]*ck[k-2]+2.0*(k-1.0)/(z*z))/ci[k-2];
      ck[k] = ckk;
   }
   for (k=2; k<=nm; k++) {
      cip[k] =  ci[k-1] - real_t(k)*ci[k]/z;
      ckp[k] = -ck[k-1] - real_t(k)*ck[k]/z;
   }
   return 0;
}


int CBessel_IK_nb(complex_t z, int &nm, CVector &ci, CVector &ck, CVector &cip, CVector &ckp)
{
   complex_t cbs, csk0, cf, cf0, cf1, ca0, cbkl, cg, cg0, cg1, cs0, cs, cr;
   real_t vt, fac, a0=abs(z);
   int k, kz, l, m, n=ci.size();
   nm = n;
   if (a0 < 1.0e-100) {
      ci = complex_t(0,0);
      ck = complex_t(std::numeric_limits<real_t>::infinity(),0);
      cip = complex_t(0,0);
      ckp = complex_t(-std::numeric_limits<real_t>::infinity(),0);
      ci[0] = complex_t(1.0,0.0);
      cip[1] = complex_t(0.5,0.0);
      return 0;
   }
   complex_t z1 = z;
   if (z.real() < 0.0)
      z1 = -z;
   if (n == 0)
      nm = 1;
   m = msta1(a0,200);
   if (m < nm)
      nm = m;
   else
      m = msta2(a0,nm,15);
   cbs = csk0 = cf0 = complex_t(0,0);
   cf1 = complex_t(1.0e-100,0.0);
   for (k=m; k>=0; k--) {
      cf = 2.0*(k+1.0)*cf1/z1 + cf0;
      if (k <=nm)
         ci[k] = cf;
      if ((k != 0) && (k == 2*(k>>1)))
         csk0 += 4.0*cf/real_t(k);
      cbs += 2.0*cf;
      cf0 = cf1;
      cf1 = cf;
   }
   cs0 = exp(z1)/(cbs-cf);
   for (k=0; k<=nm; k++)
       ci[k] *= cs0;
   if (a0 <= 9.0) {
      ck[0] = -(log(0.5*z1)+EULER_CONST)*ci[0]+cs0*csk0;
      ck[1] = (1.0/z1-ci[1]*ck[0])/ci[0];
   }
   else {
     ca0 = sqrt(M_PI_2/z1)*exp(-z1);
     if (a0 >= 200.0)
        kz = 6;
     else if (a0 >= 80.0)
        kz = 8;
     else if (a0 >= 25.0)
        kz = 10;
     else
        kz = 16;
     for (l=0; l<2; l++) {
         cbkl = complex_t(1,0);
         vt = 4.0*l;
         cr = complex_t(1,0);
         for (size_t k=1; k<=kz; k++) {
            cr *= 0.125*(vt-pow(2.0*k-1.0,2.0))/(real_t(k)*z1);
            cr *= 0.125*(vt-(2*k-1)*(2*k-1))/(real_t(k)*z1);
            cbkl += cr;
         }
         ck[l] = ca0*cbkl;
      }
   }
   cg0 = ck[0], cg1 = ck[1];
   for (k=2; k<=nm; k++) {
      cg = 2.0*(k-1.0)*cg1/z1 + cg0;
      ck[k] = cg;
      cg0 = cg1; cg1 = cg;
   }
   if (z.real() < 0.0) {
      fac = 1.0;
      for (k=0; k<=nm; k++) {
         if (z.imag() < 0.0)
            ck[k] = fac*ck[k] + complex_t(0,1)*M_PI*ci[k];
         else
            ck[k] = fac*ck[k] - complex_t(0,1)*M_PI*ci[k];
         ci[k] *= fac;
         fac = -fac;
      }
   }
   cip[0] = ci[1], ckp[0] = -ck[1];
   for (size_t k=1; k<=nm; k++) {
      cip[k] =  ci[k-1] - real_t(k)*ci[k]/z;
      ckp[k] = -ck[k-1] - real_t(k)*ck[k]/z;
   }
   return 0;
}


int CBessel_IK_v(real_t v, complex_t z, real_t &vm, CVector &civ, CVector &ckv, CVector &civp, CVector &ckvp)
{
   complex_t ca1, ca, cs, cr, ci0, cbi0, cf, cf1, cf2;
   complex_t ct, cp, cbk0, ca2, cr1, cr2, csu, cws, cb;
   complex_t cg0, cg1, cgk, cbk1, cvk, v0p, v0n, w0, gap, gan;
   int m, k, kz;
   complex_t z1=z, z2=z*z;
   int n = v;
   real_t v0=v-n, a0=abs(z);
   real_t piv = M_PI*v0, vt = 4.0*v0*v0;
   if (n == 0)
      n = 1;
   if (a0 < 1e-100) {
      for (size_t k=0; k<=n; k++) {
         civ[k] = civp[k] = complex_t(0,0);
         ckv[k] = complex_t(-std::numeric_limits<real_t>::infinity(),0);
         ckvp[k] = complex_t(std::numeric_limits<real_t>::infinity(),0);
      }
      if (v0 == 0.0) {
         civ[0] = complex_t(1,0);
         civp[1] = complex_t(0.5,0.0);
      }
      vm = v;
      return 0;
   }
   if (a0 >= 50.0)
      kz = 8;
   else if (a0 >= 35.0)
      kz = 10;
   else
      kz = 14;
   if (z.real() <= 0.0)
      z1 = -z;
   if (a0 < 18.0) {
      if (v0 == 0.0)
         ca1 = complex_t(1,0);
      else {
         gap = gamma(1.+v0);
         ca1 = pow(0.5*z1,v0)/gap;
      }
      ci0 = cr = complex_t(1,0);
      for (k=1; k<=50; k++) {
         cr *= 0.25*z2/(k*(k+v0));
         ci0 += cr;
         if (abs(cr/ci0) < OFELI_EPSMCH)
            break;
      }
      cbi0 = ci0*ca1;
   }
   else {
      ca = exp(z1)/sqrt(2.0*M_PI*z1);
      cs = cr = complex_t(1,0);
      for (k=1; k<=kz; k++) {
         cr *= -0.125*(vt-(2*k-1)*(2*k-1))/(real_t(k)*z1);
         cs += cr;
      }
      cbi0 = ca*cs;
   }
   m = msta1(a0,200);
   if (m < n)
      n = m;
   else
      m = msta2(a0,n,15);
   cf2 = complex_t(0,0), cf1 = complex_t(1.0e-100,0.0);
   for (int k=m; k>=0; k--) {
      cf = 2.0*(v0+k+1.0)*cf1/z1 + cf2;
      if (k <= n)
         civ[k] = cf;
      cf2 = cf1;
      cf1 = cf;
   }
   cs = cbi0/cf;
   for (k=0; k<=n; k++)
      civ[k] *= cs;
   if (a0 <= 9.0) {
      if (v0 == 0.0) {
         ct = -log(0.5*z1) - el;
         cs = complex_t(0,0);
         w0 = 0.0;
         cr = complex_t(1,0);
         for (size_t k=1; k<=50; k++) {
            w0 += 1.0/k;
            cr *= 0.25*z2/(real_t)(k*k);
            cp = cr*(w0+ct);
            cs += cp;
            if ((k >= 10) && (abs(cp/cs) < OFELI_EPSMCH))
               break;
         }
         cbk0 = ct + cs;
      }
      else {
         gan = gamma(1.-v0);
         ca2 = 1.0/(gan*pow(0.5*z1,v0));
         ca1 = pow(0.5*z1,v0)/gap;
         csu = ca2 - ca1;
         cr1 = cr2 = complex_t(1,0);
         cws = complex_t(0,0);
         for (size_t k=1; k<=50; k++) {
            cr1 *= 0.25*z2/(k*(k-v0));
            cr2 *= 0.25*z2/(k*(k+v0));
            csu += ca2*cr1-ca1*cr2;
            if ((k >= 10) && (abs((cws-csu)/csu) < OFELI_EPSMCH))
               break;
            cws = csu;
         }
         cbk0 = csu*M_PI_2/sin(piv);
      }
   }
   else {
      cb = exp(-z1)*sqrt(M_PI_2/z1);
      cs = cr = complex_t(1,0);
      for (size_t k=1; k<=kz; k++) {
         cr *= 0.125*(vt-(2*k-1)*(2*k-1))/(real_t(k)*z1);
         cs += cr;
      }
      cbk0 = cb*cs;
   }
   cbk1 = (1.0/z1-civ[1]*cbk0)/civ[0];
   ckv[0] = cbk0; ckv[1] = cbk1;
   cg0 = cbk0; cg1 = cbk1;
   for (size_t k=2; k<=n; k++) {
      cgk = 2.0*(v0+k-1.0)*cg1/z1 + cg0;
      ckv[k] = cgk;
      cg0 = cg1; cg1 = cgk;
   }
   if (z.real() < 0.0) {
      for (size_t k=0; k<=n; k++) {
         cvk = exp((k+v0)*M_PI*complex_t(0,1));
         if (z.imag() < 0.0) {
            ckv[k] = cvk*ckv[k] + M_PI*complex_t(0,1)*civ[k];
            civ[k] /= cvk;
         }
         else if (z.imag() > 0.0) {
            ckv[k] = ckv[k]/cvk - M_PI*complex_t(0,1)*civ[k];
            civ[k] *= cvk;
         }
      }
   }
   civp[0] = v0*civ[0]/z + civ[1];
   ckvp[0] = v0*ckv[0]/z - ckv[1];
   for (size_t k=1; k<=n; k++) {
      civp[k] = -(k+v0)*civ[k]/z + civ[k-1];
      ckvp[k] = -(k+v0)*ckv[k]/z - ckv[k-1];
   }
   vm = n + v0;
   return 0;
}

} /* namespace OFELI */
