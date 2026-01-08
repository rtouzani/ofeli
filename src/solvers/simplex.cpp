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

                      Implementation of class 'simplex'

  ==============================================================================*/

#include <iostream>
using std::ostream;
using std::endl;
using std::cerr;

#include "solvers/simplex.h"
#include "linear_algebra/Vect_impl.h"

#include <iostream>
using std::ostream;
using std::endl;

namespace OFELI {

simplex::simplex()
        :  _ret(0), _eps(DBL_EPSILON)
{
}


simplex::simplex(Vect<real_t>& A,
                 int           nv,
                 int           nb_le,
                 int           nb_ge,
                 int           nb_eq,
                 Vect<real_t>& x)
        :  _ret(0), _eps(DBL_EPSILON)
{
   set(A,nv,nb_le,nb_ge,nb_eq,x);
}


simplex::~simplex()
{
   for (int i=0; i<_nb+2; ++i)
      delete _A[i];
   delete [] _A;
}


void simplex::set(Vect<real_t>& A,
                  int           nv,
                  int           nb_le,
                  int           nb_ge,
                  int           nb_eq,
                  Vect<real_t>& x)
{
   _nv = nv;
   _nb_le = nb_le;
   _nb_ge = nb_ge;
   _nb_eq = nb_eq;
   _nl2 = _nb = nb_le + nb_ge + nb_eq;
   _A = new real_t* [_nb+3];
   for (int i=0; i<_nb+3; ++i)
      _A[i] = new real_t [_nv+2];
   for (int j=1; j<_nv+2; ++j)
      _A[1][j] = -A(1,j);
   for (int i=2; i<_nb+2; ++i) {
      _A[i][1] = A(i,1);
      for (int j=2; j<_nv+2; ++j)
         _A[i][j] = -A(i,j);
   }
   _l1.resize(_nv);
   _l2.resize(_nb);
   _l3.resize(_nb);
   _zerov.resize(_nv);
   _posv.resize(_nb);

   for (int i=0; i<_nv; ++i)
      _l1[i] = _zerov[i] = i + 1;
   for (int i=1; i<=_nb; i++) { 
      if (_A[i+1][1]<0.0) {
         throw OFELIException("In simplex::simplex(...): Constants in constraints must be nonnegative");
         _ret = 1;
      }
      _l2[i-1] = i;
      _posv[i-1] = _nv + i; 
   }
   for (int i=0; i<_nb_ge+1; ++i)
      _l3[i] = 1;
   _x = &x;
}

  
void simplex::setSolution()
{
   if (_ret)
      return;
   for (int i=1; i<=_nv; i++) {
      (*_x)(i) = 0.;
      for (int j=1; j<=_nb; j++) {
         if (_posv[j-1]==i) {
            (*_x)(i) = _A[j+1][1];
            break;
         }
      }
   }
}


void simplex::simp1(int     mm,
                    int     iabf,
                    int&    kp,
                    real_t& bmax)
{
   real_t test=0.0;
   kp = _l1[0];
   bmax = _A[mm+1][kp+1];
   if (_nv < 2)
      return;
   for (int i=1; i<_nv; ++i) {
      if (iabf==0)
         test = _A[mm+1][_l1[i]+1] - bmax;
      else
         test = fabs(_A[mm+1][_l1[i]+1]) - fabs(bmax);
      if (test>0.0) {
         bmax = _A[mm+1][_l1[i]+1]; 
         kp = _l1[i];
      }
   }
}


int simplex::simp2(int     kp,
                   real_t& qq)
{
// Locate a pivot element, taking degeneracy into account. 
   int ip = 0;
   if (_nl2 < 1)
      return ip;
   int i = 0;
   for (i=1; i<=_nl2; i++) 
      if (_A[i+1][kp+1] < -_eps)
         goto e2; 
   return ip;  // No possible pivots. Return with message. 

e2:
   qq = -_A[_l2[i-1]+1][1]/_A[_l2[i-1]+1][kp+1]; 
   ip = _l2[i-1];
   if (i+1 > _nl2)
      return ip;
   int q0=0, qp=0;
   for (i=i+1; i<=_nl2; i++) { 
      int ii = _l2[i-1];
      if (_A[ii+1][kp+1] < -_eps) { 
         real_t q = -_A[ii+1][1]/_A[ii+1][kp+1]; 
         if (q<qq)
            ip = ii, qq = q;
         else if (q==qq) {  // We have a degeneracy.
            for (int k=1; k<=_nv; ++k) { 
               real_t qp = -_A[ip+1][k+1]/_A[ip+1][kp+1];
               q0 = -_A[ii+1][k+1]/_A[ii+1][kp+1];
               if (q0 != qp)
                  goto e6; 
            }
e6:         if (q0<qp)
               ip = ii;
         }
      }
   }
   return ip;
}


void simplex::simp3(int ii,
		    int ip,
		    int kp)
{ 
   real_t pivot = 1.0/_A[ip+1][kp+1];
   if (ii >= 0) {
      for (int i=1; i<=ii+1; ++i) {
         if (i-1 != ip) {
            _A[i][kp+1] *= pivot; 
            for (int j=1; j<=ii+1; ++j) 
               if (j-1 != kp)
                  _A[i][j] -= _A[ip+1][j]*_A[i][kp+1]; 
         }
      }
   }
   for (int i=1; i<=ii+1; ++i)
      if (i-1 != kp)
         _A[ip+1][i] = -_A[ip+1][i]*pivot; 
   _A[ip+1][kp+1] = pivot;
}


int simplex::run()
{
   int kh, kp, is, ip, ir=0;
   real_t bmax, q1;
   if (_nb_ge+_nb_eq==0)
      goto e30;
   ir = 1;

// Compute the auxiliary objective function.
   for (int k=1; k<=_nv+1; ++k) { 
      q1 = 0.0;
      for (int i=_nb_le+1; i<=_nb; ++i)
         q1 += _A[i+1][k]; 
      _A[_nb+2][k] = -q1; 
   }

e10:
   simp1(_nb+1,0,kp,bmax);  // Find max. coeff. of auxiliary objective fn 
   if (bmax<=_eps && _A[_nb+2][1]<-_eps) { 
      _ret = -1;            // Auxiliary objective function is still negative and can’t be improved, 
      return _ret;          // hence no feasible solution exists.
   }
   else if (bmax<=_eps && _A[_nb+2][1]<=_eps) {
      int m=_nb_le+_nb_ge+1;
      if (m<=_nb) {
         for (ip=m; ip<=_nb; ip++) {
            if (_posv[ip-1]==ip+_nv) {    // Found an artificial variable for an equalityconstraint. 
               simp1(ip,1,kp,bmax); 
               if (bmax>_eps)         // Exchange with column corresponding to maximum
	          goto e1;            // pivot element in row.
            }                  
         }
      }
      ir = 0;
      m--;
      if (_nb_le+1>m)
	 goto e30;

//    Change sign of row for any m2 constraints
      for (int i=_nb_le+1; i<=_nb_le+_nb_ge; ++i) {
         if (_l3[i-_nb_le-1]==1)             // still present from the initial basis. 
            for (int k=1; k<=_nv+1; ++k) 
               _A[i+1][k] *= -1.0;
      }
      goto e30;           // Go to phase two. 
   }

   ip = simp2(kp,q1);       // Locate a pivot element (phase one).
   if (ip == 0) {         // Maximum of auxiliary objective function is
      _ret = -1;          // unbounded, so no feasible solution exists.
      return _ret; 
   }

e1:
   simp3(_nb+1,ip,kp); 
// Exchange a left- and a right-hand variable (phase one), then update lists. 
   if (_posv[ip-1]>=_nv+_nb_le+_nb_ge+1) { // Exchanged out an artificial variable for an 
                                           // equality constraint. Make sure it stays 
                                           // out by removing it from the l1 list.
      int k = 0;
      for (k=1; k<=_nv; k++)
         if (_l1[k-1]==kp)
	    goto e2; 
e2:
      for (int i=k; i<_nv; ++i)
         _l1[i-1] = _l1[i]; 
   }	
   else {
      if (_posv[ip-1] < _nv+_nb_le+1)
         goto e20;
      kh = _posv[ip-1] - _nb_le - _nv;
      if (_l3[kh-1]==0)
         goto e20;     // Exchanged out an m2 type constraint. 
      _l3[kh-1] = 0;   // If it’s the first time, correct the pivot column  or the minus
                       // sign and the implicit artificial variable. 
   }
   _A[_nb+2][kp+1] += 1.0;
   for (int i=1; i<=_nb+2; i++)
      _A[i][kp+1] *= -1.0;

e20:
   is = _zerov[kp-1];             // Update lists of left- and right-hand variables. 
   _zerov[kp-1] = _posv[ip-1]; 
   _posv[ip-1] = is;
   if (ir != 0)
      goto e10;       // if still in phase one, go back to 10. 
// End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it. 

e30:
   simp1(0,0,kp,bmax);             // Test the z-row for doneness. 
   if (bmax <= _eps) {             // Done. Solution found. Return with the good news. 
      _ret = 0; 
      setSolution();
      return _ret;
   }
   ip = simp2(kp,q1);                // Locate a pivot element (phase two). 
   if (ip==0) {                    // Objective function is unbounded. Report and return. 
      _ret = 1;
      return _ret;
   }
   simp3(_nb,ip,kp);               // Exchange a left- and a right-hand variable (phase two), 
   goto e20;                       // Update lists of left- and right-hand variables and 
                                   // return for another iteration.
}


real_t simplex::getObjective() const
{
   if (_ret==0)
      return _A[1][1];
   else
      return 0;
}

} /* namespace OFELI */
