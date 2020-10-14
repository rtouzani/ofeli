/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                      Implementation of class 'Simplex'

  ==============================================================================*/

#include "solvers/Simplex.h"
#include "linear_algebra/DMatrix_impl.h"
#include "linear_algebra/Vect_impl.h"
#include <iostream>
#include <cmath>

namespace OFELI {

Simplex::Simplex(const DMatrix<real_t>& A,
                 Vect<real_t>&          b,
                 const Vect<real_t>&    c,
                 Vect<real_t>&          x)
        : _nb_it(0), _b(&b), _x(&x), _max(0.), _isUnbounded(false)
{
   size_t m = A.getNbRows(), n = A.getNbColumns();
   _nr = m, _nc = n + m;
   _A.setSize(_nr,_nc);
   _c.setSize(_nc);
   for (size_t i=1; i<=n; ++i)
      _c(i) = c(i);
   for (size_t i=1; i<=_nr; ++i) {
      for (size_t j=1; j<=n; ++j)
         _A(i,j) = A(i,j);
      for (size_t j=n+1; j<=_nc; ++j)
         _A(i,j) = 0.;
      _A(i,i+n) = 1.;
   }
}


bool Simplex::solve()
{
// Check whether the table is optimal,if optimal no need to process further
   if (checkOptimality())
      return true;

// Find the column which has the pivot.
// The least coefficient of the objective function.
   size_t pc = getPivotColumn();
   if (_isUnbounded) {
      throw OFELIException("In Simple::solve(): Unbounded");
      return true;
   }

// Find the row with the pivot value.
// The least value item's row in the array b
   size_t pr = getPivotRow(pc);

// Form the next table according to the pivot value
   doPivoting(pr,pc);
   return false;
}


bool Simplex::checkOptimality()
{
// if the table has further negative constraints,then it is not optimal
   bool isOptimal = false;
   size_t positiveValueCount = 0;

// check if the coefficients of the objective function are negative
   for (size_t i=0; i<_nc; ++i) {
      real_t value = _c[i];
      if (value >= 0)
         positiveValueCount++;
   }
// if all the constraints are positive now,the table is optimal
   if (positiveValueCount == _nc)
      isOptimal = true;
   return isOptimal;
}


void Simplex::doPivoting(size_t pr,
			 size_t pc)
{
   real_t pivotValue = _A(pr+1,pc+1);
   Vect<real_t> pivotRowVals(_nc), pivotColVals(_nr), rowNew(_nc);
   _max -= (_c[pc]*((*_b)[pr]/pivotValue));
   for (size_t i=0; i<_nc; ++i)
      pivotRowVals[i] = _A(pr+1,i+1);

// get the column that has the pivot value
   for (size_t j=0; j<_nr; ++j)
     pivotColVals[j] = _A(j+1,pc+1);

//set the row values that has the pivot value divided by the pivot value and put into new row
   for (size_t k=0; k<_nc; ++k)
      rowNew[k] = pivotRowVals[k]/pivotValue;
   (*_b)[pr] /= pivotValue;

// process the other coefficients in the A array by subtracting
   for (size_t m=0; m<_nr; ++m) {
//    ignore the pivot row as we already calculated that
      if (m != pr) {
         for (size_t p=0; p<_nc; ++p) {
            real_t multiplyValue = pivotColVals[m];
            _A(m+1,p+1) -= multiplyValue*rowNew[p];
//C[p] = C[p] - (multiplyValue*C[pivotRow]);
//B[i] = B[i] - (multiplyValue*B[pivotRow]);
         }
      }
   }

// Process the values of the B array
   for (size_t i=0; i<_nr; ++i) {
      if (i != pr) {
         real_t multiplyValue = pivotColVals[i];
         (*_b)[i] -= multiplyValue*(*_b)[pr];
      }
   }

// The least coefficient of the constraints of the objective function
   real_t multiplyValue = _c[pc];
// process the C array
   for (size_t i=0; i<_nc; ++i)
      _c[i] -= multiplyValue*rowNew[i];

// Replace the pivot row in the new calculated A array
   for (size_t i=0; i<_nc; ++i)
      _A(pr+1,i+1) = rowNew[i];
}


// Find the least coefficients of constraints in the objective function's position
size_t Simplex::getPivotColumn()
{
   size_t loc = 0;
   real_t minm = _c[0];
   for (size_t i=1; i<_nc; ++i) {
      if (_c[i]<minm) {
         minm = _c[i];
         loc = i;
      }
   }
   return loc;
}


//find the row with the pivot value.The least value item's row in the B array
size_t Simplex::getPivotRow(size_t pivotColumn)
{
   Vect<real_t> positiveValues(_nr), result(_nr);
   int negativeValueCount = 0;
   for (size_t i=0; i<_nr; ++i) {
      if (_A(i+1,pivotColumn+1)>0)
         positiveValues[i] = _A(i+1,pivotColumn+1);
      else {
         positiveValues[i] = 0;
         negativeValueCount += 1;
      }
   }

// Checking the unbound condition if all the values are negative ones
   if (negativeValueCount==_nr)
      _isUnbounded = true;
   else {
      for (size_t i=0; i<_nr; ++i) {
         result[i] = 0;
         real_t value = positiveValues[i];
         if (value>0)
            result[i] = (*_b)[i]/value;
      }
   }
// Find the minimum's location of the smallest item of the array b
   real_t mm = 99999999;
   size_t loc = 0;
   for (size_t i=0; i<_nr; ++i) {
      if (result[i]>0) {
         if (result[i]<mm) {
            mm = result[i];
            loc = i;
         }
      }
   }
   return loc;
}


int Simplex::run()
{
   bool end=false;
   while (!end) {
      bool result = solve();
      _nb_it++;
      if (result)
         end = true;
   }
   for (size_t i=0; i<_nc; ++i) {
      size_t count0=0, index=0;
      for (size_t j=0; j<_nr; ++j) {
         if (_A(j+1,i+1)==0.0)
            count0 += 1;
         else if (_A(j+1,i+1)==1)
            index = j;
      }
      (*_x)[i] = 0.;
      if (count0==_nr-1)
         (*_x)[i] = (*_b)[index];
   }
   return _nb_it;
}

} /* namespace OFELI */
