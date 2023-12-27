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

                      Implementation of class 'LeastSquare'

  ==============================================================================*/

#include "solvers/LeastSquare.h"
#include "linear_algebra/Vect_impl.h"

namespace OFELI {

LeastSquare::LeastSquare()
         : _n(0), _dim(1)
{
}


LeastSquare::LeastSquare(const vector<Fct &>& f, const Vect<real_t>& x, const Vect<real_t>& y, Vect<real_t>& a)
{
   set(f,x,y,a);
}


LeastSquare::LeastSquare(const DMatrix<real_t>& B, const Vect<real_t>& y, Vect<real_t>& a)
{
   set(B,y,a);
}


LeastSquare::LeastSquare(const Vect<real_t>& x, const Vect<real_t>& y, real_t& a0, real_t& a1)
{
   set(x,y,a0,a1);
}


void LeastSquare::set(const vector<Fct &>& f, const Vect<real_t>& x, const Vect<real_t>& y, Vect<real_t>& a)
{
   _f = &f;
   _x = &x;
   _y = &y;
   _a = &a;
   _
}


void LeastSquare::set(const DMatrix<real_t>& B, const Vect<real_t>& y, Vect<real_t>& a)
{

}


int LeastSquare::run()
{

}
