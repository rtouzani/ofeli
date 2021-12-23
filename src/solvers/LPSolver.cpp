/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2022 Rachid Touzani

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

                       Implementation of class 'LPSolver'

  ==============================================================================*/

#include "solvers/LPSolver.h"
#include "solvers/simplex.h"
#include "linear_algebra/Vect_impl.h"
#include "OFELIException.h"

using std::to_string;

namespace OFELI {

LPSolver::LPSolver()
         : _ret(0)
{
}

void LPSolver::setSize(int nv,
                       int nb_le,
                       int nb_ge,
                       int nb_eq)
{
   _nv = nv;
   _nb_le = nb_le;
   _nb_ge = nb_ge;
   _nb_eq = nb_eq;
   _nc = _nb_le + _nb_ge + _nb_eq;
   _i_le = _i_ge = _i_eq = 0;
   _A.setSize(_nc+1,_nv+1);
}
		       

LPSolver::~LPSolver()
{
}


void LPSolver::set(Vect<real_t>& x)
{
   _x = &x;
}


void LPSolver::set(Setting             opt,
                   const Vect<real_t>& a,
                   real_t              b)
{
   if (opt==OBJECTIVE) {
      _A(1,1) = b;
      for (int i=2; i<=_nv+1; ++i)
         _A(1,i) = a(i-1);
   }
   else if (opt==LE_CONSTRAINT) {
      if (++_i_le>_nb_le)
         throw OFELIException("In LPSolver::set(...): Too many (<=) inequality constraints.");
      _A(_i_le+1,1) = b;
      for (int i=2; i<=_nv+1; ++i)
         _A(_i_le+1,i) = a(i-1);
   }
   else if (opt==GE_CONSTRAINT) {
      if (++_i_ge>_nb_ge)
         throw OFELIException("In LPSolver::set(...): Too many (>=) inequality constraints.");
      _A(_i_ge+_nb_le+1,1) = b;
      for (int i=2; i<=_nv+1; ++i)
         _A(_i_ge+_nb_le+1,i) = a(i-1);
   }
   else if (opt==EQ_CONSTRAINT) {
      if (++_i_eq>_nb_eq)
         throw OFELIException("In LPSolver::set(...): Too many equality constraints.");
      _A(_i_eq+_nb_le+_nb_ge+1,1) = b;
      for (int i=2; i<=_nv+1; ++i)
         _A(_i_eq+_nb_le+_nb_ge+1,i) = a(i-1);
   }
}


int LPSolver::run()
{
   _theSimplex.set(_A,_nv,_nb_le,_nb_ge,_nb_eq,*_x);
   _ret = _theSimplex.run();
   return _ret;
}


ostream& operator<<(ostream&        s,
                    const LPSolver& os)
{
   s << "\n\nSUMMARY OF LINEAR PROGRAMMING SOLVER:" << endl;
   s << "Number of optimization variables:" << setw(6) << os._nv << endl;
   s << "Number of (<=) constraints:\t\t" << setw(6) << os._nb_le << endl;
   s << "Number of (>=) constraints:\t\t" << setw(6) << os._nb_ge << endl;
   s << "Number of equality constraints:\t\t" << setw(6) << os._nb_eq << endl;
   s << "Objective:\t\t\t\t" << os.getObjective() << endl;
   return s;
}

} /* namespace OFELI */
