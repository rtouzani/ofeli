/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2024 Rachid Touzani

    This file is part of rita.

    rita is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rita is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================

                        Implementation of class 'funct'

  ==============================================================================*/

#include "funct.h"
#include <iostream>
#include "rita.h"
#include "defs.h"

namespace RITA {

funct::funct(rita *r, const string& s)
{
   _rita = r;
   name = s;
   _theFct = new OFELI::Fct;
   _theFct->setName(s);
   nb_var = 0;
   opt = 0;
   _alloc = true;
}


funct::~funct()
{
//   if (_alloc)
//      delete _theFct;
}


void funct::set(OFELI::Fct& f)
{
    delete _theFct;
   _alloc = false;
   _theFct = &f;
   name = _theFct->getName();
   nb_var = _theFct->getNbVar();
   for (size_t i=1; i<=nb_var; ++i)
      var.push_back(_theFct->getVar(i));
   setExpr(_theFct->getExpression());
}


int funct::setExpr(const string& s)
{
   if (nb_var==0)
      return 1;
   setVar();
   int ret = _theFct->set(name,s,_exvar,1);
   if (ret)
      _rita->msg("","Error in function definition");
   return ret;
}


void funct::setVar(const string& v, int n)
{
   if (v=="t" && n==1)
      opt = 1;
   _N.push_back(n);
   var.push_back(v);
   nb_var++;
}


int funct::set(const string& v, const string& s)
{
   setVar(v);
   _theFct->set(s,v);
   return 0;
}


void funct::setVar()
{
   size_t k=0;
   for (size_t i=0; i<_N.size(); ++i) {
      if (_N[i]==1)
         _exvar.push_back(var[k++]);
      else {
         for (size_t j=0; j<_N[i]; ++j)
            _exvar.push_back(var[k]+to_string(j+1));
         k++;
      }
   }
}


ostream& operator<<(ostream& s, const funct& f)
{
   s << f.name << "(" << f.var[0];
   for (size_t i=1; i<f.nb_var; ++i)
      s << "," << f.var[i];
   s << ") = " << f.getExpression() << endl;
   return s;
}

} /* namespace RITA */