/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2025 Rachid Touzani

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
      : _rita(r), _iv(0)
{
   name = s;
   _theFct = new OFELI::Fct;
   _theFct->setName(s);
   nb_var = nb_par = 0;
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
}


void funct::set(const vector<string>& pn, const vector<double>& pv)
{
   for (size_t i=1; i<pv.size(); ++i)
      _theFct->setPar(pn[i],pv[i]);
   nb_par = _theFct->getNbPar();
   setExpr(_theFct->getExpression());
}


void funct::setNoPar()
{
   nb_par = 0;
   setExpr(_theFct->getExpression());
}


int funct::setExpr(const string& s)
{
   if (nb_var==0)
      return 1;
   _theFct->setVar(_exvar);
   int ret = _theFct->set(name,s,1);
   if (ret)
      _rita->msg("","Error in function definition");
   return ret;
}


int funct::setExpr(const string& s, const vector<string>& pn, const vector<double>& pv)
{
   if (nb_var==0)
      return 1;
   _theFct->setVar(_exvar);
   for (size_t i=1; i<pv.size(); ++i) {
      var.push_back(pn[i]);
      par.push_back(pn[i]);
      _theFct->setPar(pn[i],pv[i]);
      nb_par++;
   }
   int ret = _theFct->set(s);
   if (ret)
      _rita->msg("","Error in function definition");
   return ret;
}


void funct::setVar(const string& v, int n)
{
   var.push_back(v);
   setVar(n);
}


int funct::set(const string& v, const string& s)
{
   setVar(v);
   _theFct->set(s,v);
   return 0;
}


void funct::setVar(size_t n)
{
   if (n==1) {
      _exvar.push_back(var[_iv++]);
      nb_var++;
   }
   else {
      for (size_t j=1; j<=n; ++j)
         _exvar.push_back(var[_iv]+to_string(j));
      nb_var += n;
      _iv++;
   }
}


ostream& operator<<(ostream& s, const funct& f)
{
   s << f.fct();
   return s;
}

} /* namespace RITA */