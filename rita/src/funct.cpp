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

funct::funct()
      : _iv(0)
{
}


funct::funct(const string& s)
      : _iv(0)
{
   name = s;
   _theFct = new OFELI::Fct;
   _theFct->setName(s);
   nb_var = nb_par = 0;
}


void funct::set(OFELI::Fct& f)
{
   delete _theFct;
   _theFct = &f;
   setFromFct();
}


void funct::setFromFct()
{
   name = _theFct->getName();
   nb_var = _theFct->getNbVar();
   for (size_t i=1; i<=nb_var; ++i)
      var.push_back(_theFct->getVar(i));
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
   _theFct->setVar(_Var);
   for (const auto& p: par)
      _theFct->setVar(p);
   int ret = _theFct->set(s);
   if (ret)
      cout << "Error in function definition" << endl;
   return ret;
}


void funct::setVar(const vector<string>& v)
{
   for (size_t i=0; i<v.size(); ++i)
      var.push_back(v[i]);
}


void funct::setVar(const vector<string>& v, const vector<size_t>& n)
{
   for (size_t i=0; i<v.size(); ++i) {
      var.push_back(v[i]);
      setVar(n[i]);
   }
}


void funct::setVar(const string& v, int n)
{
   var.push_back(v);
   setVar(n);
}


void funct::setPar(const string& name)
{
   par.push_back(name);
   nb_par = par.size();
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
      _Var.push_back(var[_iv++]);
      nb_var++;
   }
   else {
      for (size_t j=1; j<=n; ++j)
         _Var.push_back(var[_iv]+to_string(j));
      nb_var += n;
      _iv++;
   }
}


double funct::operator()(double x)
{
   _arg.clear();
   _arg.push_back(x);
   return getPars();
}


double funct::operator()(double x, double y)
{
   _arg.clear();
   _arg.push_back(x);
   _arg.push_back(y);
   return getPars();
}


double funct::operator()(double x, double y, double z)
{
   _arg.clear();
   _arg.push_back(x);
   _arg.push_back(y);
   _arg.push_back(z);
   return getPars();
}


double funct::operator()(double x, double y, double z, double t)
{
   _arg.clear();
   _arg.push_back(x);
   _arg.push_back(y);
   _arg.push_back(z);
   _arg.push_back(t);
   return getPars();
}


double funct::operator()(const vector<double>& x)
{
   _arg.clear();
   for (size_t i=0; i<x.size(); ++i)
      _arg.push_back(x[i]);
   return getPars();
}


double funct::D(real_t x)
{
   return _theFct->D(x);
}


double funct::D(const vector<double>& x, size_t i)
{
   return _theFct->D(x,i);
}


double funct::getPars()
{
   for (size_t i=0; i<nb_par; ++i)
      _arg.push_back(par_value[i]);
   return (*_theFct)(_arg);
}


ostream& operator<<(ostream& s, const funct& f)
{
   s << "Function Name: " << f.name << endl;
   s << "Expression: " << f.getExpression() << endl;
   if (f.nb_var) {
      s << "Variable(s): " << f.var[0];
      for (size_t i=1; i<f.nb_var; ++i)
         s << "," << f.var[i];
   }
   if (f.nb_par) {
      s << ", Parameter(s): " << f.par[0];
      for (size_t i=1; i<f.nb_par; ++i)
         s << "," << f.par[i];
      s << endl;
   }
   return s;
}

} /* namespace RITA */