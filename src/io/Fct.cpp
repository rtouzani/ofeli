/*==============================================================================

                                    O  F  E  L  I

                           Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

                            Implementation of class 'Fct'

  ==============================================================================*/

#include "io/Fct.h"
#include <ostream>
#include <iostream>

extern exprtk::parser<real_t> theParser;

namespace OFELI {

Fct::Fct()
    : name("f"), nb_var(0), exp_ok(false), var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
}


Fct::Fct(string exp)
   : name("f"), nb_var(0), exp_ok(false), var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
   set(exp);
}


Fct::Fct(string          exp,
         vector<string>& v)
    : name("f"), nb_var(0), exp_ok(false), var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
   set(exp,v);
}


Fct::Fct(string exp,
         string v)
    : name("f"), nb_var(0), exp_ok(false), var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
   set(exp,v);
}


Fct::Fct(string          n,
         string          exp,
         vector<string>& v)
    : name("f"), nb_var(0), exp_ok(false), var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
   set(n,exp,v);
}


Fct::~Fct()
{
}


string Fct::getErrorMessage()
{
   return error_message;
}


  int Fct::set(string exp,
	       int    opt)
{
   exp_ok = true;
   expr = exp;
   err = 1;
   if (var_ok) {
      add_constants(symbol_table);
      expression.register_symbol_table(symbol_table);
      err = theParser.compile(expr,expression);
      if (err==0) {
         error_message = theParser.error();
         if (opt==0)
            std::cout << "Error: " << error_message << ", Expression: " << expr << std::endl;
      }
   }
   return 1-err;
}


int Fct::set(string          n,
             string          exp,
             vector<string>& v,
             int             opt)
{
   name = n;
   return set(exp,v,opt);
}


int Fct::set(string exp,
             string v,
             int    opt)
{
   exp_ok = var_ok = true;
   nb_var = 1;
   expr = exp;
   var.push_back(v);
   xvar.push_back(0.);
   add_constants(symbol_table);
   symbol_table.add_variable(var[0],xvar[0]);
   expression.register_symbol_table(symbol_table);
   err = theParser.compile(exp,expression);
   if (err==0) {
      error_message = theParser.error();
      if (opt==0)
         std::cout << "Error: " << error_message << ", Expression: " << expr << std::endl;
   }
   return 1-err;
}


int Fct::set(string          exp,
             vector<string>& v,
             int             opt)
{
   exp_ok = var_ok = true;
   nb_var = v.size();
   expr = exp;
   for (auto it=std::begin(v); it!=std::end(v); ++it) {
      var.push_back(*it);
      xvar.push_back(0.);
   }
   add_constants(symbol_table);
   for (size_t i=0; i<nb_var; ++i)
      symbol_table.add_variable(var[i],xvar[i]);
   expression.register_symbol_table(symbol_table);
   err = theParser.compile(exp,expression);
   if (err==0) {
      error_message = theParser.error();
      if (opt==0)
         std::cout << "Error: " << error_message << ", Expression: " << expr << std::endl;
   }
   return 1-err;
}


int Fct::set(vector<string>& v,
             int             opt)
{
   var_ok = true;
   nb_var = var.size();
   for (auto it=std::begin(v); it!=std::end(v); ++it) {
      var.push_back(*it);
      xvar.push_back(0.);
   }
   add_constants(symbol_table);
   for (size_t i=0; i<nb_var; ++i)
      symbol_table.add_variable(var[i],xvar[i]);
   expression.register_symbol_table(symbol_table);
   if (exp_ok) {
      err = theParser.compile(expr,expression);
      if (err==0) {
         error_message = theParser.error();
         if (opt==0)
            std::cout << "Error: " << error_message << ", Expression: " << expr << std::endl;
      }
   }
   return 1-err;
}


real_t Fct::D(real_t x)
{
   xvar[0] = x;
   return exprtk::derivative(expression,xvar[0]);
}


real_t Fct::D(const vector<real_t>& x,
              size_t                i)
{
   xvar = x;
   if (i<=nb_var)
      return exprtk::derivative(expression,xvar[i-1]);
   else
      return 0.;
}


int Fct::check()
{
   vector<string> v;
   if (exprtk::collect_variables(expr,v))
      return 0;
   return 1;
}


real_t Fct::operator()(real_t x)
{
   xvar[0] = x;
   return expression.value();
}


real_t Fct::operator()(real_t x,
		       real_t y)
{
   xvar[0] = x, xvar[1] = y;
   return expression.value();
}


real_t Fct::operator()(real_t x,
		       real_t y,
		       real_t z)
{
   xvar[0] = x, xvar[1] = y, xvar[2] = z;
   return expression.value();
}


real_t Fct::operator()(real_t x,
		       real_t y,
		       real_t z,
		       real_t t)
{
   xvar[0] = x, xvar[1] = y, xvar[2] = z, xvar[3] = t;
   return expression.value();
}


real_t Fct::operator()(const vector<real_t>& x)
{
   xvar = x;
   return expression.value();
}

  
std::ostream& operator<<(std::ostream& s,
                         const Fct&    f)
{
   s << "Name of function: " << f.name << std::endl;
   s << "Definition: " << f.expr << std::endl;
   s << "List of variables: ";
   for (size_t i=0; i<f.nb_var-1; ++i)
      s << f.var[i] << ", ";
   s << f.var[f.nb_var-1] << std::endl;
   return s;
}

} /* namespace OFELI */
