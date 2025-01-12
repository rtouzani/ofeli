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

                            Implementation of class 'Fct'

  ==============================================================================*/

#include "io/Fct.h"
#include <ostream>
#include <iostream>
#include "io/exprtk.hpp"

namespace OFELI {

Fct::Fct()
    : _p(nullptr), _st(nullptr), _ex(nullptr), _name("f"), _nb_var(0), _nb_par(0),
      _exp_ok(false), _var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
}


Fct::Fct(const string& exp)
    : _p(nullptr), _st(nullptr), _ex(nullptr), _name("f"), _nb_var(0), _nb_par(0),
      _exp_ok(false), _var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
   set(exp);
}


Fct::Fct(const string&         exp,
         const vector<string>& v)
    : _p(nullptr), _st(nullptr), _ex(nullptr), _name("f"),  _nb_var(0), _nb_par(0),
      _exp_ok(false), _var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
   set(exp,v);
}


Fct::Fct(const string& exp,
         const string& v)
    : _p(nullptr),  _st(nullptr), _ex(nullptr), _name("f"), _nb_var(0), _nb_par(0),
      _exp_ok(false), _var_ok(false), err(1)
{
   _name = "f";
   error_message = "No error in function evaluation.";
   set(exp,v);
}


Fct::Fct(const string&         n,
         const string&         exp,
         const vector<string>& v)
    : _p(nullptr),  _st(nullptr), _ex(nullptr), _name(n), _nb_var(0), _nb_par(0),
      _exp_ok(false), _var_ok(false), err(1)
{
   error_message = "No error in function evaluation.";
   set(exp,v);
}


Fct::~Fct()
{
   if (_p!=nullptr)
      delete _p;
   if (_ex!=nullptr)
      delete _ex;
   if (_st!=nullptr)
      delete _st;
}


void Fct::add_constants()
{
   _st->add_constant("pi",OFELI_PI);
   _st->add_constant("e",OFELI_E);
   _st->add_constants();
}


string Fct::getErrorMessage()
{
   return error_message;
}


int Fct::set(const string&         n,
             const string&         exp,
             const vector<string>& v,
             int                   opt)
{
   _name = n;
   int err = set(exp,v,opt);
   _exp_ok = 1 - err;
   return 1-err;

}


int Fct::getErrorCode() const
{
   return 1-err;
}


int Fct::set(const string& exp,
             const string& v,
             int           opt)
{
   if (_p!=nullptr)
      delete _p;
   _p = new exprtk::parser<real_t>;
   if (_st!=nullptr)
      delete _st;
   _st = new exprtk::symbol_table<real_t>;
   if (_ex!=nullptr)
      delete _ex;
   _ex = new exprtk::expression<real_t>;
   _exp_ok = _var_ok = true;
   _nb_var = 1;
   _expr = exp;
   _var.push_back(v);
   _xvar.resize(1);
   add_constants();
   _st->add_variable(_var[0],_xvar[0]);
   _ex->register_symbol_table(*_st);
   err = _p->compile(exp,*_ex);
   if (err==0) {
      error_message = _p->error();
      if (opt==0)
         std::cout << "Error: " << error_message << ", Expression: " << _expr << std::endl;
   }
   checkParam();
   _exp_ok = err;
   return 1-err;
}


int Fct::set(const string& exp,
             int           opt)
{
   if (_p!=nullptr)
      delete _p;
   _p = new exprtk::parser<real_t>;
   if (_st!=nullptr)
      delete _st;
   _st = new exprtk::symbol_table<real_t>;
   if (_ex!=nullptr)
      delete _ex;
   _ex = new exprtk::expression<real_t>;
   _exp_ok = _var_ok = true;
   _expr = exp;
   if (_nb_var==0) {
      _nb_var = 4;
      _var.push_back("x");
      _var.push_back("y");
      _var.push_back("z");
      _var.push_back("t");
   }
   _xvar.resize(_nb_var);
   for (size_t i=0; i<_nb_var; ++i)
      _st->add_variable(_var[i],_xvar[i]);
   add_constants();
   _ex->register_symbol_table(*_st);
   err = _p->compile(exp,*_ex);
   if (err==0) {
      error_message = _p->error();
      if (opt==0)
         std::cout << "Error: " << error_message << ", Expression: " << _expr << std::endl;
   }
   checkParam();
   _exp_ok = err;
   return 1-err;
}


void Fct::setVar(const string& v)
{
   _var.push_back(v);
   _nb_var++;
   _var_ok = true;
}


void Fct::setPar(const string& p,
                 real_t        x)
{
   for (size_t i=0; i<_nb_par; ++i) {
      if (_par[i]==p) {
         _xvar[i] = x;
         return;
      }
   }
   _nb_par++, _nb_var++;
   _par.push_back(p);
   _var.push_back(p);
   _xpar.push_back(x);
   _xvar.push_back(x);
}


void Fct::setVar(const vector<string>& v)
{
   for (size_t i=0; i<v.size(); ++i) {
      _var.push_back(v[i]);
      _nb_var++;
   }
   _var_ok = true;
}


int Fct::set(const string&         exp,
             const vector<string>& v,
             int                   opt)
{
   if (_p!=nullptr)
      delete _p;
   _p = new exprtk::parser<real_t>;
   if (_st!=nullptr)
      delete _st;
   _st = new exprtk::symbol_table<real_t>;
   if (_ex!=nullptr)
      delete _ex;
   _ex = new exprtk::expression<real_t>;
   _exp_ok = _var_ok = true;
   _nb_var = v.size();
   _expr = exp;
   for (auto it=std::begin(v); it!=std::end(v); ++it)
      _var.push_back(*it);
   _xvar.resize(_nb_var);
   add_constants();
   for (size_t i=0; i<_nb_var; ++i)
      _st->add_variable(_var[i],_xvar[i]);
   _ex->register_symbol_table(*_st);
   err = _p->compile(exp,*_ex);
   _exp_ok = err;
   if (err==0 && opt==0)
      std::cout << "Error: " << _p->error() << ", Expression: " << _expr << std::endl;
   checkParam();
   return 1-err;
}


real_t Fct::D(real_t x)
{
   _xvar[0] = x;
   return exprtk::derivative(*_ex,_xvar[0]);
}


real_t Fct::D(const vector<real_t>& x,
              size_t                i)
{
   _xvar = x;
   if (i<=_nb_var)
      return exprtk::derivative(*_ex,_xvar[i-1]);
   else
      return 0.;
}


int Fct::check()
{
   vector<string> variable_list;
   if (exprtk::collect_variables(_expr,variable_list)) {
      for (const auto& var : variable_list)
         std::cout << var << std::endl;
   }
   vector<string> v;
   if (exprtk::collect_variables(_expr,v))
      return 0;
   return 1;
}


void Fct::checkParam()
{
   bool ok=true;
   vector<string> par;
   exprtk::collect_variables(_expr,par);
   for (auto it=std::begin(_par); it!=std::end(_par); ++it) {
      ok = false;
      for (const auto& v : par) {
         if (v==*it)
            ok = true;
      }
      if (!ok)
         _par.erase(it,it+1), it--;
   }
   _nb_par = _par.size();
   for (auto it=std::begin(_var); it!=std::end(_var); ++it) {
      ok = false;
      for (const auto& v : par) {
         if (v==*it) {
            ok = true;
            continue;
         }
      }
      if (!ok)
         _var.erase(it,it+1), it--;
   }
   _rv = _var;
   for (auto it=std::begin(_rv); it!=std::end(_rv); ++it) {
      ok = false;
      for (const auto& v : _par) {
         if (v==*it) {
            ok = true;
            continue;
         }
      }
      if (ok)
         _rv.erase(it,it+1), it--;
   }
}


void Fct::getParam(vector<string>& p)
{
   p = _par;
}


real_t Fct::operator()(real_t x)
{
   _xvar[0] = x;
   return _ex->value();
}


real_t Fct::operator()(real_t x,
                       real_t y)
{
   _xvar[0] = x, _xvar[1] = y;
   return _ex->value();
}


real_t Fct::operator()(real_t x,
                       real_t y,
                       real_t z)
{
   _xvar[0] = x, _xvar[1] = y, _xvar[2] = z;
   return _ex->value();
}


real_t Fct::operator()(real_t x,
                       real_t y,
                       real_t z,
                       real_t t)
{
   _xvar[0] = x, _xvar[1] = y, _xvar[2] = z, _xvar[3] = t;
   return _ex->value();
}


real_t Fct::operator()(const SpaceTime& p)
{
   _xvar[0] = p.x, _xvar[1] = p.y, _xvar[2] = p.z, _xvar[3] = p.t;
   return _ex->value();
}


real_t Fct::operator()(const Point<real_t>& x)
{
   _xvar[0] = x.x, _xvar[1] = x.y, _xvar[2] = x.z;
   return _ex->value();
}


real_t Fct::operator()(const Point<real_t>& x,
                       real_t               t)
{
   _xvar[0] = x.x, _xvar[1] = x.y, _xvar[2] = x.z, _xvar[3] = t;
   return _ex->value();
}


real_t Fct::operator()(const vector<real_t>& x)
{
   _xvar = x;
   return _ex->value();
}


string Fct::getExpression() const
{
   return _expr;
}


///////////////////////////////////////////////////////////////////////////////
//                           ASSOCIATED  FUNCTIONS                           //
///////////////////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& s,
                         const Fct&    f)
{
   if (!f._var_ok) {
      s << f._name << ": Undefined variable(s)" << std::endl;
      return s;
   }
   if (!f._exp_ok) {
      s << f._name << ": Function undefined." << std::endl;
      return s;
   }
   s << f._name << "(" << f._rv[0];
   for (size_t i=1; i<f._rv.size(); ++i)
      s << "," << f._rv[i];
   s << ") = " << f._expr << std::endl;
   if (f._nb_par) {
      s << "Parameter(s): " << f._par[0];
      for (size_t i=1; i<f._nb_par; ++i)
         s << ", " << f._par[i];
      s << std::endl;
   }
   return s;
}


Fct operator+(const Fct& f,
              const Fct& g)
{
   Fct s;
   return s;
}

} /* namespace OFELI */