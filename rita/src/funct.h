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

                        Definition of class 'funct'

  ==============================================================================*/

#pragma once

#include <string>
using std::string;

#include "io/Fct.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Point.h"
#include "OFELI.h"

namespace RITA {

class rita;

class funct
{

 public:

    size_t nb_var, nb_par;
    string name;
    vector<string> var, par;
    vector<double> par_value;
    funct();
    funct(const string& s);
    ~funct() { }
    void set(OFELI::Fct& f);
    void setNoPar();
    int set(const string& v, const string& s);
    int setExpr(const string& s);
    void setVar(const string& v, int n=1);
    void setVar(const vector<string>& v);
    void setVar(const vector<string>& v, const vector<size_t>& n);
    void setPar(const string& name);
    void setParValue(const double& v);
    void setFromFct();
    string getExpression() const { return _theFct->getExpression(); }
    double operator()(double x);
    double operator()(double x, double y);
    double operator()(double x, double y, double z);
    double operator()(double x, double y, double z, double t);
    double operator()(const vector<double>& x);
    double D(real_t x);
    double D(const vector<double>& x, size_t i);
    string getErrorMessage() const { return _theFct->getErrorMessage(); }
    OFELI::Fct *getFct() const { return _theFct; }

 private:

    size_t _iv;
    vector<string> _Var;
    vector<double> _arg;
    OFELI::Fct *_theFct;
    void setVar(size_t n);
    double getPars();
};

ostream& operator<<(ostream& s, const funct& f);

} /* namespace RITA */