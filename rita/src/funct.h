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

    bool _alloc;
    size_t nb_var;
    int opt;
    string name;
    vector<string> var;
    funct(rita *r, const string& s);
    ~funct();
    void set(OFELI::Fct& f);
    int set(const string& v, const string& s);
    int setExpr(const string& s);
    void setVar(const string& v, int n=1);
    string getExpression() const { return _theFct->getExpression(); }
    double operator()(double x) { return (*_theFct)(x); }
    double operator()(double x, double y) { return (*_theFct)(x,y); }
    double operator()(double x, double y, double z) { return (*_theFct)(x,y,z); }
    double operator()(double x, double y, double z, double t) { return (*_theFct)(x,y,z,t); }
    double operator()(const OFELI::Point<double>& x) { return (*_theFct)(x); }
    double operator()(const OFELI::Point<double>& x, double t) { return (*_theFct)(x,t); }
    double operator()(const vector<double>& x) { return (*_theFct)(x); }
    double D(real_t x) { return _theFct->D(x); }
    double D(const vector<double>& x, size_t i) { return _theFct->D(x,i); }
    string getErrorMessage() const { return _theFct->getErrorMessage(); }
    OFELI::Fct *getFct() const { return _theFct; }
    OFELI::Fct& fct() const { return *_theFct; }

 private:

    rita *_rita;
    vector<string> _exvar;
    vector<size_t> _N;
    OFELI::Fct *_theFct;
    void setVar();
};

ostream& operator<<(ostream& s, const funct& f);

} /* namespace RITA */