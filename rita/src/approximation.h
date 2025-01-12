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

                       Definition of class 'approximation'

  ==============================================================================*/

#pragma once

#include "OFELI.h"
#include "io/Tabulation.h"
#include "rita.h"
#include "cmd.h"
#include "configure.h"
#include "funct.h"
#include <map>

namespace RITA {
 
class approximation
{

 public:

    enum ApproxType {
       LAGRANGE,
       FITTING,
       BSPLINE,
       BEZIER,
       NURBS,
       NONE
    };

    enum FitType {
       POLYNOMIAL,
       EXPONENTIAL,
       DEFINED
    };

    approximation(rita* r, cmd* command, configure* config);
    ~approximation();
    int set();
    int run();

 private:

    rita *_rita;
    configure *_configure;
    cmd *_cmd;
    data *_data;
    string _PR = "approximation>", _pr;
    OFELI::Tabulation *_theTab;
    ApproxType _method;
    FitType ft;
    Vect<double> _x, *_y;
    int _verb, _degree, _nb_fit, _file_count, _tab_count, _basis_count, _iTab;
    vector<string> _fct_names;
    vector<Fct *> _fct;
    funct *_ffct;
    OFELI::Fct *_ff;
    string _fft, _tab, _file;
    bool _tab_alloc;
    int setBSpline();
    int setBezier();
    int setNurbs();
    int setLagrange();
    int setFitting();
    int lagrange_eval(double x);
    void setXY();
    map<ApproxType,string> rApp = {{LAGRANGE,"lagrange"},
                                   {FITTING,"fitting"},
                                   {BSPLINE,"bspline"},
                                   {BEZIER,"bezier"},
                                   {NURBS,"nurbs"}};
};

} /* namespace RITA */
