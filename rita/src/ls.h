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

                          Definition of class 'ls'

  ==============================================================================*/

#pragma once

#include "rita.h"
#include "io/Fct.h"
#include <map>

namespace RITA {
 
class ls
{
 public:

    ls(rita *r, cmd* command, configure* config);
    ~ls();
    int set();
    int run();
    void print(ostream& s) const;

    string solver_name, prec_name, name;
    OFELI::Iteration lls;
    OFELI::Preconditioner pprec;
    bool log;
    int iMat, iRHS, iSol, iInit, solved;

 private:
    
    rita *_rita;
    data *_data;
    int _ret, _verb;
    string _fn, _alg;
    const string _pr = "ls>";
    OFELI::Matrix<double> *_Mat;
    configure *_configure;
    cmd *_cmd;
};

ostream& operator<<(ostream& s, const ls& e);

} /* namespace RITA */
