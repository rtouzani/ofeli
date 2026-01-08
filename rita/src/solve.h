/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2026 Rachid Touzani

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

                            Definition of class 'solve'

  ==============================================================================*/

#pragma once

#include <string>
#include <iostream>
using std::string;

#include "rita.h"
#include "linear_algebra/Matrix.h"
#include "solvers/LinearSolver.h"
#include "solvers/OptSolver.h"
#include "solvers/EigenProblemSolver.h"
#include "io/Fct.h"


namespace RITA {

class configure;
class data;
class optim;
class eigen;

class solve
{

 public:

    solve();
    solve(rita *r, cmd *command, configure *config);
    ~solve() { }
    void setVerbose(int verb) { _verb = verb; }
    int getVerbose() const { return _verb; }
    void set(cmd* com) { _cmd = com; }
    int run();

 private:

    rita *_rita;
    data *_data;
    configure *_configure;
    cmd *_cmd;
    const string _pr = ">solve>";
    int _verb, _nb_pb, _nb_ls, _nb_ae, _nb_pde, _nb_opt, _nb_eigen;
    vector<string> _pb;
    analysis_type _analysis;    
    OFELI::Fct _theFct;
    optim *_optim;
    eigen *_eigen;
    ae *_ae;
    pde *_pde;
    ls *_ls;

    int run_pb();
    int run_steady();
    int run_optim();
    int run_eigen();
    DataType checkPb(string p);
    void getPDEError(double &e2, double &eI);
};

} /* namespace RITA */
