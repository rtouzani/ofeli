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

                          Definition of class 'optim'

  ==============================================================================*/

#pragma once

#include "solvers/OptSolver.h"
#include "rita.h"
#include "solve.h"
#include "funct.h"
#include <map>

namespace RITA {

class optim
{
 public:

    optim(rita *r, cmd* command, configure* config);
    ~optim();
    int set();
    int run();
    OFELI::OptSolver::OptMethod Alg;
    string name, opt_name;
    vector<string> var_name;
    int solved, size, nb_eqc, nb_lec, nb_gec, igrad, ihess, iincons, ieqcons, verbose;
    funct *J_Fct;
    bool G_ok, H_ok, log, lp;
    double penal, b, obj;
    vector<funct *> G_Fct, H_Fct, inC_Fct, eqC_Fct;
    vector<double> init;
    vector<OFELI::Vect<double> *> a_le, a_ge, a_eq;
    OFELI::Vect<double> *opt_var, lb, ub, a, b_eq, b_ge, b_le;
    void print(ostream& s) const;

 private:
    
    rita *_rita;
    int _verb;
    string _fn, _alg;
    const string _pr = ">optimization>";
    configure *_configure;
    cmd *_cmd;
    data *_data;

    map<string,OFELI::OptSolver::OptMethod> Nopt = {{"gradient",OFELI::OptSolver::GRADIENT},
                                                    {"truncated-newton",OFELI::OptSolver::TRUNCATED_NEWTON},
                                                    {"simulated-annealing",OFELI::OptSolver::SIMULATED_ANNEALING},
                                                    {"nelder-mead",OFELI::OptSolver::NELDER_MEAD},
                                                    {"newton",OFELI::OptSolver::NEWTON}};
};

ostream& operator<<(ostream& s, const optim& o);

} /* namespace RITA */