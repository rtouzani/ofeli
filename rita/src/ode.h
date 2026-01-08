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

                          Definition of class 'ode'

  ==============================================================================*/

#pragma once

#include "rita.h"
#include "funct.h"
#include <map>

namespace RITA {

#define NO_ODE  cout << "No ordinary differential equation defined." << endl;

class ode
{
 public:

    ode(rita *r, cmd* command, configure* config);
    ~ode() { }
    int set();
    int run();
    void print(ostream& s) const;

    bool isSet, log, isFct;
    DataType type;
    vector<int> iFct;
    vector<string> vars, analytic, analytic_f;
    vector<funct *> SolFct;
    OFELI::Vect<string> J;
    int size, solved, ind_fct, every, nb_var, get_err;
    size_t nb_vars;
    string nls, name, scheme, err;
    vector<string> var_name, phase;
    vector<int> ivect, iphase;
    double init_time, time_step, final_time, adapted_time_step;
    OFELI::TimeScheme Scheme;


 private:

    rita *_rita;
    data *_data;
    int _ret, _verb;
    string _fn, _alg;
    configure *_configure;
    cmd *_cmd;
    const string _pr = ">ode>";

    map<string,OFELI::TimeScheme> _sch = {{"forward-euler",OFELI::FORWARD_EULER},
                                          {"backward-euler",OFELI::BACKWARD_EULER},
                                          {"crank-nicolson",OFELI::CRANK_NICOLSON},
                                          {"heun",OFELI::HEUN},
                                          {"newmark",OFELI::NEWMARK},
                                          {"leap-frog",OFELI::LEAP_FROG},
                                          {"AB2",OFELI::ADAMS_BASHFORTH},
                                          {"RK4",OFELI::RK4},
                                          {"RK3-TVD",OFELI::RK3_TVD},
                                          {"BDF2",OFELI::BDF2},
                                          {"builtin",OFELI::BUILTIN}};

};

ostream& operator<<(ostream& s, const ode& e);

} /* namespace RITA */
