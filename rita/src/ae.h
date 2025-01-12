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

                          Definition of class 'ae'

  ==============================================================================*/

#pragma once

#include "rita.h"
#include "funct.h"
#include <map>

namespace RITA {

class ae
{
 public:

    ae(rita *r, cmd* command, configure* config);
    ~ae() { }
    int set();
    int run();
    void print(ostream& s) const;

    bool isSet, log;
    vector<string> analytic, vars;
    vector<int> iFct;
    size_t nb_vars;
    OFELI::Vect<string> J;
    int solved, size, ind_fct, every, isFct;
    vector<int> ivect;
    string nls, fn, name;
    vector<string> var_name;
    NonLinearIter nnls;

 private:
    
    rita *_rita;
    data *_data;
    int _verb;
    configure *_configure;
    cmd *_cmd;
    string _PR = "", _pr;
};

ostream& operator<<(ostream& s, const ae& e);

} /* namespace RITA */
