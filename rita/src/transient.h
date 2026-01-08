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

                        Definition of class 'transient'

  ==============================================================================*/

#pragma once

#include "mesh/Mesh.h"
#include "solvers/LinearSolver.h"
#include "solvers/TimeStepping.h"
#include "OFELI.h"
#include "rita.h"
#include "solve.h"

namespace RITA {

class ls;
class ae;
class ode;
class pde;

class transient
{

 public:

    transient(rita *r, vector<string> &pb);
    ~transient() { }
    int run(const vector<string>& pb);
    int ret;

 private:

    rita *_rita;
    data *_data;
    const string _pr = ">transient>";
    int _nb_ls, _nb_ae, _nb_ode, _nb_pde;
    ae *_ae;
    ode *_ode;
    pde *_pde;
    vector<string> _pb;
    vector<OFELI::ODESolver> _ode_eq;
    vector<OFELI::NLASSolver> _nlas_eq;
    vector<OFELI::TimeStepping> _ts_eq;
    void getPDEError(double &e2, double &eI);
};

} /* namespace RITA */
