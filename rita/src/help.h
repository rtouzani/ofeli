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

                          Definition of class 'help'

  ==============================================================================*/

#pragma once

#include "rita.h"
#include "cmd.h"

#include <string>
#include <iostream>
#include <vector>

namespace RITA {

class help
{
 public:

const string H0 = "\n"
                  "DATA COMMANDS: vector, matrix, tabulation, mesh, grid, function, sample, plot, history, data, list,\n"
                  "               remove, rename, print, =\n"
                  "SOLVER COMMANDS: approximation, integration, ls, algebraic, ode, pde, stationary, transient, \n"
                  "                 optim, eigen, solve\n"
                  "GENERAL PURPOSE COMMANDS: help or ?, license, set, load, !, #, end or <, exit or quit\n";
const string H1 = "rita is an open-source numerical computational software for numerical analysis and simulation.\n"
                  "rita enables solving most numerical analysis problems using scientific computing tools.\n"
                  "More specifically, it handles interpolation, approximation, algebraic equations, ordinary \n"
                  "differential and partial differential equations.\n"
                  "The problems to solve in rita are defined by data and numerical tools to treat them.\n"
                  "SETTING COMMANDS: help, license, set\n"
                  "DATA COMMANDS: vector, matrix, tabulation, mesh, grid, function, sample, plot, history, data,\n"
                  "               remove, rename, list, print, =\n"
                  "SOLVER COMMANDS: approximation, integration, ls, algebraic, ode, pde, optim, eigen, solve\n"
                  "GENERAL PURPOSE COMMANDS: help or ?, license, set, load, !, #, end or <, exit or quit\n";
const string H2 = "rita is an open-source numerical computational software for numerical analysis and simulation.\n"
                  "rita enables solving most numerical analysis problems using scientific computing tools.\n"
                  "More specifically, it handles interpolation, approximation, algebraic equations, ordinary \n"
                  "differential and partial differential equations.\n"
                  "The problems to solve in rita are defined by data and numerical tools to treat them.\n\n"
                  "DATA COMMANDS:\n" 
                  "vector:     Define or read a vector\n"
                  "matrix:     Define or read a matrix\n"
                  "tabulation: Define a tabulation: (xi,yi) array\n"
                  "mesh:       Define a finite element mesh\n"
                  "grid:       Define a grid\n"
                  "function:   Define analytically a function\n"
                  "sample:     Sample a function in a given vector\n"
                  "plot:       Plot some types of data\n"
                  "history:    Define a history vector to store a time dependent solution of an ode or a pde\n"
                  "data:       List all data\n"
                  "list:       List a specific data type\n"
                  "remove:     Remove (delete) specified data or entity\n"
                  "rename:     Rename specified data or entity\n"
                  "print:      Display specific data\n"
                  "=:          Display specific data\n\n"
                  "SOLVER COMMANDS:\n"
                  "approximation: Define and solve a data approximation problem\n"
                  "integration:   Define and solve a numerical integration problem\n"
                  "ls:            Define and solve a linear system of equations\n"
                  "algebraic:     Define an algebraic equation (or set of equations)\n"
                  "ode:           Define an ordinary differential equation (or set of equations)\n"
                  "pde:           Define a partial differential equation (or set of equations)\n"
                  "optim:         Define an optimization problem\n"
                  "eigen:         Define an eigenvalue problem\n"
                  "solve:         Solve defined problem\n\n"
                  "SETTING COMMANDS: \n"
                  "help or ?:     Display this text or an abriged version\n"
                  "license:       Display rita license information\n"
                  "set:           Set configuration parameters\n"
                  "load:          Read commands in a given script file\n"
                  "!:             The string following this sign is executed as a shell command\n"
                  "#:             A line starting with # is considered as a comment line\n"
                  "end or <:      Go back to command level, if available\n"
                  "exit or quit:  End execution\n";

    help(rita* r, cmd* command);
    ~help() { }
    int run(int opt);

 private:
    rita *_rita;
    cmd *_cmd;
    size_t _nb_args;
    const string _pr = "help>";
    int Calc();
    int Vector();
    int Matrix();
    int mesh();
    int ls();
    int AE();
    int ODE();
    int PDE();
    int Opt();
    int Eig();
    int solve();
    int plot();
    int set1D();
    int setRectangle();
    int setCube();
    int setPoint();
    int setCurve();
    int setSurface();
    int setVolume();
    int setContour();
    int setCode();
    int MeshRead();
};

} /* namespace RITA */
