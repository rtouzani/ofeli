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

                           Definition of class 'pde'

  ==============================================================================*/

#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <map>
using std::map;

#include "rita.h"
#include "data.h"

#include "LinearPDE.h"
#include "Laplace.h"
#include "Therm.h"
#include "Solid.h"
#include "Fluid.h"
#include "Electromagnetics.h"
#include "funct.h"

using namespace OFELI;

namespace RITA {

class cmd;
class rita;

#define MAX_NB_VECTORS 6
#define NO_PDE  cout << "No partial differential equation defined." << endl;

class pde
{

 public:

    struct Log {
      bool pde, vect, spd, ls, nl, mesh;
      Log() { pde = vect = spd = ls = nl = mesh = false; }
      bool fail() const { return (pde || vect || spd || ls || nl || mesh); }
    };

    struct PdeData {
       map<int,string> cexp;
       string exp, ft_name, in_file, out_file;
       int size;
       PdeData() { exp=""; in_file=""; out_file=""; size=0; };
    };

    struct VectorData {
       vector<Vect<double> *> theVect;
       int                    vect, nb_dof;
       string                 fn;
       data::DataSize         ds; 
    };

    pde(rita* r, cmd* command, configure* config);
    ~pde();
    int run();
    int set_nls(string nls);
    int set_ls(string ls, string prec);

    enum pde_eq {
       LINEAR_PDE,
       LAPLACE,
       HEAT,
       WAVE,
       TRANSPORT,
       LINEAR_ELASTICITY,
       TRUSS,
       BEAM,
       INCOMPRESSIBLE_NAVIER_STOKES,
       COMPRESSIBLE_EULER,
       COMPRESSIBLE_NAVIER_STOKES,
       INCOMPRESSIBLE_POROUS_1PHASE,
       EDDY_CURRENTS,
       MAXWELL,
       HELMHOLTZ,
       EIKONAL
    };

    enum sdm {
       FD,
       FE_P1,
       FE_P2,
       FE_Q1,
       FV,
       DG
    };

    int solved, nb_vectors, nb_dof, ieq, verbose, every, get_err;
    sdm Sdm;
    string name, eq, ls, nls, file, lsolv, lprec, prec, scheme, e2, eI;
    bool axi;
    PdeData in_data, bc_data, bf_data, sf_data;
    OFELI::Iteration lls;
    OFELI::Preconditioner pprec;
    double final_time, time_step;
    vector<string> analytic, analytic_f;
    vector<funct *> SolFct;
    Mesh *theMesh;
    Grid *theGrid;
    VectorData fd[MAX_NB_VECTORS];
    Equa *theEquation;
    analysis_type analysis;
    void setVectors();
    void set(string e);
    int setSpD(string spd);
    int set(Grid* gr);
    int setEq();
    void set();
    void set(data *d);
    int setCoef();
    int getIn();
    int getBC();
    int getBF();
    int getSF();
    int setIn();
    int setBC();
    int setBF();
    int setSF();
    void check();
    void set(cmd* cmd) { _cmd = cmd; }
    int setNodeBC(int code, string exp, double t, Vect<double>& v);
    int setSize(Vect<double>& v, data::DataSize s);
    Log log;
    bool set_u, set_bc, set_bf, set_sf, set_in, set_coef;
    Vect<double> u, b, bc, bf, sf;
    OFELI::TimeScheme Scheme;
    string regex_u;
    map<string,sdm> pde_sdm {{"fd",FD},
                             {"feP1",FE_P1},
                             {"feP2",FE_P2},
                             {"feQ1",FE_Q1},
                             {"fv",FV},
                             {"dg",DG}};
    map<string,OFELI::TimeScheme> _Sch = {{"forward-euler",OFELI::FORWARD_EULER},
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
    void print(ostream& s) const;

 private:

    rita *_rita;
    configure *_configure;
    cmd *_cmd;
    int _verb, _dim, _ret;
    bool _c00_set, _c10_set, _c01_set, _c20_set, _c02_set;
    bool _rho_set, _Cp_set, _kappa_set, _mu_set, _sigma_set, _Mu_set, _epsilon_set, _omega_set;
    bool _beta_set, _v_set, _young_set, _poisson_set;
    const vector<string> _var {"t","x","y","z"};
    map<string,int> pde_map = {{"linear-pde",LINEAR_PDE},
                               {"laplace",LAPLACE},
                               {"heat",HEAT},
                               {"wave",WAVE},
                               {"transport",TRANSPORT},
                               {"linear-elasticity",LINEAR_ELASTICITY},
                               {"truss",TRUSS},
                               {"beam",BEAM},
                               {"incompressible-navier-stokes",INCOMPRESSIBLE_NAVIER_STOKES}, 
                               {"compressible-euler",COMPRESSIBLE_EULER},
                               {"incompressible-porous-1phase",INCOMPRESSIBLE_POROUS_1PHASE}};

    const vector<string> _kw = {"expression","value","file","save"};
    const string _pr = "pde>";
    data *_data;
    string _c00, _c10, _c01, _c20, _c02;
    string _rho_exp, _Cp_exp, _kappa_exp, _mu_exp,_sigma_exp, _Mu_exp, _epsilon_exp, _omega_exp;
    string _beta_exp, _v_exp, _young_exp, _poisson_exp;
    vector<string> _vect;
    OFELI::Fct _theFct;
    void setPDECoef();
};

ostream& operator<<(ostream& s, const pde& e);

} /* namespace RITA */
