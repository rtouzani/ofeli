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

                            Definition of class 'rita'

  ==============================================================================*/

#pragma once

#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>

#include <iostream>
using std::ostream;
using std::cin;
using std::cout;
using std::endl;
using std::cerr;
using std::string;

#include "../muparserx/mpParser.h"
#include "../muparserx/mpDefines.h"
#include "../muparserx/mpTest.h"

#include "config.h"
#include "help.h"
#include "ritaException.h"
#include "data.h"
#include "OFELI.h"
#include "defs.h"

const string sPrompt = "rita>";

namespace RITA {
/*!
 *  \addtogroup RITA
 *  @{
 */

class cmd;
class calc;
class configure;
class optim;
class integration;
class approximation;
class ls;
class ae;
class ode;
class pde;
class eigen;
class help;
class solve;
class mesh;

enum class analysis_type {
   NONE          = 0,
   APPROXIMATION = 1,
   INTEGRATION   = 2,
   STEADY_STATE  = 3,
   TRANSIENT     = 4,
   EIGEN         = 5,
   OPTIMIZATION  = 6
};

enum objective_type { ANALYTIC_FUNCTION, PDE_BASED };


class rita
{

 public:

    rita();
    ~rita();
    int run();
    void setRelease(string r, string yr);
    void setVerbose(int verb) { _verb = verb; }
    void setInput(string file, int opt=1);
    void initConfig();
    bool getEcho() const { return _echo; }
    void setEcho(size_t nb_args);
    bool meshOK, dataOK;
    help *hh;
    string release, year, p2s;
    ofstream *ofh, *ofl, ocf;
    void finish();

    friend class configure;
    friend class help;
    friend class mesh;
    friend class solve;
    friend class transient;
    friend class optim;
    friend class eigen;
    friend class integration;
    friend class approximation;
    friend class data;
    friend class calc;
    friend class ls;
    friend class ae;
    friend class ode;
    friend class pde;
    friend class funct;

 private:

   rita *_rita;
   bool _echo, _load, _obj_analytic;
   string _script_file, _scheme, _sLine;
   ifstream _icf, *_in;
   cmd *_cmd;
   int _verb, _opt;
   double _init_time, _time_step, _final_time;
   int _adapted_time_step;
   size_t _nb_args;
   bool _analysis_ok;
   int _dim;
   analysis_type _analysis_type;
   configure *_configure;
   data *_data;
   calc *_calc;
   mesh *_mesh;
   ae *_ae;
   ode *_ode;
   ls *_ls;
   pde *_pde;
   solve *_solve;
   optim *_optim;
   eigen *_eigen;
   approximation *_approx;
   integration *_integration;

   void set(cmd* com) { _cmd = com; }
   void setDim(int dim) { _dim = dim; }
   void set(data *d) { _data = d; }
   int runLS();
   int runPDE();

   void getLicense();
   int Load();
   void setStationary();
   void setTransient();
   int findVector(const string& s);
   int CheckKeywords();
   void ListConst();
   void ListVar();
   void ListExprVar();
   void msg(const string& loc, const string& m1, const string& m2="", int c=0);

   const vector<string> _rita_kw {"load","eigen","optim","approx$imation","integ$ration","alg$ebraic","ae",
                                  "ode","pde","ls","solve"};
   const vector<string> _gkw {"?","help","lic$ense","set","end","<","echo"};
   const vector<string> _data_kw {"grid","mesh","vect$or","tab$ulation","func$tion","matr$ix","sample",
                                  "save","del$ete","rem$ove","desc$ription","hist$ory","data","list",
                                  "print","=","ren$ame","plot","%"};

   map<string,analysis_type> _ant = {{"",analysis_type::NONE},
                                     {"approximation",analysis_type::APPROXIMATION},
                                     {"integration",analysis_type::INTEGRATION},
                                     {"stationary",analysis_type::STEADY_STATE},
                                     {"transient",analysis_type::TRANSIENT},
                                     {"eigen",analysis_type::EIGEN},
                                     {"optimization",analysis_type::OPTIMIZATION}};

};

} /* namespace RITA */
