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

                          Implementation of class 'rita'

  ==============================================================================*/

#include <fstream>
#include "rita.h"
#include "data.h"
#include "calc.h"
#include "cmd.h"
#include "mesh.h"
#include "ls.h"
#include "ae.h"
#include "ode.h"
#include "pde.h"
#include "solve.h"
#include "optim.h"
#include "eigen.h"
#include "integration.h"
#include "approximation.h"
#include "configure.h"

using std::cout;
using std::exception;

namespace RITA {

rita::rita()
     : meshOK(false), dataOK(false), _echo(false), _load(false), _script_file(""),
       _in(nullptr), _verb(1), _analysis_type(analysis_type::NONE), _mesh(nullptr),
       _ae(nullptr), _ode(nullptr), _ls(nullptr), _pde(nullptr), _solve(nullptr),
       _optim(nullptr), _eigen(nullptr), _approx(nullptr), _integration(nullptr)
{
   _rita = this;
   _cmd = new cmd(this);
   _configure = new configure(this,_cmd);
   _data = new data(this,_cmd,_configure);
   _calc = new calc(this,_cmd);
   hh = new help(this,_cmd);
   ofl = _configure->getOStreamLog();
   ofh = _configure->getOStreamHistory();
}


rita::~rita()
{
   _data->remove_temp();
   if (_in!=nullptr)
      delete _in;
   delete _cmd;
   delete _data;
   delete _calc;
   delete hh;
   delete _configure;
   if (_solve!=nullptr)
      delete _solve;
   if (_mesh!=nullptr)
      delete _mesh;
   if (_ls!=nullptr)
      delete _ls;
   if (_ae!=nullptr)
      delete _ae;
   if (_ode!=nullptr)
      delete _ode;
   if (_pde!=nullptr)
      delete _pde;
   if (_optim!=nullptr)
      delete _optim;
   if (_eigen!=nullptr)
      delete _eigen;
   if (_integration!=nullptr)
      delete _integration;
   if (_approx!=nullptr)
      delete _approx;
}


void rita::setRelease(string r,
                      string yr)
{
   release = r;
   year = yr;
}


void rita::setInput(string file,
                    int    opt)
{
   _script_file = file;
   _opt = opt;
   Load();
}


int rita::run()
{
   int key=0, ret=0;
   string fn="", td="";
   size_t nb_args=0;
   for (;;) {
     if (_cmd->readline(sPrompt+" ")<0)
         continue;
      nb_args = _cmd->getNbArgs();
      key = _cmd->getKW(_rita_kw,_gkw,_data_kw);
      if (key>=200) {
         _data->setDataExt(key);
         continue;
      }
      switch (key) {

         case   0:
            Load();
            break;

         case 201:
            _mesh = new mesh(this,_cmd,_configure);
            _mesh->setVerbose(_verb);
            ret = _mesh->run();
            if (ret) {
               delete _mesh;
               _mesh = nullptr;
            }
            break;

         case   1:
            _eigen = new eigen(this,_cmd,_configure);
            _eigen->run();
            break;

         case   2:
            _optim = new optim(this,_cmd,_configure);
            _optim->run();
            break;

         case   3:
            _approx = new approximation(this,_cmd,_configure);
            _approx->run();
            break;

         case   4:
            _integration = new integration(this,_cmd,_configure);
            if (!_integration->run())
               _integration->go();
            break;

         case   5:
         case   6:
            _ae = new ae(this,_cmd,_configure);
            ret = _ae->run();
            break;

         case   7:
            _ode = new ode(this,_cmd,_configure);
            ret = _ode->run();
            break;

         case   8:
            _pde = new pde(this,_cmd,_configure);
            ret = _pde->run();
            break;

         case   9:
            _ls = new ls(this,_cmd,_configure);
            ret = _ls->run();
            break;

         case  10:
            _solve = new solve(this,_cmd,_configure);
            _solve->setVerbose(_verb);
            ret = _solve->run();
            if (ret>=500)
               return ret;
            break;

         case 100:
            hh->run(0);
            break;

         case 101:
            hh->run(1);
            break;

         case 102:
            getLicense();
            break;

         case 103:
            ret = _configure->run();
            _verb = _configure->getVerbose();
            break;

         case 104:
         case 105:
            msg("","No higher level available");
            break;

         case 106:
            setEcho(nb_args);
            break;

         default:
            DEFAULT_KW
      }
   }
   return 0;
}


void rita::finish()
{
   *ofh << "exit" << endl;
   exit(0);
}


void rita::getLicense()
{
   _cmd->setNbArg(0);
   cout << "Copyright (C) 2024\n";
   cout << "rita is free software: you can redistribute it and/or modify\n";
   cout << "it under the terms of the GNU General Public License as published by\n";
   cout << "the Free Software Foundation, either version 3 of the License, or\n";
   cout << "(at your option) any later version.\n\n";
   cout << "rita is distributed in the hope that it will be useful,\n";
   cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
   cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n";
   cout << "GNU General Public License for more details." << endl;
}


int rita::Load()
{
   int ret=0;
   if (_script_file.size()==0) {
      if (_cmd->setNbArg(1,"Give input script file.")) {
         msg("load>","Missing script file to load.","",1);
         return 1;
      }
      _cmd->get(_script_file);
   }
   if (_in!=nullptr)
      delete _in, _in = nullptr;
   _in = new ifstream(_script_file);
   if (_in->is_open())
      _cmd->setIFStream(_in);
   else {
      msg("load>","Unable to open file: "+_script_file);
      if (_opt==0)
         exit(EXIT_SUCCESS);
      ret = 1;
   }
   _script_file = "";
   _load = true;
   ret = run();
   return ret;
}


void rita::setEcho(size_t nb_args)
{
   static vector<string> kw {"on","off"};
   _cmd->set(kw);
   int k = _cmd->getArg();
   if (nb_args==0 || k==0) {
      _echo = true;
      cout << "Echo is set." << endl;
      *ofh << " echo on" << endl;
   }
   else if (k==1) {
      _echo = false;
      cout << "Echo is unset." << endl;
      *ofh << " echo off" << endl;
   }
   else
      cout << "Echo status is unchanged" << endl;
}


void rita::msg(const string& loc,
               const string& m1,
               const string& m2,
               int           c)
{
   if (c==0) {
      cout << m1 << endl;
      if (m2!="")
         cout << m2 << endl;
   }
   *ofl << "In " + sPrompt+loc + ": " << m1 << endl;
}

} /* namespace RITA */
