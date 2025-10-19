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

                         Implementation of class ode

  ==============================================================================*/

#include "ode.h"
#include "rita.h"
#include "data.h"
#include "calc.h"
#include "cmd.h"
#include "configure.h"
#include "defs.h"
#include "helps.h"

using std::cout;
using std::exception;

namespace RITA {

ode::ode(rita*      r,
         cmd*       command,
         configure* config)
    : _rita(r), _data(r->_data), _configure(config), _cmd(command)
{
   _verb = _rita->_verb;
   every = 1;
   solved = 0;
   nb_vars = 0;
   err = "";
}


int ode::run()
{
   bool vector_ok=false;
   double xx=0.;
   string an="", str="", Phase="";
   var_name.push_back("y");
   int ret=0, nb=0, count_fct=0, count_vector=0, count_def=0, count_init=0, ind=-1, an_count=0, an_ret=0;
   int key=0, IPhase=1;
   vector<string> Var;
   vector<size_t> nbv;
   init_time=0., time_step=0.1, final_time=1., adapted_time_step = 0;
   scheme = "forward-euler";
   analytic.clear();
   analytic_f.clear();
   err = "";
   log = false;
   size = 1;
   get_err = 0;
   phase.clear();
   iphase.clear();

   vector<string> def, fct_name, var;
   vector<double> init;
   static const vector<string> kw {"size","func$tion","def$inition","var$iable","vect$or","init$ial",
                                   "final$-time","time-step","scheme","phase","analytic","analytic-function","err$or"};

   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();

   for (size_t k=0; k<nb_args; ++k) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << ODE_help << endl;
            return 0;

         case 101:
            cout << ODE_Help << endl;
            return 0;

         case 1:
            NO_VALUE_ARG(_pr)
            fct_name.push_back(_cmd->string_token(0));
            count_fct++, vector_ok = true;
            break;

         case 2:
            NO_VALUE_ARG(_pr)
            def.push_back(_cmd->string_token(0));
            fct_name.push_back("");
            count_def++;
            break;

         case 3:
         case 4:
            NO_VALUE_ARG(_pr)
            var_name[0] = _cmd->string_token(0);
            for (int i=1; i<nb; ++i) {
               var_name.push_back(_cmd->string_token(i));
               nb_vars++;
            }
            vector_ok = true;
            count_vector++;
            break;

         case 5:
            NO_VALUE_ARG(_pr)
            ret = _data->getPar(0,_pr,xx);
            init.push_back(xx);
            count_init++;
            break;

         case 6:
            NO_VALUE_ARG(_pr)
            ret = _data->getPar(0,_pr,final_time);
            break;

         case 7:
            NO_VALUE_ARG(_pr)
            ret = _data->getPar(0,_pr,time_step);
            break;

         case 8:
            NO_VALUE_ARG(_pr)
            scheme = _cmd->string_token(0);
            break;

         case 10:
            if (int(analytic.size())==1) {
               an_ret = 1;
               continue;
            }
            NO_VALUE_ARG(_pr)
            an = _cmd->string_token(0);
            analytic.push_back(an);
            get_err = 1;
            an_count++;
            break;

         case 11:
            if (int(analytic_f.size())==1) {
               an_ret = 1;
               continue;
            }
            NO_VALUE_ARG(_pr)
            an = _cmd->string_token(0);
            analytic_f.push_back(an);
            get_err = 2;
            an_count++;
            break;

         case 12:
            NO_VALUE_ARG(_pr)
            err = _cmd->string_token(0);
            break;

         default:
            UNKNOWN_ARG(_pr)
      }
   }

   if (nb_args>0) {
      CHK_MSGR(count_fct && count_vector,_pr,"Function already defined.")
      CHK_MSGR(count_fct && count_def,_pr,"Function already defined.")
      CHK_MSGR(count_fct>1 || count_def>1,_pr,"Number of function names is larger than system size.")
      CHK_MSGR(_data->nb_ode>0 && count_vector==0,_pr,"No variable defined as unknown for new ode system.")
      CHK_MSGR(!vector_ok,_pr,"Missing a variable name")
      CHK_MSGR(count_init>1,_pr,"Number of initial conditions is larger than system size.")
      CHK_MSGR(an_ret,_pr,"Too many analytic solutions given.")
      if (count_init<1)
         init.push_back(0.);
      CHK_MSGR(an_count>1,_pr,"Analytic solutions can be given only once.")
      *_rita->ofh << "ode";
      isSet = false;
      name = "";
      iFct.clear();
      isFct = false;
      if (count_fct>0) {
         isFct = true;
         CHK_MSGR(count_fct<=0,_pr,"Number of function names is lower than system size.")
         int n = _data->FctLabel[fct_name[0]];
         CHK_MSGR(n==0,_pr,"Undefined function "+fct_name[0])
         iFct.push_back(n);
         *_rita->ofh << " function=" << fct_name[0];
      }
      else {
         *_rita->ofh << " var=" << var_name[0];
         Var.clear();
         Var.push_back("t");
         Var.push_back(var_name[0]);
         nbv.clear();
         nbv.push_back(1), nbv.push_back(1);
         int k = _data->addFunction(fct_name[0],Var,def[0],nbv);
         FCT_NOT_DEFINED("",fct_name[0])
         FCT_ALREADY_DEFINED("",fct_name[0])
         iFct.push_back(k);
         *_rita->ofh << " definition=" << def[0];
      }
      isSet = true;
      ivect.push_back(_data->addVector(var_name[0],0.,1,""));
      *_rita->ofh << " init=" << init[0];
      Scheme = _sch[scheme];
      *_rita->ofh << " scheme=" << scheme;
      *_rita->ofh << " time-step=" << time_step << " final-time=" << final_time;
      if (get_err==1) {
         *_rita->ofh << " analytic=" << analytic[0];
         string na="";
         Var.clear();
         Var.push_back("t");
         nbv.clear();
         nbv.push_back(1);
         int k = _data->addFunction(na,Var,analytic[0],nbv);
         FCT_NOT_DEFINED("",fct_name[0])
         FCT_ALREADY_DEFINED("",fct_name[0])
         SolFct.push_back(_data->theFct[k]);
      }
      if (get_err==2) {
         *_rita->ofh << " analytic-function=";
         int k = _data->FctLabel[analytic_f[0]];
         CHK_MSGR(k==0,_pr,"Undefined function "+analytic_f[0])
         *_rita->ofh << analytic_f[0];
         SolFct.push_back(_data->theFct[k]);
      }
      if (err!="")
         *_rita->ofh << " err=" << err;
      *_rita->ofh << endl;
      type = DataType::ODE;
      _data->addODE(this,name);
      if (_verb)
         cout << "Ordinary Differential Equation " << name << " created." << endl;
   }

   else {
      *_rita->ofh << "ode" << endl;
      for (;;) {
         READLINE
         switch (key) {

            case   0:
               ret = _cmd->setNbArg(1,"Size of differential system to be given.");
               CHK_MSG1B(ret,_pr+"size>","Missing system size.","",1)
               ret = _cmd->get(size);
               if (!ret)
                  *_rita->ofh << "  size " << size << endl;
               break;

            case   1:
               if (count_def>0)
                  ret = 1;
               CHK_MSGB(ret,_pr+"function>","Function already defined by an expression.")
               CHK_MSGB(count_fct==size,_pr+"function>","Too many functions defining ODE.")
               ret = _cmd->setNbArg(1,"Function F to define equation y'(t) = F(t,y(t)) to be given.");
               CHK_MSG1B(!ret,_pr+"function>","Missing function expression.","",1)
               CHK_MSGB(_data->FctLabel.count(str)==0,_pr+"function>","Non defined function "+str)
               CHK_MSGB(_data->dn[str].active==false,_pr+"function>","Function "+str+" has been removed")
               ret = _cmd->get(fct_name[count_fct]);
               if (!ret) {
                  *_rita->ofh << "  function " << fct_name[count_fct++] << endl;
                  vector_ok = true;
                  count_fct++;
               }
               break;

            case   2:
               CHK_MSGB(count_fct>0,_pr+"definition>","Function already defined by its name.")
               CHK_MSGB(count_def==size,_pr+"definition>","Too many functions defining ODE.")
               CHK_MSG1B(_cmd->setNbArg(1,"Function F to define equation y'(t) = F(t,y(t)) to be given."),_pr+"definition>","Missing function expression","",1)
               ret = _cmd->get(str);
               if (!ret) {
                  def.push_back(str);
                  *_rita->ofh << "  definition \"" << str << "\"" << endl;
                  fct_name.push_back("");
                  count_def++;
               }
               break;

            case   3:
            case   4:
               if (log)
                  ret = 1;
               CHK_MSGB(ret,_pr+"variable>","ode must be set first.")
               CHK_MSG1B(_cmd->setNbArg(1,"Give name of associated variable (vector).",1),_pr+"variable>","Missing name of associated variable (vector).","",1)
               nb_vars = _cmd->getNbArgs();
               for (size_t i=0; i<nb_vars; ++i) {
                  ret = _cmd->get(str);
                  CHK_MSGB(ret,_pr+"variable>","Unknown variable "+str)
                  vector_ok = true, count_vector++;
                  if (i==0) {
                     var_name[0] = str;
                     *_rita->ofh << "  variable " << str;
                  }
                  else {
                     var_name.push_back(str);
                     *_rita->ofh << " " << str;
                  }
                  *_rita->ofh << endl;
               }
               break;

            case   5:
               CHK_MSG1B(_cmd->setNbArg(size,"Initial conditions to be given."),_pr+"initial>","Missing initial conditions.","",1)
               init.resize(size);
               count_init++;
               for (int i=0; i<size; ++i)
                  ret += _cmd->get(init[i]);
               CHK_MSGB(ret,_pr+"initial>","Error in initial data.")
               *_rita->ofh << "  initial  ";
               for (const auto& v: init)
                  *_rita->ofh << v << "  ";
               *_rita->ofh << endl;
               break;

            case   6:
               CHK_MSG1B(_cmd->setNbArg(1,"Final time to be given."),_pr+"final-time>","Missing final time value.","",1)
               ret = _cmd->get(final_time);
               if (!ret)
                  *_rita->ofh << "  final-time " << final_time << endl;
               break;

            case   7:
               CHK_MSG1B(_cmd->setNbArg(1,"Missing time step."),_pr+"time-step>","Missing time step value.","",1)
               ret = _cmd->get(time_step);
               if (!ret)
                  *_rita->ofh << "  time-step " << time_step << endl;
               break;

            case   8:
               CHK_MSG1B(_cmd->setNbArg(1,"Time integration scheme to be given."),_pr+"scheme>","Missing time integration scheme.","",1)
               ret = _cmd->get(scheme);
               if (!ret)
                  *_rita->ofh << "  scheme " << scheme << endl;
               break;

            case   9:
               if (size==1) {
                  CHK_MSG1B(_cmd->setNbArg(1,"Name of history vector to contain phase."),_pr+"phase>","Missing Argument","",1)
                  ret = _cmd->get(Phase);
                  if (!ret) {
                     phase.push_back(Phase);
                     iphase.push_back(1);
                     *_rita->ofh << "  phase " << Phase << endl;
                  }
               }
               else {
                  CHK_MSG1B(_cmd->setNbArg(2,"Name of history vector to contain phase and solution component."),_pr+"phase>","Missing Argument","",1)
                  ret  = _cmd->get(Phase);
                  ret += _cmd->get(IPhase);
                  if (!ret) {
                     phase.push_back(Phase);
                     iphase.push_back(IPhase);
                     _data->setHistory(Phase,1);
                     *_rita->ofh << "  phase " << Phase << " " << IPhase << endl;
                  }
               }
               break;

            case  10:
               CHK_MSG1B(_cmd->setNbArg(1,"Expression of the analytical solution of the ODE."),_pr+"analytic>","Missing analytical solution expression.","",1)
               ret = _cmd->get(an);
               if (!ret) {
                  CHK_MSGR(int(analytic.size())==size,_pr+"analytic>","Too many analytic solutions given.")
                  CHK_MSGR(an_count==size,_pr+"analytic>","Too many analytic solutions given.")
                  analytic.push_back(an);
                  get_err = 1;
                  an_count++;
                  *_rita->ofh << "  analytic " << an << endl;
               }
               break;

            case  11:
               CHK_MSG1B(_cmd->setNbArg(1,"Function defininig the analytical solution of the ODE."),_pr+"analytic-function>","Missing analytical solution function.","",1)
               ret = _cmd->get(an);
               if (!ret) {
                  CHK_MSGR(int(analytic_f.size())==size,_pr+"analytic-function>","Too many analytic solutions given.")
                  CHK_MSGR(an_count==size,_pr+"analytic-function>","Too many analytic solutions given.")
                  analytic_f.push_back(an);
                  get_err = 2;
                  an_count++;
                  *_rita->ofh << "  analytic-function " << an << endl;
               }
               break;

            case  12:
               if (!_cmd->get(err))
                  *_rita->ofh << "  err " << err << endl;
               break;

            case 100:
               cout << "Available arguments:\n" << ODE_help << endl;
               break;

            case 101:
               cout << ODE_HHelp << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
                _configure->run();
               break;

            case 104:
            case 105:
cout<<"**1**"<<endl;
               if (ret) {
                  NO_ODE
                  *_rita->ofh << "  end" << endl;
                  return 1;
               }
cout<<"**2**"<<endl;
               if (count_fct==0 && count_def==0) {
                  NO_ODE
                  *_rita->ofh << "  end" << endl;
                  return 0;
               }
cout<<"**3**"<<endl;
               if ((count_fct>0 && count_fct<size) || (count_fct==0 && count_def<size)) {
                  *_rita->ofh << "  end" << endl;
                  MSGB(_pr+"end>","Insufficient number of functions defining system.")
               }
cout<<"**4**"<<endl;
               if (!vector_ok) {
                  _rita->msg(_pr+"end>","No variable defined for ode system.");
                  *_rita->ofh << "  end" << endl;
                  break;
               }
cout<<"**5**"<<endl;
               if (count_fct && count_vector)
                  MSGR("ode>","Function already defined.")
               CHK_MSGR(an_count>1,_pr,"Commands analytic and analytic-function cannot be given simultaneously.")
               *_rita->ofh << "  end" << endl;
               isSet = false;
               name = "";
               if (count_fct>0) {
                  isFct = true;
                  for (int i=0; i<size; ++i) {
                     int n = _data->FctLabel[fct_name[i]];
                     CHK_MSGR(n==0,_pr,"Undefined function "+fct_name[i])
                     iFct.push_back(n);
                  }
               }
               else {
                  isFct = false;
                  for (int i=0; i<size; ++i) {
                     Var.clear();
                     Var.push_back("t");
                     for (size_t j=0; j<nb_vars; ++j)
                        Var.push_back(var_name[j]);
                     nbv.clear();
                     nbv.push_back(1);
                     nbv.push_back(size);
                     int k = _data->addFunction(fct_name[i],Var,def[i],nbv);
                     FCT_NOT_DEFINED("end>",fct_name[i])
                     FCT_ALREADY_DEFINED("end>",fct_name[i])
                     iFct.push_back(k);
                  }
                  ind_fct = ind;
               }
               if (!count_init) {
                  init.resize(size);
                  for (int i=0; i<size; ++i)
                     init.push_back(0.);
               }
               ind_fct = ind;
               isFct = count_fct;
               for (size_t i=0; i<nb_vars; ++i)
                  ivect.push_back(_data->addVector(var_name[i],0.,size,""));
               _data->theVector[ivect[0]]->resize(size);
               Scheme = _sch[scheme];
               for (size_t i=0; i<phase.size(); ++i) {
                  if (phase[i]!="")
                     _data->addVector(phase[i],0,size,"");
               }
               for (int j=0; j<size; ++j)
                  _data->theVector[ivect[0]]->at(j) = init[j];
               isSet = true;
               log = false;
               isFct = false;
               if (count_fct)
                  isFct = true;
               if (get_err==1) {
                  for (int i=0; i<size; ++i) {
                     string nm="";
                     Var.clear();
                     Var.push_back("t");
                     for (size_t j=1; j<nb_vars; ++j)
                        Var.push_back(var_name[j]);
                     nbv.clear();
                     nbv.push_back(size);
                     int k = _data->addFunction(nm,Var,analytic[i],nbv);
                     FCT_NOT_DEFINED("end>",nm)
                     FCT_ALREADY_DEFINED("end>",nm)
                     SolFct.push_back(_data->theFct[k]);
                  }
               }
               if (get_err==2) {
                  for (int i=0; i<size; ++i) {
                     int k = _data->FctLabel[analytic_f[i]];
                     CHK_MSGR(k==0,_pr+"end>","Function "+analytic_f[i]+" already defined.")
                     SolFct.push_back(_data->theFct[k]);
                  }
               }
               _data->addODE(this,name);
               if (_verb)
                  cout << "Ordinary Differential Equation " << name << " created." << endl;
               return 0;

            case 106:
               _rita->setEcho(nb_args);
               break;

            case -2:
            case -3:
            case -4:
               break;

            default:
               DEFAULT_KW(_rita)
         }
      }
   }
   return ret;
}


void ode::print(ostream& s) const
{
   s << "Ordinary Differential Equation Name: " << name << endl;
   s << "Size: " << size << endl;
   if (_data->dn[var_name[0]].active==false)
      s << "WARNING: Vector " << var_name[0] << " has been removed." << endl;
   s << "Name of variable(s): ";
   s << var_name[0];
   for (size_t i=1; i<nb_vars; ++i)
      s << ", " << var_name[i];
   s << endl;
   s << "Time integration scheme: " << scheme << endl;
   s << "Initial time: " << init_time << endl;
   s << "Final time: " << final_time << endl;
   s << "Time step: " << time_step << endl;
   s << "ODE defined by function(s):\n";
   for (int i=0; i<size; ++i)
      s << _data->theFct[iFct[i]]->getExpression() << endl;
   if (nls!="")
      s << "Nonlinear iteration procedure: " << nls << endl;
   if (solved==0) {
      s << "ODE unsolved !" << endl;
      return;
   }
   s << "Solution of the ordinary differential equation (system):\n" << *_data->theVector[ivect[0]];
}


ostream& operator<<(ostream& s, const ode& e)
{
   e.print(s);
   return s;
}

} /* namespace RITA */