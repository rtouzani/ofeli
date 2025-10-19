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

                         Implementation of class ae

  ==============================================================================*/

#include "ae.h"
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

ae::ae(rita*      r,
       cmd*       command,
       configure* config)
   : _rita(r), _data(r->_data), _configure(config), _cmd(command)
{
   _verb = _rita->_verb;
   log = true;
   nb_vars = 0;
   solved = 0;
}


int ae::run()
{
   _pr = _PR + ">algebraic>";
   name = "";
   string str="", fct_name="";
   double xx=0.;
   nls = "newton";
   nnls = _data->NLs[nls];
   int var_ok=0, k=0, nb=0, ind=-1, ret=0, key=0;
   int count_def=0, count_J=0, count_init=0;
   size = 1;
   vector<string> def, fname, var;
   vector<double> init;
   vector<size_t> nbv;
   isSet = false;
   J.setSize(1,1);
   isFct = 0;
   log = false;

   const static vector<string> kw {"size","func$tion","def$inition","jacob$ian",
                                   "init$","var$iable","vect$or","nls"};

   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   for (size_t k=0; k<nb_args; ++k) {

      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << AE_help << endl;
            RET(0)

         case 101:
            cout << AE_Help << endl;
            RET(0)

         case 1:
            NO_VALUE_ARG(_pr)
            fname.push_back(_cmd->string_token(0));
            isFct++;
            break;

         case 2:
            NO_VALUE_ARG(_pr)
            def.push_back(_cmd->string_token(0));
            fname.push_back(" ");
            count_def++;
            break;

         case 4:
            NO_VALUE_ARG(_pr)
            ret = _data->getPar(0,_pr,xx);
            init.push_back(xx);
            count_init++;
            break;

         case 5:
         case 6:
            NO_VALUE_ARG(_pr)
            CHK_MSG1R(var_ok,_pr,"Variable(s) already defined "+_cmd->getToken(),"",0)
            for (int i=0; i<nb; ++i) {
               var_name.push_back(_cmd->string_token(i));
               nb_vars++;
            }
            var_ok++;
            break;

         case 7:
            NO_VALUE_ARG(_pr)
            nls = _cmd->string_token(0);
            break;

         default:
            UNKNOWN_ARG(_pr)
      }
   }

   if (nb_args>0) {
      CHK_MSGR(isFct && var_ok,_pr,"Function already defined.")
      CHK_MSGR(var_ok>1,_pr,"Only one variable must be defined for an algebraic system.")
      CHK_MSGR(!isFct && !count_def,_pr,"No algebraic function defined.")
      CHK_MSGR(isFct && count_def,_pr,"Function already defined.")
      CHK_MSGR(isFct>size || count_def>size,_pr,"Number of function names is larger than system size.")
      CHK_MSGR(!var_ok && !isFct,_pr,"Missing variable name")
      CHK_MSGR(count_init>size,_pr,"Number of initial guesses is larger than system size.")
      CHK_MSGR(isFct>0 && count_def<size-1,_pr,"Number of function definitions is larger than system size.")
      if (count_init<1) {
         for (int i=count_init; i<size; ++i)
            init.push_back(0.);
      }
      *_rita->ofh << "algebraic";
      if (isFct) {
         k = _data->FctLabel[fname[0]];
         CHK_MSGR(k==0,_pr,"Non defined function "+fname[0])
         iFct.push_back(k);
         *_rita->ofh << " function=" << fname[0];
         nb_vars = _data->theFct[k]->nb_var;
         for (size_t i=0; i<nb_vars; ++i)
            var_name.push_back(_data->theFct[k]->var[i]);
      }
      else {
         *_rita->ofh << " var=" << var_name[0];
         for (size_t i=1; i<var_name.size(); ++i)
            *_rita->ofh << "," << var_name[i];
         fn = var_name[0];
         var.clear();
         var.push_back(var_name[0]);
         nbv.resize(1);
         nbv[0] = size;
         k = _data->addFunction(fct_name,var,def[0],nbv);
         FCT_NOT_DEFINED("",fct_name);
         FCT_ALREADY_DEFINED("",fct_name);
         _data->theFct[k]->setExpr(def[0]);
         iFct.push_back(k);
         ind_fct = ind;
         *_rita->ofh << " definition=" << def[0];
      }
      for (size_t i=0; i<nb_vars; ++i)
         ivect.push_back(_data->addVector(var_name[i],0.,size,""));
      J.setSize(1,1);
      isSet = true;
      log = false;
      for (const auto& v: init)
         *_rita->ofh << " init=" << v;
      for (size_t i=0; i<init.size(); ++i)
         (*_data->theVector[ivect[0]])[i] = init[i];
      *_rita->ofh << " nls=" << nls << endl;
      nnls = _data->NLs[nls];
      _data->addAE(this,name);
      if (_verb)
         cout << "Algebraic Equation " << name << " created." << endl;
   }

   else {
      *_rita->ofh << "algebraic" << endl;
      for (;;) {
         READLINE
         switch (key) {

            case   0:
               CHK_MSG1B(_cmd->setNbArg(1,"Size of algebraic system to be given."),_pr+"size>","Missing system size.","",1)
               if (!_cmd->get(size))
                  *_rita->ofh << "  size " << size << endl;
               break;

            case   1:
               CHK_MSGR(isFct==size,_pr+"function>","Number of functions is larger than system size.")
               CHK_MSG1B(_cmd->setNbArg(1,"Name of function to be given."),_pr+"function>","Missing function name","",1)
               str = _cmd->string_token();
               CHK_MSGB(_data->FctLabel.count(str)==0,_pr+"function>","Non defined function "+str)
               CHK_MSGB(_data->dn[str].active==false,_pr+"function>","Function "+str+" has been removed")
               if (!ret) {
                  *_rita->ofh << "  function " << str << endl;
                  fname[isFct++] = str;
                  var_ok++;
                  k = _data->FctLabel[str];
                  iFct.push_back(k);
                  nb_vars = _data->theFct[k]->nb_var;
                  for (size_t i=0; i<nb_vars; ++i)
                     var_name.push_back(_data->theFct[k]->var[i]);
               }
               break;

            case   2:
               CHK_MSGB(isFct>0,_pr+"definition>","Function already defined by its name.")
               CHK_MSGB(count_def==size,_pr+"definition>","Too many functions defining system.")
               CHK_MSG1B(_cmd->setNbArg(1,"Function F to define equation F(x)=0 to be given."),
                         _pr+"definition>","Missing function expression","",1)
               CHK_MSGB(_cmd->get(str),_pr+"definition>","Error in equation definition.")
               def.push_back(str);
               count_def++;
               *_rita->ofh << "  definition " << str << endl;
               break;

            case   3:
               CHK_MSGB(count_J==size,_pr+"jacobian>","Too many functions defining jacobian.")
               CHK_MSG1B(_cmd->setNbArg(size,"Partial derivatives of function defining equation to be given."),
                         _pr+"jacobian>","Missing partial derivatives of function defining equation.","",1)
               for (int i=1; i<=size; ++i)
                  ret = _cmd->get(J(count_J+1,i));
               if (ret==0) {
                  *_rita->ofh << "  jacobian " << endl;
                  for (int i=1; i<=size; ++i)
                     *_rita->ofh << J(count_J+1,i) << "  ";
                  count_J++;
               }
               break;

            case   4:
               CHK_MSG1B(_cmd->setNbArg(size,"Initial guesses to be given."),_pr+"initial>","Missing initial guesses.","",1)
               ret = 0;
               init.resize(size);
               count_init++;
               for (int i=0; i<size; ++i)
                  ret += _cmd->get(init[i]);
               CHK_MSGB(ret,_pr+"initial>","Error in initial guess data.")
               *_rita->ofh << "  initial  ";
               for (const auto& v: init)
                  *_rita->ofh << v << "  ";
               *_rita->ofh << endl;
               break;

            case   5:
            case   6:
               CHK_MSGB(log,_pr+"variable>","Equation must be defined first.")
               CHK_MSG1B(_cmd->setNbArg(1,"Give name of associated vector.",1),_pr+"variable>","Missing name of associated variable.","",1)
               nb_vars = _cmd->getNbArgs();
               for (size_t i=0; i<nb_vars; ++i) {
                  ret = _cmd->get(str);
                  CHK_MSGB(ret,_pr+"variable>","Unknown variable "+str)
                  var_ok = true;
                  if (i==0) {
                     var_name.push_back(str);
                     *_rita->ofh << "  variable " << str;
                  }
                  else {
                     var_name.push_back(str);
                     *_rita->ofh << " " << str;
                  }
                  *_rita->ofh << endl;
               }
               break;

            case   7:
               CHK_MSGB(log,_pr+"nls>","equation must be set first.")
               CHK_MSG1B(_cmd->setNbArg(1,"Nonlinear solver to be supplied.",1),_pr+"nls>","Missing nonlinear solver data.","",1)
               ret = _cmd->get(nls);
               if (!ret)
                  *_rita->ofh << "  nls " << nls << endl;
               break;

            case 100:
               cout << "Available commands:\n" << AE_hhelp << endl;
               break;

            case 101:
               cout << AE_HHelp << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               _configure->run();
               break;

            case 104:
            case 105:
               _cmd->setNbArg(0);
               CHK_MSGR(!var_ok || (!isFct && !count_def),_pr+"end>","Algebraic equation incompletely defined.")
               CHK_MSGB((isFct && isFct<size) || (!isFct && count_def<size),_pr+"end>","Insufficient number of functions defining system.")
               CHK_MSGB(isFct && var_ok,_pr+"end>","Function already defined.")
               *_rita->ofh << "  end" << endl;
               if (!count_init) {
                  init.resize(size);
                  for (int i=0; i<size; ++i)
                     init[i] = 0.;
               }
               if (!count_J) {
                  J.setSize(size,size);
                  for (int i=1; i<=size; ++i)
                     J(i,i) = 1.;
               }
               ind_fct = ind;
               ivect.push_back(_data->addVector(var_name[0],0.,size,"",SetCalc::SET));
               for (size_t i=1; i<nb_vars; ++i) {
                  int k = _data->VectorLabel[var_name[i]];
                  CHK_MSGR(k==0,_pr+"end>","Variable "+var_name[i]+" undefined")
                  ivect.push_back(k);
               }
               fn = var_name[0];
               for (int j=0; j<size; ++j)
                  (*_data->theVector[ivect[0]])[j] = init[j];
               if (!isFct) {
                  for (int j=0; j<size; ++j) {
                     fct_name = "";
                     var.clear();
                     var.push_back(var_name[0]);
                     nbv.resize(1);
                     nbv[0] = size;
                     int k = _data->addFunction(fct_name,var,def[j],nbv);
                     FCT_NOT_DEFINED("end>",fct_name);
                     FCT_ALREADY_DEFINED("end>",fct_name);
                     _data->theFct[k]->setVar(var_name[0],size);
                     iFct.push_back(k);
                  }
               }
               isSet = true;
               log = false;
               _data->addAE(this,name);
               if (_verb)
                  cout << "Algebraic Equation " << name << " created." << endl;
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
   RET(0)
}


void ae::print(ostream& s) const
{
   s << "Algebraic Equation Name: " << name << endl;
   if (_data->dn[var_name[0]].active==false)
      s << "WARNING: Vector " << var_name[0] << " has been removed." << endl;
   if (size>1)
      s << "Size: " << size << endl;
   if (isFct) {
      if (size==1)
         cout << "Equation defined by function: " << _data->theFct[iFct[0]]->name << endl;
      else {
         for (int j=0; j<size; ++j)
            s << "Equation: " << j+1 << ", defined by function: " << _data->theFct[iFct[j]]->name;
      }
      cout << "\nUnknown variable: " << var_name[0] << endl;
   }
   else {
      if (size==1) {
         s << "Equation defined by: " << _data->theFct[iFct[0]]->getExpression() << "\n  Variable(s):      " << var_name[0];
         for (size_t i=1; i<var_name.size(); ++i)
            s << "," << var_name[i];
         s << ". Unknown variable: " << var_name[0] << endl;
      }
      else {
         for (int j=0; j<size; ++j)
            s << "Equation: " << j+1 << ", defined by: " << _data->theFct[iFct[j]]->getExpression() << endl;
         s << "Variable(s):      " << var_name[0];
         for (size_t i=1; i<var_name.size(); ++i)
            s << "," << var_name[i];
         cout << ". Unknown variable: " << var_name[0] << endl;
      }
   }
   s << "Nonlinear iteration procedure: " << nls << endl;
   if (solved==0) {
      s << "Algebraic equation unsolved !" << endl;
      return;
   }
   s << "Solution of the Algebraic equation:\n" << *_data->theVector[ivect[0]];
}


ostream& operator<<(ostream& s, const ae& e)
{
   e.print(s);
   return s;
}

} /* namespace RITA */