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

                        Implementation of class 'optim'

  ==============================================================================*/

#include "optim.h"
#include "configure.h"
#include "cmd.h"
#include "data.h"
#include "rita.h"
#include "calc.h"
#include "defs.h"
#include "helps.h"

namespace RITA {

optim::optim(rita*      r,
             cmd*       command,
             configure* config)
      : _rita(r), _configure(config), _cmd(command)
{
   _verb = _rita->_verb;
   log = true;
   nb_lec = nb_gec = nb_eqc = 0;
   G_ok = H_ok = lp = false;
   penal = 1./OFELI_TOLERANCE;
   _data = _rita->_data;
   solved = 0;
}


optim::~optim()
{
}


int optim::run()
{
   _rita->_analysis_type = analysis_type::OPTIMIZATION;
   size = 1;
   b = 0.0;
   opt_name = "";
   _alg = "gradient";
   int ret=0, ind=0, nn=0, k=0;
   int count_fct=0, count_obj=0, count_grad=0, count_hess=0, count_init=0, count_lp=0;
   int count_lec=0, count_gec=0, count_eqc=0, count_vector=0, penal_ok=0, nb=0, key=0;
   double in=0., xx=0.;
   vector<string> grad, hess, var, le_cons, eq_cons;
   vector<size_t> nbv;
   OFELI::Vect<double> *A_le, *A_ge, *A_eq;
   string str, J, name="", constr="";
   string fn="";
   var_name.push_back("x");
   init.clear();
   grad.clear();
   hess.clear();
   map<int,double> llb, uub;
   static const vector<string> kw {"size","func$tion","obj$ective","lp","grad$ient","hess$ian","low$-bound",
                                   "up$-bound","ge$-constraint","le$-constraint","eq$-constraint",
                                   "penal$ty","var$iable","vect$or","init$ial","algo$rithm"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   for (size_t k=0; k<nb_args; ++k) {

      switch (_cmd->getArgs(nb)) {

         case 0:
            ret = _data->getPar(0,_pr,size);
            if (ret)
               MSGR(_pr,"Illegal parameter value")
            break;

         case 1:
            NO_VALUE_ARG(_pr)
            name = _cmd->string_token(0);
            count_fct++;
            break;

         case 2:
            NO_VALUE_ARG(_pr)
            J = _cmd->string_token(0);
            count_obj++;
            break;

         case 3:
            lp = true;
            break;

         case 4:
            NO_VALUE_ARG(_pr)
            for (int i=0; i<nb; ++i) {
               grad.push_back(_cmd->string_token(i));
               count_grad++;
            }
            break;

         case 5:
            NO_VALUE_ARG(_pr)
            for (int i=0; i<nb; ++i) {
               hess.push_back(_cmd->string_token(i));
               count_hess++;
            }
            break;

         case 6:
            NO_VALUE_ARG(_pr)
            ret = _data->getPar(0,_pr,nn);
            ret  = _data->getPar(1,_pr,xx);
            llb[nn] = xx;
            break;

         case 7:
            NO_VALUE_ARG(_pr)
            ret = _data->getPar(0,_pr,nn);
            ret  = _data->getPar(1,_pr,xx);
            uub[nn] = xx;
            break;

         case 12:
         case 13:
            NO_VALUE_ARG(_pr)
            var_name[0] = _cmd->string_token(0);
            count_vector++;
            break;

         case 14:
            NO_VALUE_ARG(_pr)
            for (int i=0; i<nb; ++i) {
               ret  = _data->getPar(0,_pr,in);
               init.push_back(in);
               count_init++;
            }
            break;

         case 15:
            NO_VALUE_ARG(_pr)
            _alg = _cmd->string_token(0);
            break;

         case 100:
            cout << "Available arguments: " << Opt_hhelp << endl;
            break;

         case 101:
            cout << Opt_Help << endl;
            break;

         default:
            UNKNOWN_ARG(_pr)
      }
   }

   if (nb_args>0) {
      CHK_MSGR(size<=0,_pr,"Illegal size value.")
      CHK_MSGR(!count_obj && !count_fct,_pr,"Missing objective function.")
      if (count_fct) {
         CHK_MSGR(count_vector,_pr,"Function already defined.")
         CHK_MSGR(count_fct && count_obj,_pr,"Function already defined.")
      }
      CHK_MSGR(count_fct>1 || count_obj>1,_pr,"Too many objective functions defined.")
      CHK_MSGR(count_obj && !count_vector,_pr,"Missing a variable name.")
      CHK_MSGR(count_init>size,_pr,"Too many initial guesses given.")
//      if (_data->VectorLabel.count(str)>0 && _data->dn[var_name[0]].active) {
//         _rita->msg(_pr,"Variable name "+var_name+" already used.");
//         return 1;
//      }
      *_rita->ofh << "optimization";
      if (count_fct>0) {
         ind = _data->FctLabel[name];
         CHK_MSGR(ind<=0,_pr,"Non defined function "+name)
         J_Fct = _data->theFct[ind];
//         CHK_MSGR(J_Fct->set(name,_data->theFct[ind]->expr,_data->theFct[ind]->var,1),_pr,"Error in function evaluation: "+J_Fct->getErrorMessage())
         *_rita->ofh << " function=" << name;
      }
      else {
         ind = _data->addVector(var_name[0],0.,size,"",SetCalc::SET);
         opt_var = _data->theVector[ind];
         if (size==1)
            var.push_back(var_name[0]);
         else {
            for (int i=0; i<size; ++i)
               var.push_back(var_name[0]+to_string(i+1));
         }
         opt_var = _data->theVector[ind];
         *_rita->ofh << " var=" << var_name[0];
         k = _data->addFunction(fn,var,J);
         FCT_NOT_DEFINED("",fn)
         FCT_ALREADY_DEFINED("",fn)
         J_Fct = _data->theFct[k];
         if (count_grad) {
            CHK_MSGR(count_grad!=size,_pr,"Illegal number of gradient components given.")
            G_ok = true;
            for (int i=0; i<size; ++i) {
               k = _data->addFunction(fn,var,grad[i]);
               FCT_NOT_DEFINED("",fn)
               FCT_ALREADY_DEFINED("",fn)
               igrad = k;
               G_Fct.push_back(_data->theFct[k]);
            }
            *_rita->ofh << " gradient=" << grad[0];
            for (int i=1; i<size-1; ++i)
               *_rita->ofh << grad[i] << ",";
            *_rita->ofh << grad[size-1];
         }
         if (count_hess) {
            CHK_MSGR(count_hess!=size*size,_pr,"Illegal number of hessian components given.")
            H_ok = true;
            for (int i=0; i<size; ++i) {
               for (int j=0; j<size; ++j) {
                  k = _data->addFunction(fn,var,hess[size*i+j]);
                  FCT_NOT_DEFINED("",fn)
                  FCT_ALREADY_DEFINED("",fn)
                  ihess = k;
                  H_Fct.push_back(_data->theFct[k]);
               }
            }
            *_rita->ofh << " hessian=" << hess[0];
            for (int i=1; i<size*size-1; ++i)
               *_rita->ofh << hess[i] << ",";
            *_rita->ofh << hess[size*size-1];
         }
         for (int i=0; i<size; ++i) {
            lb.push_back(-std::numeric_limits<real_t>::max());
            ub.push_back( std::numeric_limits<real_t>::max());
         }
         for (const auto& v: llb) {
            lb[v.first-1] = v.second;
            *_rita->ofh << " low-bound=" << v.first << "," << v.second;
         }
         for (const auto& v: uub) {
            ub[v.first-1] = v.second;
            *_rita->ofh << " up-bound=" << v.first << "," << v.second;
         }
         Alg = Nopt[_alg];
         *_rita->ofh << " method=" << _alg;
         for (int i=count_init; i<size; ++i)
            init.push_back(0.);
         *_rita->ofh << " init=" << init[0];
         for (int i=1; i<size; ++i)
            *_rita->ofh << "," << init[i];
         for (int i=0; i<size; ++i)
            (*opt_var)[i] = init[i];
         _data->addOpt(this,opt_name);
         if (_verb)
            cout << "Optimization problem " << opt_name << " created." << endl;
         *_rita->ofh << endl;
      }
      log = false;
   }

   else {
      while (1) {
         READLINE
         switch (key) {

            case   0:
               ret = 1;
               CHK_MSGB(size<=0,_pr+"size>","Illegal size value.")
               if (_cmd->get(size))
                  break;
               ret = 0;
               break;

            case   1:
               ret = 1;
               CHK_MSGB(lp,_pr+"function>","Argument not compatible with linear programming.")
               CHK_MSG1B(_cmd->setNbArg(1,"Function name to be given."),_pr+"function>","Missing function name.","",1)
               if (_cmd->get(name))
                  break;
               count_fct++;
               ret = 0;
               break;

            case   2:
               if (lp) {
                  CHK_MSGB(count_lp==size+1,_pr+"objective>","Too many objective function data.")
                  while (!_data->getPar(-1,_pr+"objective>",xx)) {
                     if (count_lp==size)
                        b = xx;
                     else
                        a.push_back(xx);
                     count_lp++;
                  }
               }
               else {
                  ret = 1;
                  CHK_MSG1B(_cmd->setNbArg(1,"Objective function to be given.",1),_pr+"objective>","Missing objective function expression.","",1)
                  if (_cmd->get(J))
                     break;
                  ret = 0;
                  count_obj++;
               }
               ret = 0;
               break;

            case   3:
               ret = 1;
               CHK_MSGB(count_fct || count_grad || count_hess || penal_ok,_pr+"lp>","Linear programming option incompatible with other arguments.")
               ret = 0;
               lp = true;
               count_lec = count_gec = count_eqc = 0;
               break;

            case   4:
               ret = 1;
               CHK_MSGB(lp,_pr+"gradient>","Argument not compatible with linear programming.")
               CHK_MSGB(count_grad==size,_pr+"gradient>","Too many gradient definitions.")
               CHK_MSG1B(_cmd->setNbArg(size,"Partial derivative of objective function to be given."),_pr+"gradient>","Missing partial derivatives of objective function.","",1)
               ret = _cmd->get(str);
               grad.push_back(str);
               count_grad++;
               ret = 0;
               break;

            case   5:
               ret = 1;
               CHK_MSGB(lp,_pr+"hessian>","Argument not compatible with linear programming.")
               CHK_MSGB(count_hess==size*size,_pr+"hessian>","Too many functions defining hessian.")
               CHK_MSG1B(_cmd->setNbArg(size,"Partial derivatives of gradient to be given."),_pr+"hessian>","Missing partial derivatives of gradient.","",1)
               ret = _cmd->get(str);
               hess.push_back(str);
               count_hess++;
               ret = 0;
               break;

            case   6:
               ret = 1;
               CHK_MSGB(lp,_pr+"low-bound>","Argument not compatible with linear programming.")
               CHK_MSG1B(_cmd->setNbArg(2,"Lower value to enforce to be given."),_pr+"low-bound>","Lower value to enforce to be given.","",1)
               ret  = _cmd->get(ind);
               ret += _cmd->get(xx);
               llb[ind] = xx;
               break;

            case   7:
               ret = 1;
               CHK_MSGB(lp,"optimization>up-bound>","Argument not compatible with linear programming.")
               CHK_MSG1B(_cmd->setNbArg(2,"Upper value to enforce to be given."),_pr+"up-bound>","Upper value to enforce to be given.","",1)
               ret = _cmd->get(ind);
               if (!_data->getPar(-1,_pr+"up-bound>",xx))
                  uub[ind] = xx;
               break;

            case   8:
               ret = 1;
               CHK_MSGB(!lp,_pr+"ge-constraint>","This argument is valid for linear programming problems only.")
               ret = 0;
               count_gec = 0;
               while (!_data->getPar(-1,_pr+"ge-constraint>",xx)) {
                  if (count_gec==0)
                     A_ge = new OFELI::Vect<double>(size);
                  if (count_gec==size) {
                     a_ge.push_back(A_ge);
                     b_ge.push_back(xx);
                  }
                  else
                     (*A_ge)[count_gec++] = xx;
               }
               break;

            case   9:
               if (lp) {
                  count_lec = 0;
                  while (!_data->getPar(-1,_pr+"le-constraint>",xx)) {
                     if (count_lec==0)
                        A_le = new OFELI::Vect<double>(size);
                     if (count_lec==size) {
                        a_le.push_back(A_le);
                        b_le.push_back(xx);
                     }
                     else
                        (*A_le)[count_lec++] = xx;
                  }
                  break;
               }
               ret = _cmd->get(str);
               le_cons.push_back(str);
               count_lec++;
               break;

            case  10:
               if (lp) {
                  count_eqc = 0;
                  while (!_data->getPar(-1,_pr+"eq-constraint>",xx)) {
                     if (count_eqc==0)
                        A_eq = new OFELI::Vect<double>(size);
                     if (count_eqc==size) {
                        a_eq.push_back(A_eq);
                        b_eq.push_back(xx);
                     }
                     else
                        (*A_eq)[count_eqc++] = xx;
                  }
                  break;
               }
               ret = _cmd->get(str);
               eq_cons.push_back(str);
               count_eqc++;
               break;

            case  11:
               CHK_MSGB(lp,_pr+"penal>","Argument not compatible with linear programming.")
               if (!_data->getPar(-1,_pr+"penal>",penal))
                  penal_ok++;
               break;
 
            case  12:
            case  13:
               CHK_MSG1B(_cmd->setNbArg(1,"Give name of associated vector."),_pr+"variable>","Missing name of associated vector.","",1)
               _cmd->get(str);
//               if (_data->VectorLabel.count(str)>0 && _data->dn[var_name].active) {
//                  _rita->msg(_pr+"vector>","Variable name "+var_name+" already used.");
//                  break;
//               }
               var_name[0] = str;
               count_vector++;
               break;

            case  14:
               CHK_MSGB(lp,_pr+"initial>","Argument not necessary for linear programming.")
               CHK_MSGB(count_init==size,_pr+"initial>","Too many initial solutions definitions.")
               CHK_MSG1B(_cmd->setNbArg(size,1,"Initial guess to be given.",1),_pr+"initial>","Missing initial guess.","",1)
               if (!_data->getPar(-1,_pr+"initial>",xx))
                  init.push_back(xx), count_init++;
               break;

            case  15:
               CHK_MSG1B(_cmd->setNbArg(1,"Optimization algorithm to be supplied.",1),_pr+"algorithm>","Missing optimization algorithm.","",1)
               ret = _cmd->get(_alg);
               CHK_MSGB(Nopt.find(_alg)==Nopt.end(),_pr+"algorithm>","Unknown or nonimplemented optimization algorithm: "+_alg)
               break;

            case 100:
               cout << "Available commands:\n" << Opt_hhelp << endl;
               break;

            case 101:
               cout << "Available commands:\n" << Opt_HHelp << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               ret = _rita->_configure->run();
               break;

            case 104:
            case 105:
               _cmd->setNbArg(0);
               if (lp) {
                  CHK_MSGR(!count_lp && lp,_pr+"end>","Missing objective function.")
                  CHK_MSGR(!count_vector,_pr+"end>","Missing a variable name.")
                  *_rita->ofh << "optimization\n  size " << size << endl << "  lp" << endl;
                  ind = _data->addVector(var_name[0],0.,size,"",SetCalc::SET);
                  opt_var = _data->theVector[ind];
                  if (size==1)
                     var.push_back(var_name[0]);
                  else {
                     for (int i=0; i<size; ++i)
                        var.push_back(var_name[0]+to_string(i+1));
                  }
                  *_rita->ofh << "  variable " << var_name[0] << endl;
                  *_rita->ofh << "  objective ";
                  for (int i=0; i<size; ++i)
                     *_rita->ofh << a[i] << " ";
                  *_rita->ofh << b << endl;
                  nb_lec = a_le.size(), nb_gec = a_ge.size(), nb_eqc = a_eq.size();
                  for (int i=0; i<nb_eqc; ++i) {
                     *_rita->ofh << "  eq-constraint ";
                     for (int j=0; j<size; ++j) 
                        *_rita->ofh << (*a_eq[i])[j] << " ";
                     *_rita->ofh << b_eq[i] << endl;
                  }
                  for (int i=0; i<nb_lec; ++i) {
                     *_rita->ofh << "  le-constraint ";
                     for (int j=0; j<size; ++j)
                        *_rita->ofh << (*a_le[i])[j] << " ";
                     *_rita->ofh << b_le[i] << endl;
                  }
                  for (int i=0; i<nb_gec; ++i) {
                     *_rita->ofh << "  ge-constraint ";
                     for (int j=0; j<size; ++j)
                        *_rita->ofh << (*a_ge[i])[j] << " ";
                     *_rita->ofh << b_ge[i] << endl;
                  }
               }
               else {
                  CHK_MSGR(!count_obj && !count_fct,_pr+"end>","Missing objective function.")
                  CHK_MSGR(count_fct && count_vector,_pr+"end>","Function already defined.")
                  CHK_MSGR(count_fct && count_obj,_pr+"end>","Function already defined.")
                  CHK_MSGR(count_fct>1 || count_obj>1,_pr+"end>","Too many objective functions defined.")
                  CHK_MSGB(count_obj && !count_vector,_pr+"end>","Missing a variable name.")
                  *_rita->ofh << "optimization\n  size " << size << endl;
                  if (count_fct>0) {
                     ind = _data->FctLabel[name];
                     ret = 1;
                     CHK_MSGB(ind==-1,_pr+"end>","Non defined function "+name)
                     J_Fct = _data->theFct[ind];
                     ret = 0;
                     *_rita->ofh << "  function " << name << endl;
                  }
                  else {
                     ind = _data->addVector(var_name[0],0.,size,"",SetCalc::SET);
                     if (size==1)
                        var.push_back(var_name[0]);
                     else {
                        for (int i=1; i<=size; ++i)
                           var.push_back(var_name[0]+to_string(i));
                     }
                     opt_var = _data->theVector[ind];
                     *_rita->ofh << "  variable " << var_name[0] << endl;
                     *_rita->ofh << "  objective \"" << J << "\"" << endl;
                     int k = _data->addFunction(fn,var,J);
                     FCT_NOT_DEFINED("end>",fn)
                     FCT_ALREADY_DEFINED("end>",fn)
                     J_Fct = _data->theFct[k];
                  }
                  if (count_grad) {
                     CHK_MSGB(count_grad!=size,_pr+"end>","Illegal number of gradient components given.")
                     G_ok = true;
                     for (int i=0; i<size; ++i) {
                        int k = _data->addFunction(fn,var,grad[i]);
                        FCT_NOT_DEFINED("end>",fn)
                        FCT_ALREADY_DEFINED("end>",fn)
                        igrad = k;
                        G_Fct.push_back(_data->theFct[k]);
                     }
                     *_rita->ofh << "  gradient  ";
                     for (int i=0; i<size; ++i)
                        *_rita->ofh << "\"" << grad[i] << "\" ";
                     *_rita->ofh << endl;
                  }
                  nb_lec = count_lec, nb_eqc = count_eqc;
                  for (int i=0; i<nb_lec; ++i) {
                     int k = _data->addFunction(fn,var,le_cons[i]);
                     FCT_NOT_DEFINED("end>",fn)
                     FCT_ALREADY_DEFINED("end>",fn)
                     iincons = k;
                     inC_Fct.push_back(_data->theFct[k]);
                     *_rita->ofh << "  le-constraint \"" << le_cons[i] << "\"" << endl;
                  }
                  for (int i=0; i<nb_eqc; ++i) {
                     int k = _data->addFunction(fn,var,eq_cons[i]);
                     FCT_NOT_DEFINED("end>",fn)
                     FCT_ALREADY_DEFINED("end>",fn)
                     ieqcons = k;
                     eqC_Fct.push_back(_data->theFct[k]);
                     *_rita->ofh << "  eq-constraint \"" << eq_cons[i] << "\"" << endl;
                  }
                  if (penal_ok)
                     *_rita->ofh << "  penalty " << penal << endl;
                  if (count_hess) {
                     log = true;
                     CHK_MSGR(count_hess!=size*size,_pr+"end>","Illegal number of hessian components given.")
                     log = false;
                     H_ok = true;
                     *_rita->ofh << "  hessian ";
                     for (int i=0; i<size; ++i) {
                        for (int j=0; j<size; ++j) {
                           int k = _data->addFunction(fn,var,hess[size*i+j]);
                           FCT_NOT_DEFINED("end>",fn)
                           FCT_ALREADY_DEFINED("end>",fn)
                           ihess = k;
                           H_Fct.push_back(_data->theFct[k]);
                           *_rita->ofh << "\"" << hess[size*i+j] << "\" ";
                        }
                     }
                     *_rita->ofh << endl;
                  }
                  for (int i=0; i<size; ++i) {
                     lb.push_back(-std::numeric_limits<real_t>::max());
                     ub.push_back( std::numeric_limits<real_t>::max());
                  }
                  if (llb.size()) {
                     for (const auto& v: llb) {
                        lb[v.first-1] = v.second;
                        *_rita->ofh << "  low-bound " << v.first << " " << v.second << endl;
                     }
                  }
                  if (uub.size()) {
                     for (const auto& v: uub) {
                        ub[v.first-1] = v.second;
                        *_rita->ofh << "  up-bound " << v.first << " " << v.second << endl;
                     }
                  }
                  Alg = Nopt[_alg];
                  *_rita->ofh << "  algorithm " << _alg << endl;
                  for (int i=count_init; i<size; ++i)
                     init.push_back(0.);
                  *_rita->ofh << "  init  ";
                  for (int i=0; i<size; ++i) {
                     (*opt_var)[i] = init[i];
                     *_rita->ofh << init[i] << "  ";
                  }
               }
               *_rita->ofh << "  end" << endl;
               _data->addOpt(this,opt_name);
               if (_verb)
                  cout << "Optimization problem " << opt_name << " created." << endl;
               log = false;
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


void optim::print(ostream& s) const
{
   cout << "Optimization problem name:        " << opt_name << endl;
   cout << "Number of optimization variables: " << size << endl;
   cout << "Objective function:               " << J_Fct->getExpression() << endl;
   cout << "Number of equality constraints:   " << nb_eqc << endl;   
   cout << "Number of <= constraints:         " << nb_lec << endl;   
   cout << "Number of >= constraints:         " << nb_gec << endl;
   cout << "Linear programming problem ?      " << lp << endl;
   cout << "Solution vector:                  " << var_name[0] << endl;
   if (_data->dn[var_name[0]].active==false)
      s << "WARNING: Vector " << var_name[0] << " has been removed." << endl;
   if (solved) {
      s << "Optimization vector: " << *opt_var;
      return;
   }
   else {
      s << "Optimization problem unsolved !" << endl;
      return;
   }

}


ostream& operator<<(ostream& s, const optim& o)
{
   o.print(s);
   return s;
}

} /* namespace RITA */
