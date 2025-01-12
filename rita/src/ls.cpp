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

                         Implementation of class ls

  ==============================================================================*/

#include "ls.h"
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

ls::ls(rita*      r,
       cmd*       command,
       configure* config)
    : _rita(r), _data(r->_data), _configure(config), _cmd(command)
{
   _verb = _rita->_verb;
   log = true;
   solved = 0;
}


ls::~ls()
{
}


int ls::run()
{
   string matrix_name="", rhs_name="", sol_name="", init_name="";
   solver_name="direct";
   prec_name="diag", init_name="";
   bool prec_ok=false;
   int nb=0;
   const static vector<string> kw {"matrix","rhs","solu$tion","solv$er","prec$onditioner","init$ial"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(_pr)
   for (size_t k=0; k<nb_args; ++k) {

      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << LS_help << endl;
            return 0;

         case 101:
            cout << LS_Help << endl;
            return 0;
   
         case 0:
            NO_VALUE_ARG(_pr)
            matrix_name = _cmd->string_token(0);
            break;

         case 1:
            NO_VALUE_ARG(_pr)
            rhs_name = _cmd->string_token(0);
            break;

         case 2:
            NO_VALUE_ARG(_pr)
            sol_name = _cmd->string_token(0);
            break;

         case 3:
            MSGR(_pr,"Argument 'solver' not implemented yet")
            NO_VALUE_ARG(_pr)
            solver_name = _cmd->string_token(0);
            break;

         case 4:
            MSGR(_pr,"Argument 'preconditioner' not implemented yet")
            NO_VALUE_ARG(_pr)
            prec_name = _cmd->string_token(0);
            break;

         case 5:
            NO_VALUE_ARG(_pr)
            init_name = _cmd->string_token(0);
            break;

         default:
            UNKNOWN_ARG(_pr)
      }
   }

   if (nb_args>0) {
      CHK_MSGR(matrix_name=="",_pr,"No matrix given.")
      CHK_MSGR(rhs_name=="",_pr,"No right-hand side vector given.")
      CHK_MSGR(prec_ok && solver_name=="direct",_pr,"No preconditioner needed for a direct solver.")
      lls = _data->Ls[solver_name];
      CHK_MSGR(lls==0,_pr,"Unknown or non implemented solver "+solver_name)
      pprec = _data->Pr[prec_name];
      CHK_MSGR(pprec==0,_pr,"Unknown or non implemented preconditioner "+prec_name)
      string storage = "dense";
      if (solver_name!="direct")
         storage = "sparse";
      iMat = _data->MatrixLabel[matrix_name];
      CHK_MSGR(iMat==0,_pr,"Undefined reference to matrix "+matrix_name)
      _Mat = _data->theMatrix[iMat];
      int nr=_Mat->getNbRows();
      CHK_MSGR(nr!=int(_Mat->getNbColumns()),_pr,"Matrix "+matrix_name+" is not a square matrix")
      iRHS = _data->VectorLabel[rhs_name];
      CHK_MSGR(iRHS==0,_pr,"Undefined reference to right-hand side vector "+rhs_name)
      CHK_MSGR(nr!=int(_data->theVector[iRHS]->size()),_pr,"Matrix "+matrix_name+" and right-hand side "
               +rhs_name+" have incompatible sizes.")
      *_rita->ofh << "ls matrix=" << matrix_name << " rhs=" << rhs_name;
      iSol = iRHS;
      if (sol_name!="") {
         iSol = _data->addVector(sol_name,0.,nr,"",0);
         *_rita->ofh << " sol=" << sol_name;
      }
      iInit = iSol;
      if (init_name!="") {
         iInit = _data->addVector(init_name,0.,nr,"",0);
         *_rita->ofh << " initial=" << init_name;
      }
      name = "";
      _data->addLS(this,name);
   }
   return 0;
}


void ls::print(ostream& s) const
{
   s << "Linear System Name:     " << name << endl;
   s << "Matrix:                 " << _data->MatrixName[iMat] << endl;
   s << "Right-hand Side vector: " << _data->VectorName[iRHS] << endl;
   s << "Solution vector:        " << _data->VectorName[iSol] << endl;
   if (solver_name!="direct")
      s << "Initialization vector:  " << _data->VectorName[iInit] << endl;
   s << "Solver name:            " << solver_name << endl;
   s << "Preconditioner name:    " << prec_name << endl;
   if (_data->dn[_data->MatrixName[iMat]].active==false)
      s << "WARNING: Vector " << _data->MatrixName[iMat] << " has been removed." << endl;
   if (_data->dn[_data->VectorName[iRHS]].active==false)
      s << "WARNING: Vector " << _data->VectorName[iRHS] << " has been removed." << endl;
   if (solved==0) {
      cout << "Linear system unsolved !" << endl;
      return;
   }
}


ostream& operator<<(ostream& s, const ls& e)
{
   e.print(s);
   return s;
}

} /* namespace RITA */
