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

                      Implementation of class 'eigen'

  ==============================================================================*/

#include "data.h"
#include "configure.h"
#include "eigen.h"
#include "linear_algebra/Matrix_impl.h"
#include "defs.h"
#include "helps.h"

namespace RITA {

eigen::eigen(rita*      r,
             cmd*       command,
             configure* config)
      : _rita(r), _configure(config), _cmd(command)
{
   _verb = _rita->_verb;
   nb_eigv = 0;
   eig_vec = false;
   verbose = 0;
   _data = _rita->_data;
   solved = 0;
   log = true;
}


eigen::~eigen()
{
}


int eigen::run()
{
   _rita->_analysis_type = analysis_type::EIGEN;
   string mat_name="", method="qr", name="";
   evect = "";
   symm = eig_vec = false;
   nb_eigv = 0;
   int nb=0;

   static const vector<string> kw {"matrix","symm$etric","method","nb","eigv","evect"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(_pr)
   for (size_t k=0; k<nb_args; ++k) {

      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Eig_help << endl;
            return 0;

         case 101:
            cout << Eig_Help << endl;
            return 0;

         case 0:
            NO_VALUE_ARG(_pr)
            mat_name = _cmd->string_token(0);
            break;

         case 1:
            symm = true;
            break;

         case 2:
            NO_VALUE_ARG(_pr)
            method = _cmd->string_token(0);
            break;

         case 3:
            NO_VALUE_ARG(_pr)
            CHK_MSGR(_data->getPar(0,_pr,nb_eigv),_pr,"Illegal parameter value")
            break;

         case 4:
            eig_vec = true;
            break;

         case 5:
            NO_VALUE_ARG(_pr)
            evect = _cmd->string_token(0);
            break;

         default:
            log = true;
            UNKNOWN_ARG(_pr)
      }
   }

   if (nb_args>0) {
      CHK_MSGR(mat_name=="",_pr,"No matrix given.")
      int k = _data->MatrixLabel[mat_name];
      CHK_MSGR(k==0,_pr,"Matrix "+mat_name+" not defined.")
      M = _data->theMatrix[k];
      int nr=M->getNbRows(), nc=M->getNbColumns();
      CHK_MSGR(nr!=nc,_pr,"Matrix "+mat_name+" must be a square matrix.")
      size = nr;
      CHK_MSGR(nb_eigv<0,_pr,"Illegal number of eigen values: "+to_string(nb_eigv))
      if (nb_eigv==0)
         nb_eigv = size;
      Alg = meth[method];
      if (Alg==OFELI::SUBSPACE)
         eig_vec = true;
      CHK_MSGR(Alg!=OFELI::SUBSPACE && Alg!=OFELI::QR,_pr,"Method "+to_string(Alg)+" not available")
      *_rita->ofh << "eigen matrix=" << mat_name;
      if (symm)
         *_rita->ofh << " symmetric";
      if (nb_eigv<size)
         *_rita->ofh << " nb=" << nb_eigv;
      *_rita->ofh << " method=" << method;
      if (evect=="")
         evect = mat_name+"ev_";
      if (eig_vec) {
         for (int i=0; i<nb_eigv; ++i) {
            evectR.push_back(evect+"R"+to_string(i+1));
            evectI.push_back(evect+"I"+to_string(i+1));
            _data->addVector(evectR[i],0.,size,"",SetCalc::SET);
            _data->addVector(evectI[i],0.,size,"",SetCalc::SET);
         }
         *_rita->ofh << " eigv";
      }
      evR = mat_name + "ev_R", evI = mat_name + "ev_I";
      _data->addVector(evR,0.,nb_eigv,"",SetCalc::SET);
      _data->addVector(evI,0.,nb_eigv,"",SetCalc::SET);
      *_rita->ofh << endl;
      _data->addEig(this,name);
      if (_verb)
         cout << "Eigen problem " << name << " created." << endl;
   }

   else {
   }
   return 0;
}


void eigen::print(ostream& s) const
{
   string nm="";
   s << "Eigen problem solver" << endl;
   s << "Matrix: " << " " << ", size: " << endl;

   for (int i=1; i<=nb_eigv; ++i) {
      if (symm) {
         nm = (*_data->theVector[_data->VectorLabel[evR]])(i);
         if (_data->dn[nm].active==false)
            s << "WARNING: Vector " << nm << " has been removed." << endl;
      }
      else {
         nm = (*_data->theVector[_data->VectorLabel[evR]])(i);
         if (_data->dn[nm].active==false)
            s << "WARNING: Vector " << nm << " has been removed." << endl;
         nm = (*_data->theVector[_data->VectorLabel[evI]])(i);
         if (_data->dn[nm].active==false)
            s << "WARNING: Vector " << nm << " has been removed." << endl;
      }
   }
   if (eig_vec) {
      for (int i=1; i<=nb_eigv; ++i) {
         int k=_data->VectorLabel[evectR[i-1]], l=_data->VectorLabel[evectI[i-1]];
         if (symm) {
            for (int j=1; j<=size; ++i) {
               nm = (*_data->theVector[k])(j);
               if (_data->dn[nm].active==false)
                  s << "WARNING: Vector " << nm << " has been removed." << endl;
            }
         }
         else {
            for (int j=1; j<=size; ++i) {
               nm = (*_data->theVector[k])(j);
               if (_data->dn[nm].active==false)
                  s << "WARNING: Vector " << nm << " has been removed." << endl;
               nm = (*_data->theVector[l])(j);
               if (_data->dn[nm].active==false)
                  s << "WARNING: Vector " << nm << " has been removed." << endl;
            }
         }
      }
   }
   if (solved==0) {
      s << "Eigenvalue problem unsolved !" << endl;
      return;
   }
   s << "Number of computed eigenvalues: " << nb_eigv << endl;
   for (int i=1; i<=nb_eigv; ++i) {
      s << "Eigenvalue #" << i << ": ";
      if (symm) {
         nm = (*_data->theVector[_data->VectorLabel[evR]])(i);
         if (_data->dn[nm].active==false)
            s << "WARNING: Vector " << nm << " has been removed." << endl;
         s << nm << endl;
      }
      else {
         s << (*_data->theVector[_data->VectorLabel[evR]])(i) << " + " 
           << (*_data->theVector[_data->VectorLabel[evI]])(i) << " I" << endl;
      }
   }
   if (!eig_vec || verbose==1)
      return;
   s << "Eigenvectors" << endl;
   for (int i=0; i<nb_eigv; ++i) {
      s << "Eigenvector " << i+1 << ": ";
      int k=_data->VectorLabel[evectR[i]], l=_data->VectorLabel[evectI[i]];
      if (symm)
         for (int j=1; j<=size; ++j)
            s << (*_data->theVector[k])(j) << endl;
      else
         for (int j=1; j<=size; ++j)
            s << (*_data->theVector[k])(j) << " + " << (*_data->theVector[l])(j) << " I" << endl;
   }
   if (solved==0) {
      cout << "Eigenproblem unsolved !" << endl;
      return;
   }
}


ostream& operator<<(ostream& s, const eigen& es)
{
   es.print(s);
   return s;
}

} /* namespace RITA */
