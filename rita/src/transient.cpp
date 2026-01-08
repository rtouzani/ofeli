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

                      Implementation of class 'transient'

  ==============================================================================*/

#include "OFELI_Config.h"
#include "util/util.h"
#include "io/IOField.h"
#include "io/saveField.h"
#include <iostream>

#include "transient.h"
#include "calc.h"
#include "ls.h"
#include "ae.h"
#include "ode.h"
#include "pde.h"
#include "defs.h"
#include "solvers/ODESolver.h"
#include "solvers/NLASSolver.h"

namespace RITA {

transient::transient(rita *r, vector<string> &pb)
          : _rita(r), _data(r->_data), _nb_ls(_data->nbAllLS()), _nb_ae(_data->nbAllAE()),
            _nb_ode(_data->nbAllODE()), _nb_pde(_data->nbAllPDE())
{
   _ode_eq.resize(_nb_ode+1);
   _nlas_eq.resize(_nb_ae+1);
   _ts_eq.resize(_nb_pde+1);
   _pb = pb;
   ret = 0;
}


int transient::run(const vector<string>& pb)
{
   Verbosity = 1;
   try {
      theStep = 1;
      for (auto const& p: pb) {
         if (_data->dn.count(p)) {

//          Algebraic Equation
            if (_data->dn[p].dt==DataType::AE) {
               int e = _data->AELabel[p];
               _ae = _data->theAE[e];
               _nlas_eq[e].set(_ae->nnls);
               _nlas_eq[e].setNbEq(_ae->size);
               for (int i=0; i<_ae->size; ++i)
                  _nlas_eq[e].setf(*_data->theFct[_ae->iFct[i]]->getFct());
               if (_ae->size==1)
                  cout << "Solution to algebraic equation: " << _data->theVector[_ae->ivect[0]]->at(0) << endl;
            }

//          Ordinary Differential Equation
            if (_data->dn[p].dt==DataType::ODE) {
               int e = _data->ODELabel[p];
               _ode = _data->theODE[e];
               theTimeStep = _ode->time_step;
               theFinalTime = _ode->final_time;
               _ode_eq[e].set(_ode->Scheme);
               _ode_eq[e].setNbEq(_ode->size);
               int f = _ode->ivect[0];
               if (_ode->size==1)
                  _ode_eq[e].setInitial(_data->theVector[f]->at(0));
               else
                  _ode_eq[e].setInitial(*_data->theVector[f]);
               string fh = _data->vect_hist[_ode->var_name[0]];
               if (fh!="%$§&")
                  _data->theHVector[_data->HVectorLabel[fh]]->set(*_data->theVector[f],theTime);
               for (int i=0; i<_ode->size; ++i)
                  _ode_eq[e].setF(*_data->theFct[_ode->iFct[i]]->getFct());
            }

//          Partial Differential Equation
            if (_data->dn[p].dt==DataType::PDE) {
               int e = _data->PDELabel[p];
               _pde = _data->thePDE[e];
               theTimeStep = _pde->time_step;
               theFinalTime = _pde->final_time;
               theTime = 0.;
               _ts_eq[e].set(_pde->Scheme,theTimeStep,theFinalTime);
               _ts_eq[e].setPDE(*_pde->theEquation);
               _ts_eq[e].setLinearSolver(_pde->lls,_pde->pprec);
               for (int i=0; i<_pde->nb_vectors; ++i) {
                  int f = _pde->fd[i].vect;
                  if (i==0) {
                     Vect<double> *u = _data->theVector[_pde->fd[0].vect];
                     if (u->withRegex(1))
                        u->set(_pde->in_data.exp);
                     _ts_eq[e].setInitial(*u);
                     _pde->theEquation->setInput(EType::INITIAL,*u);
                  }
                  if (_pde->eq=="incompressible-navier-stokes" && _pde->Sdm==pde::FE_P1)
                     _pde->theEquation->setInput(EType::PRESSURE,*_data->theVector[_pde->fd[1].vect]);
                  string fh = _data->vect_hist[_pde->fd[i].fn];
                  _data->theVector[f]->setTime(theTime);
                  _data->theVector[f]->setName(_data->VectorName[f]);
                  if (fh!="")
                     _data->theHVector[_data->HVectorLabel[fh]]->set(*_data->theVector[f],theTime);
               }
            }
         }
      }

//    LOOP ON TIME STEPS
      TimeLoop {

         if (_rita->_verb)
            cout << "Performing time step " << theStep << ", Time = " << theTime << endl;

         for (auto const& p: pb) {
            if (_data->dn.count(p)) {

//             Algebraic Equation
               if (_data->dn[p].dt==DataType::AE) {
                  int e = _data->AELabel[p];
                  _ae = _data->theAE[e];
                  MSGR(_pr,"No algebraic equation solver implemented.")
               }

//             Ordinary Differential Equation
               if (_data->dn[p].dt==DataType::ODE) {
                  int e = _data->ODELabel[p];
                  _ode = _data->theODE[e];
                  int f = _ode->ivect[0];
                  _ode_eq[e].runOneTimeStep();
                  if (_ode->size==1)
                     _data->theVector[f]->at(0) = _ode_eq[e].get();
                  _data->setVectorValue(f,_data->theVector[f]);
                  string fh = _data->vect_hist[_ode->var_name[0]];
                  if (fh!="%$§&")
                     _data->theHVector[_data->HVectorLabel[fh]]->set(*_data->theVector[f],theTime);
                  if (_ode->phase.size()) {
                     Vect<double> dydt(_ode->size), dy(1);
                     _ode_eq[e].getTimeDerivative(dydt);
                     for (size_t i=0; i<_ode->phase.size(); ++i) {
                        int j = _ode->iphase[i];
                        dy[0] = dydt[j];
                        _data->add2History(_ode->phase[i],dy,_data->theVector[_ode->ivect[0]]->at(j));
                     }
                  }
               }

//             Partial Differential Equation
               if (_data->dn[p].dt==DataType::PDE) {
                  int e = _data->PDELabel[p];
                  _pde = _data->thePDE[e];
                  if (_pde->set_bf) {
                     _pde->bf.setTime(theTime);
                     if (_pde->bf.withRegex(1))
                        _pde->bf.set(_pde->bf_data.exp);
                     _ts_eq[e].setRHS(_pde->bf);
                     _pde->theEquation->setInput(EType::BODY_FORCE,_pde->bf);
                  }

                  if (_pde->set_bc) {
                     _pde->bc.setTime(theTime);
                     if (_pde->bc.withRegex(1)) {
                        for (auto const& v: _pde->bc_data.cexp)
                           _pde->setNodeBC(v.first,v.second,theTime,_pde->bc);
                     }
                     _ts_eq[e].setBC(_pde->bc);
                  }

                  if (_pde->set_sf) {
                     _pde->sf.setTime(theTime);
                     if (_pde->sf.withRegex(1)) {
                        for (auto const& v: _pde->sf_data.cexp)
                           _pde->setNodeBC(v.first,v.second,theTime,_pde->sf);
                     }
//                  ts.setSF(_data->sf[i]);
                  }

                  _ts_eq[e].runOneTimeStep();
                  for (int i=0; i<_pde->nb_vectors; ++i) {
                     int f = _pde->fd[i].vect;
                     _data->setVectorValue(f,_data->theVector[f]);
                     string fh = _data->vect_hist[_pde->fd[i].fn];
                     _data->theVector[f]->setTime(theTime);
                     _data->theVector[f]->setName(_data->VectorName[f]);
                     if (fh!="")
                        _data->theHVector[_data->HVectorLabel[fh]]->set(*_data->theVector[f],theTime);
                  }
               }
            }
         }
      }

//    Problems are solved
      for (auto const& p: pb) {
         if (_data->dn.count(p)) {
            if (_data->dn[p].dt==DataType::AE)
               _data->theAE[_data->AELabel[p]]->solved = 1;
            if (_data->dn[p].dt==DataType::ODE)
               _data->theODE[_data->ODELabel[p]]->solved = 1;
            if (_data->dn[p].dt==DataType::PDE)
               _data->thePDE[_data->PDELabel[p]]->solved = 1;
         }
      }

//    Compute error if analytical solution is given
      for (auto const& p: pb) {
         if (_data->dn.count(p)) {
            if (_data->dn[p].dt==DataType::ODE) {
               _ode = _data->theODE[_data->ODELabel[p]];
               if (_ode->get_err>0) {
                  CHK_MSGB(!_ode->solved,_pr,"Problem has not been solved.")
                  double err=0.;
                  for (int i=0; i<_ode->size; ++i) {
                     double u=(*_data->theVector[_ode->ivect[0]])[i], v=(*_ode->SolFct[i])(_ode->final_time);
                     err = fmax(err,fabs(u-v));
                  }
                  cout << "Error: " << err << endl;
                  if (_ode->err!="")
                     _data->addParam(_ode->err,err,SetCalc::SET);
               }
            }
            if (_data->dn[p].dt==DataType::PDE) {
               _pde = _data->thePDE[_data->PDELabel[p]];
               if (_pde->get_err>0) {
                  CHK_MSGR(!_pde->solved,_pr,"Problem has not been solved.")
                  double err2=0., errI=0.;
                  getPDEError(err2,errI);
                  cout << "L2 and Max Errors: " << err2 << "  " << errI << endl;
                  if (_pde->e2!="") {
                     _data->addParam(_pde->e2,err2,SetCalc::SET);
                     _data->addParam(_pde->eI,errI,SetCalc::SET);
                  }
               }
            }
         }
      }

   } CATCH

   return 0;
}


void transient::getPDEError(double &e2, double &eI)
{
   e2 = eI = 0.;
   OFELI::Vect<double> *u = _data->theVector[_pde->fd[0].vect];
   OFELI::Mesh *ms = _pde->theMesh;
   OFELI::Grid *gr = _pde->theGrid;
   int nb_dof=_pde->fd[0].nb_dof;
   if (ms!=nullptr) {
      if (ms->NodesAreDOF()) {
         int nb=ms->getNbNodes(), k=0;
         for (int i=1; i<=nb; ++i) {
            OFELI::Point<double> x = (*ms)[i]->getCoord();
            for (int j=0; j<nb_dof; ++j) {
               double ui=(*u)[k++], vi=(*_pde->SolFct[j])(x.x,x.y,x.z,_pde->final_time);
               e2 += (ui-vi)*(ui-vi);
               eI = fmax(eI,fabs(ui-vi));
            }
         }
         e2 = sqrt(e2/k);
      }
      else if (ms->ElementsAreDOF()) {
         cout << "The case of element based solutions is not implemented yet." << endl;
         *_rita->ofl << "The case of element based solutions is not implemented yet." << endl;
         return;
      }
      else if (ms->SidesAreDOF()) {
         cout << "The case of side based solutions is not implemented yet." << endl;
         *_rita->ofl << "The case of side based solutions is not implemented yet." << endl;
         return;
      }
   }
   else if (gr!=nullptr) {
      for (size_t i=1; i<=gr->getNx(); ++i) {
         for (size_t j=1; j<=gr->getNy(); ++j) {
            for (size_t k=1; k<=gr->getNz(); ++k) {
               OFELI::Point<double> p = gr->getXYZ(i,j,k);
               double ui=(*u)(i,j,k), vi=(*_pde->SolFct[0])(p.x,p.y,p.z,_pde->final_time);
               e2 += (ui-vi)*(ui-vi);
               eI = fmax(eI,fabs(ui-vi));
            }
         }
      }
      e2 = sqrt(e2/(gr->getNx()*gr->getNy()*gr->getNz()));
   }
}

} /* namespace RITA */
