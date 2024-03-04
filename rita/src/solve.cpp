/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2024 Rachid Touzani

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

                         Implementation of class 'solve'

  ==============================================================================*/

#include "rita.h"
#include "solve.h"
#include "configure.h"
#include "data.h"
#include "calc.h"
#include "ls.h"
#include "ae.h"
#include "ode.h"
#include "pde.h"
#include "optim.h"
#include "eigen.h"
#include "transient.h"
#include "defs.h"
#include "util/macros.h"
#include "helps.h"

namespace RITA {

solve::solve()
{
}


solve::solve(rita *r, cmd *command, configure *config)
      : _rita(r), _data(r->_data), _configure(config), _cmd(command), _verb(1),
        _nb_pb(0), _nb_ls(_data->nbAllLS()), _nb_ae(_data->nbAllAE()),
        _nb_pde(_data->nbAllPDE()), _nb_opt(_data->nbAllOpt()),
        _nb_eigen(_data->nbAllEig()), _analysis(analysis_type::NONE)
{
}


int solve::run()
{
   int ret=0;
   string fn="", a="", p="";
   vector<string> var;
   static const vector<string> kw {"prob$lem","an$alysis"};

   int nb_pb = _data->nb_ls+_data->nb_ae+_data->nb_ode+_data->nb_pde+_data->nb_opt+_data->nb_eig;
   CHK_MSGR(nb_pb==0,_pr,"No problem to solve.")

// Case of a unique defined problem
   int nb=0;
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   for (size_t k=0; k<nb_args; ++k) {

      size_t n = _cmd->getArgs(nb);
      switch (n) {

         case   0:
            for (int i=0; i<nb; ++i) {
               p = _cmd->string_token(i);
               _pb.push_back(p);
               _nb_pb++;
            }
            break;

         case   1:
            a = _cmd->string_token(0);
            break;

         case 100:
            cout << "Available arguments: " << Solve_help << endl;
            return 0;

         case 101:
            cout << Solve_Help << endl;
            return 0;

         default:
            MSGR(_pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   if (nb_args>0) {
      CHK_MSGR(_nb_pb>nb_pb,_pr,"Number of problems to solve is larger than number of defined problems.")
      for (auto const& p: _pb)
         CHK_MSGR(CHK_PB(p),_pr,p+" does not define a problem to solve")
      if (_nb_pb==0) {
         _pb.push_back(_data->p2s);
         _nb_pb = 1;
      }
      *_rita->ofh << "solve " << _pb[0];
      for (auto const& p: _pb)
         *_rita->ofh << "," << p;
      _analysis = _rita->_ant[a];
      if (a!="")
         *_rita->ofh << " analysis=" << a;
      *_rita->ofh << endl;
      run_pb();
   }
   else {
      *_rita->ofh << "solve" << endl;
      if (nb_pb==1) {
         _pb.push_back(_data->p2s);
         _nb_pb = 1;
         return run_pb();
      }
      while (1) {
         nb = _cmd->readline(sPrompt+"solve> ");
         if (nb<0)
            continue;
         int key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw);
         if (key>=200) {
            _data->setDataExt(key);
            continue;
         }
         switch (key) {

            case   0:
               ret = _cmd->get(p);
               if (!ret) {
                  CHK_MSGB(_nb_pb==nb_pb,_pr,"Number of problems to solve is larger than number of defined problems.")
                  CHK_MSGB(CHK_PB(p),_pr,p+" does not define a problem to solve")
                  CHK_MSGB(_data->dn.count(p)==0,_pr,"Problem "+p+" undefined.")
                  CHK_MSGB(_data->dn[p].active==false,_pr,"Problem "+p+" has been removed.")
                  _pb.push_back(p);
                  _nb_pb++;
               }
               break;

            case   1:
               ret = _cmd->get(a);
               if (!ret) {
                  CHK_MSGB(_rita->_ant.count(a)==0,_pr,"Unknown analysis type "+a)
                   _analysis = _rita->_analysis_type = _rita->_ant[a];
                  *_rita->ofh << "  analysis " << a << endl;
               }
               break;

            case 100:
               cout << "Available Commands: " << Solve_help << endl;
               break;

            case 101:
               cout << "\nAvailable Commands:\n" << Solve_HHelp << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               _configure->run();
               break;

            case 104:
            case 105:
               *_rita->ofh << "  end" << endl;
               return run_pb();

            case 106:
               _rita->setEcho(nb_args);
               break;

            case -2:
               break;

            default:
               DEFAULT_KW
         }
      }
   }
   return 0;
}


int solve::run_pb()
{
   if (_analysis==analysis_type::NONE) {
      for (auto const& p: _pb) {
         if (_data->dn.count(p) && _data->dn[p].active) {
            if (_data->dn[p].dt==DataType::LS)
               _analysis = _rita->_analysis_type = analysis_type::STEADY_STATE;
            else if (_data->dn[p].dt==DataType::AE)
               _analysis = _rita->_analysis_type = analysis_type::STEADY_STATE;
            else if (_data->dn[p].dt==DataType::ODE)
               _analysis = _rita->_analysis_type = analysis_type::TRANSIENT;
            else if (_data->dn[p].dt==DataType::PDE)
               _analysis = _data->thePDE[_data->PDELabel[p]]->analysis;
            else if (_data->dn[p].dt==DataType::OPTIM)
               _analysis = analysis_type::OPTIMIZATION;
            else if (_data->dn[p].dt==DataType::EIGEN)
               _analysis = analysis_type::EIGEN;
         }
      }
   }

   if (_verb) {
      cout << "\nSolving problem(s): " << _pb[0];
      for (int i=1; i<_nb_pb; ++i)
         cout << "," << _pb[i];
      cout << " ..." << endl;
   }

   int ret = 0;
   switch (_analysis) {

      case analysis_type::NONE:
         break;

      case analysis_type::STEADY_STATE:
         ret = run_steady();
         break;

      case analysis_type::TRANSIENT:
         {
            transient ts(_rita,_pb);
            ret = ts.run(_pb);
            break;
         }

      case analysis_type::OPTIMIZATION:
         {
            int e = _data->OptLabel[_pb[0]];
            if (e>0)
               _optim = _data->theOpt[e];
            ret = run_optim();
            break;
         }

      case analysis_type::EIGEN:
         {
            if (_verb)
               cout << "Running eigen problem solver ..." << endl;
            int e = _data->EigLabel[_pb[0]];
            if (e>0)
               _eigen = _data->theEig[e];
            ret = run_eigen();
            break;
         }

      default:
         break;

   }
   return ret;
}


DataType solve::checkPb(string p)
{
   if (_data->LSLabel.count(p))
      return DataType::LS;
   else if (_data->AELabel.count(p))
      return DataType::AE;
   else if (_data->ODELabel.count(p))
      return DataType::ODE;
   else if (_data->PDELabel.count(p))
      return DataType::PDE;
   else if (_data->OptLabel.count(p)) {
      _rita->_analysis_type = analysis_type::OPTIMIZATION;
      return DataType::OPTIM;
   }
   else if (_data->EigLabel.count(p)) {
      _rita->_analysis_type = analysis_type::EIGEN;
      return DataType::EIGEN;
   }
   return DataType::NOTHING;
}


int solve::run_steady()
{
   double err2=0., errI=0.;
   int ret = 0;
   try {
      for (auto const& p: _pb) {
         if (_data->dn.count(p)) {

//          Linear system
//          -------------

            if (_data->dn[p].dt==DataType::LS) {
               _ls = _data->theLS[_data->LSLabel[p]];
               LinearSolver lsolv;
               lsolv.setMatrix(_data->theMatrix[_ls->iMat]);
               lsolv.setRHS(*_data->theVector[_ls->iRHS]);
               lsolv.setSolution(*_data->theVector[_ls->iSol]);
               ret = lsolv.solve(_ls->lls,_ls->pprec);
               _ls->solved = 1 - ret;
            }

//          Algebraic equation
//          ------------------

            if (_data->dn[p].dt==DataType::AE) {
               _ae = _data->theAE[_data->AELabel[p]];
               int f = _ae->ivect[0];
               _data->setVectorValue(f,_data->theVector[f]);
               NLASSolver nls(_ae->nnls,_ae->size);
               if (_ae->size==1)
                  nls.setInitial(_data->theVector[f]->at(0));
               else
                  nls.setInitial(*_data->theVector[f]);
               for (int i=0; i<_ae->size; ++i)
                  nls.setf(*_data->theFct[_ae->iFct[i]]->getFct());
               ret = nls.run();
               _ae->solved = 1 - ret;
               _data->setVectorValue(f,_data->theVector[f]);
               if (_ae->size==1)
                  cout << "Solution to algebraic equation: " << _data->theVector[f]->at(0) << endl;
            }

//          Partial Differential Equation
//          -----------------------------

            if (_data->dn[p].dt==DataType::PDE) {
               _pde = _data->thePDE[_data->PDELabel[p]];

//             Set solution vector and linear solver
               for (int i=0; i<_pde->nb_vectors; ++i)
                  _pde->theEquation->setInput(EType::SOLUTION,*_data->theVector[_pde->fd[i].vect]);
               _pde->theEquation->setSolver(_pde->lls,_pde->pprec);

//             Set boundary condition
               if (_pde->set_bc) {
                  if (_pde->bc.withRegex(1)) {
                     for (auto const& v: _pde->bc_data.cexp)
                        _pde->bc.setNodeBC(v.first,v.second);
                  }
                  _pde->theEquation->setInput(EType::BOUNDARY_CONDITION,_pde->bc);
               }

//             Set body force
               if (_pde->set_bf) {
                  if (_pde->bf.withRegex(1))
                     _pde->bf.set(_pde->bf_data.exp);
                  _pde->theEquation->setInput(EType::BODY_FORCE,_pde->bf);
               }

//             Set surface force (Neumann)
               if (_pde->set_sf) {
                  if (_pde->sf.withRegex(1)) {
                     for (auto const& v: _pde->sf_data.cexp)
                        _pde->sf.setSideBC(v.first,v.second);
                  }
                  _pde->theEquation->setInput(EType::BOUNDARY_FORCE,_pde->sf);
               }

//             Solve problem
               _pde->solved = 1 - _pde->theEquation->run();

//             Compute error if analytical solution is given
               if (_pde->get_err>0) {
                  CHK_MSGB(!_pde->solved,_pr,"Problem has not been solved.")
                  getPDEError(err2,errI);
                  cout << "L2 and Max Errors: " << err2 << "  " << errI << endl;
                  if (_pde->e2=="") {
                     _pde->e2 = "e2_" + to_string(_data->iParam+1);
                     _pde->eI = "eI_" + to_string(_data->iParam+1);
                  }
                  _data->addParam(_pde->e2,err2,false);
                  _data->addParam(_pde->eI,errI,false);
               }
            }
         }
      }
   } CATCH
   return ret;
}


int solve::run_eigen()
{
   CHK_MSGR(!_eigen->log,_pr+"run_eigen>","Eigenproblem undefined or improperly defined.")
   try {
      EigenProblemSolver es;
      es.setMatrix(_eigen->M);
      if (_eigen->symm) {
         es.set(OFELI::SUBSPACE,true);
         es.setNbEigv(_eigen->nb_eigv);
      }
      else {
         es.set(OFELI::QR,false);
         if (_eigen->eig_vec)
            es.setEigenVectors();
      }
      es.run();
      Vect<double> vvR, vvI;
      for (int i=1; i<=_eigen->nb_eigv; ++i) {
         vvR.push_back(es.getEigenValue(i,1));
         vvI.push_back(es.getEigenValue(i,2));
      }
      _data->setVectorValue(_data->VectorLabel[_eigen->evR],&vvR);
      _data->setVectorValue(_data->VectorLabel[_eigen->evI],&vvI);
      for (int i=0; i<_eigen->nb_eigv; ++i) {
         if (_eigen->eig_vec) {
            int kR=_data->VectorLabel[_eigen->evectR[i]], kI=_data->VectorLabel[_eigen->evectI[i]];
            es.getEigenVector(i+1,*(_data->theVector[kR]),*(_data->theVector[kI]));
            _rita->_calc->setVectorValue(_eigen->evectR[i],_data->theVector[kR]);
            _rita->_calc->setVectorValue(_eigen->evectI[i],_data->theVector[kI]);
         }
      }
      cout << "Eigenvalues stored in vectors: " << _eigen->evR << ", " << _eigen->evI << endl;
      if (_eigen->eig_vec)
         cout << "Eigenvectors stored in vectors: " << _eigen->evectR[0] << " .. " << _eigen->evectR[_eigen->nb_eigv-1]
              << ", " << _eigen->evectI[0] << " .. " << _eigen->evectI[_eigen->nb_eigv-1] << endl;
      _eigen->solved = 1;
   } CATCH
   return 0;
}


int solve::run_optim()
{
   CHK_MSGR(_optim->log,_pr+"run_optim>","Optimization problem undefined or improperly defined.")
   try {
      int ret=0, size=_optim->size;
      if (_optim->lp) {
         LPSolver s;
         s.setSize(size,_optim->nb_lec,_optim->nb_gec,_optim->nb_eqc);
         s.set(*_optim->opt_var);
         s.set(OFELI::LPSolver::OBJECTIVE,_optim->a,_optim->b);
         for (int i=0; i<_optim->nb_eqc; ++i)
            s.set(OFELI::LPSolver::EQ_CONSTRAINT,*_optim->a_eq[i],_optim->b_eq[i]);
         for (int i=0; i<_optim->nb_lec; ++i)
            s.set(OFELI::LPSolver::LE_CONSTRAINT,*_optim->a_le[i],_optim->b_le[i]);
         for (int i=0; i<_optim->nb_gec; ++i)
            s.set(OFELI::LPSolver::GE_CONSTRAINT,*_optim->a_ge[i],_optim->b_ge[i]);
         ret = s.run();
         if (ret==0)
            _optim->obj = s.getObjective();
      }
      else {
         OptSolver s(*_optim->opt_var);
         s.setOptMethod(_optim->Alg);
         s.setObjective(*_optim->J_Fct->getFct());
         if (_optim->G_ok) {
            for (int i=0; i<size; ++i)
               s.setGradient(*_data->theFct[_optim->igrad+i]->getFct(),i+1);
         }
         if (_optim->H_ok) {
            for (int i=0; i<size*size; ++i)
               s.setHessian(*_data->theFct[_optim->ihess+i]->getFct(),i+1);
         }
         for (int i=0; i<_optim->nb_lec; ++i)
            s.setIneqConstraint(*_optim->inC_Fct[i]->getFct(),_optim->penal);
         for (int i=0; i<_optim->nb_eqc; ++i)
            s.setEqConstraint(*_optim->eqC_Fct[i]->getFct(),_optim->penal);
         s.setLowerBounds(_optim->lb);
         s.setUpperBounds(_optim->ub);
         ret = s.run();
         if (ret==0)
            _optim->obj = s.getObjective();
      }
      if (ret==0) {
         _data->setVectorValue(_data->VectorLabel[_optim->var_name[0]],_optim->opt_var);
         cout << "Optimization variable stored in vector: " << _optim->var_name[0] << endl;
      }
      _optim->solved = 1 - ret;
   } CATCH
   return 0;
}


void solve::getPDEError(double &e2, double &eI)
{
   e2 = eI = 0.;
   Vect<double> &u = *_data->theVector[_pde->fd[0].vect];
   int nb=0, nb_dof=_pde->fd[0].nb_dof, k=0;
   if (_pde->theMesh!=nullptr) {
      Mesh &ms = *_pde->theMesh;
      if (ms.NodesAreDOF()) {
         nb = ms.getNbNodes(), k = 0;
         for (int i=1; i<=nb; ++i) {
            Point<double> x=ms[i]->getCoord();
            for (int j=0; j<nb_dof; ++j) {
               funct &f = *_pde->SolFct[j];
               double ui=u[k++], vi=f(x,0.);
               e2 += (ui-vi)*(ui-vi);
               eI = fmax(eI,fabs(ui-vi));
            }
         }
      }
      else if (ms.ElementsAreDOF())
         nb = ms.getNbElements(), k = 0;
      else if (ms.SidesAreDOF())
         nb = ms.getNbSides(), k = 0;
      e2 = sqrt(e2/(nb_dof*nb));
   }
   else if (_pde->theGrid!=nullptr) {
      Grid &gr = *_pde->theGrid;
      for (size_t i=1; i<=gr.getNx(); ++i) {
         for (size_t j=1; j<=gr.getNy(); ++j) {
            for (size_t k=1; k<=gr.getNz(); ++k) {
               double ui=u(i,j,k), vi=(*_pde->SolFct[0])(gr.getXYZ(i,j,k),_pde->final_time);
               e2 += (ui-vi)*(ui-vi);
               eI = fmax(eI,fabs(ui-vi));
            }
         }
      }
      e2 = sqrt(e2/(gr.getNx()*gr.getNy()*gr.getNz()));
   }
}

} /* namespace RITA */