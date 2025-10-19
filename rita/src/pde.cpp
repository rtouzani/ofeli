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

                       Implementation of class 'pde'

  ==============================================================================*/

#include "configure.h"
#include "pde.h"
#include "cmd.h"
#include "rita.h"
#include "calc.h"
#include "defs.h"
#include "helps.h"

namespace RITA {

pde::pde(rita*      r,
         cmd*       command,
         configure* config)
    : _rita(r), _configure(config), _cmd(command), _data(r->_data)
{
   theMesh = nullptr;
   theGrid = nullptr;
   eq = "laplace";
   nb_vectors = 0;
   ls = "direct";
   nls = "cg";
   prec = "dilu";
   lls = _data->Ls[ls];
   pprec = _data->Pr[prec];
   axi = false;
   _verb = _rita->_verb;
   ieq = pde_map["laplace"];
   _c00_set = _c10_set = _c01_set = _c20_set = _c02_set = false;
   _rho_set = _Cp_set = _kappa_set = _mu_set = _sigma_set = _Mu_set = false;
   _epsilon_set = _omega_set = _beta_set = _v_set = _young_set = _poisson_set = false;
   set_u = set_bc = set_bf = set_sf = set_in = set_coef = false;
   every = 1;
   e2 = eI = "";
   file = "";
   solved = 0;
}


pde::~pde()
{
   if (theEquation!=nullptr)
      delete theEquation;
}


void pde::set(string e)
{
   ieq = pde_map[e];
   eq = e;
}


int pde::run()
{
   string pde_id, name="", fn="", lin_solv="", an="";
   bool space_d=false;
   e2 = eI = "";
   int ret=0, every=0, size=0;
   CHK_MSG1R(_cmd->setNbArg(1,"Give PDE name."),_pr,"Missing pde name.","",1)
   static const vector<string> pn {"linear","laplace","heat","wave","transport","linear-elasticity",
                                   "truss","beam","incompressible-navier-stokes"};
   get_err = 0;
   ret = _cmd->get(pn,pde_id);
   CHK_MSG1R(ret<0,_pr,"Unknown pde "+pde_id,"Unknown pde name. Available pde's:\n"
                       "linear, laplace, heat, wave, transport, linear-elasticity, truss, beam,"
                       " incompressible-navier-stokes",1)
   set(pde_id);
   analysis = analysis_type::STEADY_STATE;
   *_rita->ofh << "pde " << pde_id << endl;

   size_t nb=0;
   int nbv=0, nb_args=0, key=0;
   bool vector_ok = false;
   string str = "", str1 = "", ff="";
   scheme = "backward-euler";
   set(_cmd);
   log.vect = true;
   CHK_MSGR(_data->nb_meshes==0&&_data->nb_grids==0,_pr,"No available mesh or grid")
   CHK_MSGR(_data->theMesh[_data->iMesh]->getNbNodes()==0,_pr,"Empty mesh")
   log.mesh = false;
   ls = OFELI::CG_SOLVER;
   lsolv = "cg", lprec = "dilu";
   prec = OFELI::DILU_PREC;
   string spd = "feP1";
   static const vector<string> kw {"var$iable","vect$or","coef","axi","in$it","bc","bf","source","sf",
                                   "traction","space","ls","nls","final$-time","time-step","scheme",
                                   "name","analytic","analytic-function","err","save-every"};

   while (1) {
      READLINE
      switch (key) {

         case   0:
         case   1:
            CHK_MSG1B(_cmd->setNbArg(1,"Give name of an associated variable/vector."),_pr+"variable>",
                                       "Missing name of an associated vector.","",1)
            if (!_cmd->get(str)) {
               CHK_MSGB(_data->VectorLabel.count(str)>0 && _data->dn[str].active,_pr+"variable>",
                        "Variable name "+str+" already used.")
               _vect.push_back(str);
               nbv++;
               vector_ok = true;
               *_rita->ofh << "  vector " << str << endl;
            }
            break;

         case   2:
            set_coef = true;
            setCoef();
            break;

         case   3:
            axi = true;
            break;

         case   4:
            ret = getIn();
            break;

         case   5:
            ret = getBC();
            break;

         case   6:
         case   7:
            ret = getBF();
            break;

         case   8:
         case   9:
            ret = getSF();
            break;

         case  10:
            if (_cmd->setNbArg(1,"Give space discretization method.")) {
               _rita->msg(_pr+"space>","Missing space discretization method.","",1);
               cout << "Available Arguments\n";
               cout << "fd:   Finite Differences\n";
               cout << "feP1: P1 finite elements\n";
               cout << "feP2: P2 finite elements\n";
               cout << "feQ1: Q1 finite elements\n";
               cout << "fv:   Finite volumes\n";
               cout << "dg:   Discontinuous Galerkin" << endl;
               break;
            }
            _cmd->get(spd);
            ret = 1;
            CHK_MSGB(pde_sdm.find(spd)==pde_sdm.end(),_pr+"space>",
                     "Unknown or unimplemented space discretization method: "+spd)
            ret = 0;
            *_rita->ofh << "  space " << spd << endl;
            space_d = true;
            break;

         case  11:
            if (_cmd->setNbArg(1,"Linear solver and optional preconditioner to be supplied.",1)) {
               _rita->msg(_pr+"ls>","Missing linear solver data.","",1);
               ret = 1;
               break;
            }
            nb = _cmd->getNbArgs();
            if (nb==0)
               _rita->msg(_pr+"ls>","Missing linear solver data.");
            ret = _cmd->get(str);
            str1 = "ident";
            if (nb>1)
               ret += _cmd->get(str1);
            if (!ret) {
               *_rita->ofh << "  ls " << str << " " << str1 << endl;
               if (!set_ls(str,str1)) {
                  lsolv = str;
                  lprec = str1;
                  ls = str;
                  prec = str1;
                  lls = _data->Ls[ls];
                  pprec = _data->Pr[prec];
               }
            }
            else
               log.ls = true;
            break;

         case  12:
            if (_cmd->setNbArg(1,"Nonlinear solver to be supplied.",1)) {
               _rita->msg(_pr+"nls>","Missing nonlinear solver data.","",1);
               ret = 1;
               break;
            }
            ret = _cmd->get(str);
            if (!ret) {
               *_rita->ofh << "  nls " << str << endl;
               ret = set_nls(str);
               if (!ret)
                  nls = str;
               else
                  log.nl = true;
            }
            break;

         case  13:
            ret = _cmd->get(final_time);
            if (!ret) {
               *_rita->ofh << "  final-time " << final_time << endl;
               analysis = analysis_type::TRANSIENT;
            }
            break;

         case  14:
            ret = _cmd->get(time_step);
            if (!ret) {
               *_rita->ofh << "  time-step " << time_step << endl;
               analysis = analysis_type::TRANSIENT;
            }
            break;

         case  15:
            ret = _cmd->get(scheme);
            if (!ret) {
               *_rita->ofh << "  scheme " << scheme << endl;
               analysis = analysis_type::TRANSIENT;
            }
            break;

         case  16:
            ret = _cmd->get(name);
            if (!ret)
               *_rita->ofh << "  name " << name << endl;
            break;
 
         case 17:
            ret = _cmd->get(an);
            if (!ret) {
               *_rita->ofh << "  analytic " << an << endl;
               analytic.push_back(an);
               get_err = 1;
            }
            break;

         case 18:
            ret = _cmd->get(an);
            if (!ret) {
               *_rita->ofh << "  analytic-function " << an << endl;
               analytic_f.push_back(an);
               get_err = 2;
            }
            break;

         case  19:
            ret  = _cmd->get(e2);
            ret += _cmd->get(eI);
            if (!ret)
               *_rita->ofh << "  err " << e2 << " " << eI << endl;
            break;
            
         case  20:
            ret = _cmd->get(every);
            if (!ret)
               *_rita->ofh << "  save-every " << every << endl;
            break;

         case 100:
            cout << PDE_help << endl;
            break;

         case 101:
            cout << PDE_Help << endl;
            break;

         case 102:
            _rita->getLicense();
            break;

         case 103:
            ret = _configure->run();
            break;

         case 104:
         case 105:
            _cmd->setNbArg(0);
            if (!space_d) {
               _rita->msg(_pr+"end>","No space discretization method given. No PDE data created.");
               *_rita->ofh << "end" << endl;
               return ret;
            }
            if (ret<0) {
               _rita->msg(_pr+"end>","No PDE data created.");
               *_rita->ofh << "end" << endl;
               return ret;
            }
            if (!vector_ok) {
               _rita->msg(_pr+"end>","No vector(s) defined for PDE.");
               break;
            }
            setSpD(spd);
            setEq();
            *_rita->ofh << "  space " << spd << endl;
            CHK_MSGB(nbv>nb_vectors,_pr+"end>","Too many vectors for defined PDE.")
            CHK_MSGB(nbv<nb_vectors,_pr+"end>","Not enough vectors for defined PDE.")
            Scheme = _Sch[scheme];
            ff = spd.substr(0,2);
            _data->addPDE(this,name);
            if (_verb)
               cout << "Partial Differential Equation " << name << " created." << endl;
            for (int i=0; i<nb_vectors; ++i) {
               if (ff=="fd")
                  fd[i].vect = _data->addGridVector(_data->GridName[_data->iGrid],_vect[i],fd[i].nb_dof,SetCalc::SET);
               else if (ff=="fe" || ff=="fv" || ff=="dg")
                  fd[i].vect = _data->addMeshVector(_data->MeshName[_data->iMesh],_vect[i],data::DataSize::NODES,fd[i].nb_dof,SetCalc::SET);
            }
            log.vect = false;
            b.setSize(_data->theMesh[_data->iMesh]->getNbEq());
            if (set_in)
               setIn();
            if (set_bc)
               setBC();
            if (set_sf)
               setSF();
            if (set_bf)
               setBF();
            size = std::max(analytic.size(),analytic_f.size());
            if (get_err==1) {
               for (int i=0; i<size; ++i) {
                  string nm="";
                  int k = _data->addFunction(nm,_var,analytic[i]);
                  FCT_NOT_DEFINED("end>",nm)
                  FCT_ALREADY_DEFINED("end>",nm)
                  SolFct.push_back(_data->theFct[k]);
               }
            }
            if (get_err==2) {
               for (int i=0; i<size; ++i) {
                  int k = _data->FctLabel[analytic_f[i]];
                  CHK_MSGR(k==0,_pr+"end>","Undefined function "+analytic_f[i])
                  SolFct.push_back(_data->theFct[k]);
               }
            }
            *_rita->ofh << "  end" << endl;
            log.pde = false;
            set();
            return 0;

         case 106:
            _rita->setEcho(nb_args);
            break;

         default:
            DEFAULT_KW(_rita)
      }
   }
   return 0;
}


int pde::set_nls(string nls)
{
   static const vector<string> lnl {"bisection","regula-falsi","picard","secant","newton"};
   bool found = false;
   for (auto const& l: lnl) {
      if (nls==l) { found = true; break; }
   }
   CHK_MSGR(found==false,_pr+"ls>","Unknown nonlinear iterative solver: "+nls)
   return 0;
}


int pde::set_ls(string ls,
                string prec)
{
   static const vector<string> lls {"direct","cg","cgs","bicg","bicg-stab","gmres"};
   static const vector<string> lpr {"ident","diag","dilu","ilu","ssor"};
   bool found = false;
   for (auto const& l: lls) {
      if (ls==l) { found = true; break; }
   }
   CHK_MSGR(found==false,_pr+"pde>ls>","Unknown linear solver: "+ls)
   found = false;
   for (auto const& l: lpr) {
      if (prec==l) { found = true; break; }
   }
   CHK_MSGR(!found,_pr+"ls>","Unknown linear preconditioner: "+prec);
   return 0;
}


int pde::setSpD(string spd)
{
   Sdm = pde_sdm[spd];
   if (Sdm==FE_P1 || Sdm==FE_P2 || Sdm==FE_Q1 || Sdm==FV || Sdm==DG) {
      CHK_MSGR(_data->nb_meshes==0,_pr,"No defined mesh.")
      theMesh = _data->theMesh[_data->iMesh];
      _dim = theMesh->getDim();
   }
   else if (Sdm==FD) {
      CHK_MSGR(_data->nb_grids==0,_pr,"No defined grid.")
      theGrid = _data->theGrid[_data->iGrid];
      _dim = theGrid->getDim();
   }
   return 0;
}


void pde::setVectors()
{
   nb_vectors = 1;
   fd[0].fn = _vect[0];
   fd[0].nb_dof = 1;
   fd[0].ds = data::DataSize::NODES;

   switch (ieq) {

      case LINEAR_PDE:
         break;

      case LAPLACE:
         break;

      case HEAT:
         break;

      case WAVE:
         break;

      case TRANSPORT:
         break;

      case LINEAR_ELASTICITY:
         fd[0].nb_dof = _dim;
         break;

      case TRUSS:
         fd[0].nb_dof = 2;
         break;

      case BEAM:
         if (_dim==2)
            fd[0].nb_dof = 2;
         else if (_dim==3)
            fd[0].nb_dof = 6;
         break;

      case INCOMPRESSIBLE_NAVIER_STOKES:
         nb_vectors = 2;
         fd[1].fn = _vect[1];
         fd[0].nb_dof = _dim, fd[1].nb_dof = 1;
         break;
         
      case COMPRESSIBLE_EULER:
         nb_vectors = 4;
         fd[1].fn = _var[1], fd[2].fn = _var[2], fd[3].fn = _vect[3];
         fd[0].nb_dof = _dim, fd[1].nb_dof = fd[2].nb_dof = fd[3].nb_dof = 1;
         break;
         
      case COMPRESSIBLE_NAVIER_STOKES:
         nb_vectors = 4;
         fd[1].fn = _vect[1], fd[2].fn = _vect[2], fd[3].fn = _vect[3];
         fd[0].nb_dof = _dim, fd[1].nb_dof = fd[2].nb_dof = fd[3].nb_dof = 1;
         break;

      case INCOMPRESSIBLE_POROUS_1PHASE:
         break;

      case EDDY_CURRENTS:
         if (_dim==3)
            fd[0].nb_dof = _dim;
         break;

      case MAXWELL:
         fd[0].nb_dof = _dim;
         break;

      case HELMHOLTZ:
         break;
   }
}


int pde::setSize(Vect<double>& v, data::DataSize s)
{
   nb_dof = 1;
   if (theMesh!=nullptr) {
      if (theMesh->getDOFSupport()==NODE_DOF)
         nb_dof = theMesh->getNbDOF()/theMesh->getNbNodes();
      else if (theMesh->getDOFSupport()==SIDE_DOF)
         nb_dof = theMesh->getNbDOF()/theMesh->getNbSides();
      else if (theMesh->getDOFSupport()==ELEMENT_DOF)
         nb_dof = theMesh->getNbDOF()/theMesh->getNbElements();
      else if (theMesh->getDOFSupport()==EDGE_DOF)
         nb_dof = theMesh->getNbDOF()/theMesh->getNbEdges();

      if (s==data::DataSize::NODES) {
         CHK_MSGR(theMesh->getNbNodes()==0,_pr,"Mesh has no nodes")
         v.setMesh(*theMesh,NODE_DOF,nb_dof);
      }
      else if (s==data::DataSize::ELEMENTS) {
         CHK_MSGR(theMesh->getNbElements()==0,_pr,"Mesh has no elements")
         v.setMesh(*theMesh,ELEMENT_DOF,nb_dof);
      }
      else if (s==data::DataSize::SIDES) {
         CHK_MSGR(theMesh->getNbSides()==0,_pr,"Mesh has no sides")
         v.setMesh(*theMesh,SIDE_DOF,nb_dof);
      }
      else if (s==data::DataSize::EDGES) {
         CHK_MSGR(theMesh->getNbEdges()==0,_pr,"Mesh has no edges")
         v.setMesh(*theMesh,EDGE_DOF,nb_dof);
      }
   }
   else if (theGrid!=nullptr) {
      nb_dof = theGrid->getNbDOF()/theGrid->getNbNodes();
      v.setGrid(*theGrid);
   }
   else
      MSGR(_pr,"No available mesh or grid.")
   return 0;
}


int pde::getIn()
{
   int ret=0, k=0;
   bool val_ok=false, file_ok=false, save_ok=false;
   static const vector<string> kw {"val$ue","file","save"};
   _cmd->set(kw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,_pr+"initial>","No argument given")
   for (size_t i=0; i<nb_args; ++i) {
      size_t n = _cmd->getArg("=");
      switch (n) {

         case 0:
            in_data.exp = _cmd->string_token();
            val_ok = true;
            break;

         case 1:
            in_data.in_file = _cmd->string_token();
            file_ok = true;
            break;

         case 2:
            in_data.out_file = _cmd->string_token();
            save_ok = true;
            break;

         default:
            UNKNOWN_ARG(_pr+"initial>")
      }
   }
   CHK_MSGR(!val_ok,_pr+"initial>","No value or expression given for initial condition.")
   *_rita->ofh << "  in  value=" << in_data.exp;
   in_data.ft_name = "";
   k = _data->addFunction(in_data.ft_name,_var,in_data.exp);
   FCT_NOT_DEFINED("initial>",in_data.ft_name)
   if (file_ok)
      *_rita->ofh << "  file=" << in_data.in_file;
   if (save_ok) {
      if (_verb)
         cout << "Initial condition saved in file: " << in_data.out_file << endl;
      *_rita->ofh << "  save=" << in_data.out_file;
   }
   *_rita->ofh << endl;
   in_data.size++;
   set_in = true;
   return ret;
}


int pde::setIn()
{
   CHK_MSGR(in_data.size==0,_pr+"initial>","No defined initial data")
   if (theMesh!=nullptr) {
      nb_dof = theMesh->getNbDOF()/theMesh->getNbNodes();
      u.setMesh(*theMesh,NODE_DOF,nb_dof);
   }
   Vect<double> *sol = _data->theVector[fd[0].vect];
   sol->setRegex(1);
   if (in_data.in_file.size()) {
      OFELI::IOField ffi(in_data.in_file,OFELI::IOField::IN);
      ffi.get(*sol);
   }
   if (in_data.out_file.size()) {
      OFELI::IOField ffo(in_data.out_file,OFELI::IOField::OUT);
      ffo.put(*sol);
   }
   set_u = true;
   return 0;
}


int pde::getBC()
{
   string val="";
   bool code_ok=false, val_ok=false, file_ok=false, save_ok=false;
   int code=0, nb=0, ret=0;
   static const vector<string> kw {"code","val$ue","file","save"};
   _cmd->set(kw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,_pr+"bc>","No arguments")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 0:
            DEF_PAR_R(0,_pr,code)
            code_ok = true;
            break;

         case 1:
            val = _cmd->string_token(0);
            val_ok = true;
            break;

         case 2:
            bc_data.in_file = _cmd->string_token(0);
            file_ok = true;
            break;

         case 3:
            bc_data.out_file = _cmd->string_token(0);
            save_ok = true;
            break;

         default:
            UNKNOWN_ARG(_pr+"bc>")
      }
   }
   CHK_MSGR(!code_ok,_pr+"bc>","No code given for boundary condition.")
   CHK_MSGR(!val_ok,_pr+"bc>","No value or expression given for boundary condition.")
   CHK_MSGR(code<=0,_pr+"bc>","Illegal boundary condition for a nonpositive code "+to_string(code))
   bc_data.cexp[code] = val;
   bc_data.size++;
   if (_verb>1)
      cout << "Nodes with code " << code << " have prescribed value by the expression: " << val << endl;
   *_rita->ofh << "  bc  code=" << code << "  value=" << val;
   if (file_ok)
      *_rita->ofh << "  file=" << bc_data.in_file;
   if (save_ok) {
      *_rita->ofh << "  save=" << bc_data.out_file;
      if (_verb)
         cout << "Boundary condition saved in file: " << bc_data.out_file << endl;
   }
   *_rita->ofh << endl;
   set_bc = true;
   return ret;
}


int pde::setBC()
{
   CHK_MSGR(bc_data.size==0,_pr+"bc>","No defined boundary condition")
   theMesh = _data->theMesh[_data->iMesh];
   int ret = setSize(bc,data::DataSize::NODES);
   if (ret)
      return ret;
   bc.setRegex(1);
   if (bc_data.in_file.size()) {
      OFELI::IOField ffi(bc_data.in_file,OFELI::IOField::IN);
      ffi.get(bc);
   }
   if (bc_data.out_file.size()) {
      OFELI::IOField ffo(bc_data.out_file,OFELI::IOField::OUT);
      ffo.put(bc);
   }
   return 0;
}


int pde::getSF()
{
   int code=0, nb=0, ret=0;
   string val="";
   bool code_ok=false, val_ok=false, file_ok=false, save_ok=false;
   static const vector<string> kw {"code","val$ue","file","save"};
   _cmd->set(kw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,_pr+"sf>","No arguments")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 0:
            ret = _data->getPar(0,_pr,code);
            code_ok = true;
            break;

         case 1:
            val = _cmd->string_token(0);
            val_ok = true;
            break;

         case 2:
            sf_data.in_file = _cmd->string_token(0);
            file_ok = true;
            break;

         case 3:
            sf_data.out_file = _cmd->string_token(0);
            save_ok = true;
            break;

         default:
            UNKNOWN_ARG(_pr+"sf>")
      }
   }
   CHK_MSGR(!code_ok,_pr+"sf>","No code given.")
   CHK_MSGR(!val_ok,_pr+"sf>","No value or expression given for surface force.")
   CHK_MSGR(code<=0,_pr+"sf>","Illegal value of code: "+to_string(code))
   sf_data.cexp[code] = val;
   sf_data.size++;
   *_rita->ofh << "  sf  code=" << code << "  value=" << val;
   if (file_ok)
      *_rita->ofh << "  file=" << sf_data.in_file;
   if (save_ok)
      *_rita->ofh << "  save=" << sf_data.out_file;
   *_rita->ofh << endl;
   set_sf = true;
   return ret;
}


int pde::setSF()
{
   CHK_MSGR(sf_data.size==0,_pr+"sf>","No defined boundary force")
   theMesh = _data->theMesh[_data->iMesh];
   int ret = setSize(sf,data::DataSize::SIDES);
   if (ret)
      return ret;
   sf.setRegex(1);
   if (sf_data.in_file.size()) {
      OFELI::IOField ffi(sf_data.in_file,OFELI::IOField::IN);
      ffi.get(sf);
   }
   if (sf_data.out_file.size()) {
      OFELI::IOField ffo(sf_data.out_file,OFELI::IOField::OUT);
      ffo.put(sf);
   }
   return 0;
}


int pde::getBF()
{
   bool val_ok=false, file_ok=false, save_ok=false;
   static const vector<string> kw {"val$ue","file","save"};
   int nb=0;
   _cmd->set(kw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,_pr+"bf>","No arguments")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 0:
            bf_data.exp = _cmd->string_token(0);
            val_ok = true;
            break;

         case 1:
            bf_data.in_file = _cmd->string_token(0);
            file_ok = true;
            break;

         case 2:
            bf_data.out_file = _cmd->string_token(0);
            save_ok = true;
            break;

         default:
            UNKNOWN_ARG(_pr+"source>")
      }
   }
   CHK_MSGR(!val_ok,_pr+"source>","No value or expression given for source.")
   *_rita->ofh << "  source value=" << bf_data.exp;
   bf_data.ft_name = "";
   int k = _data->addFunction(bf_data.ft_name,_var,bf_data.exp);
   FCT_NOT_DEFINED("source>",bf_data.ft_name)
   if (file_ok)
      *_rita->ofh << " file=" << bf_data.in_file;
   if (save_ok)
      *_rita->ofh << " save=" << bf_data.out_file;
   *_rita->ofh << endl;
   bf_data.size++;
   set_bf = true;
   return 0;
}


int pde::setBF()
{
   CHK_MSGR(bf_data.size==0,_pr+"bf>","No defined body force")
   theMesh = _data->theMesh[_data->iMesh];
   int ret = setSize(bf,data::DataSize::NODES);
   if (ret)
      return ret;
   bf.setRegex(1);
   if (bf_data.in_file.size()) {
      IOField ffi(bf_data.in_file,IOField::IN);
      ffi.get(bf);
   }
   if (bf_data.out_file.size()) {
      IOField ffo(bf_data.out_file,IOField::OUT);
      ffo.put(bf);
   }
   return 0;
}


int pde::setNodeBC(int code, string exp, double t, Vect<double>& v)
{
   for (int i=1; i<=nb_dof; ++i)
      v.setNodeBC(code,exp,i);
   return 0;
}


int pde::setEq()
{
   setVectors();
   switch (ieq) {

//    Generic Linear PDE
      case LINEAR_PDE:
         switch (_dim) {
  
            case 1:
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",pde)

            case 2:
               CHK_MSG_PDE(Sdm!=FE_P1 && !axi,"Approximation of equation not implemented in rita.",spd)

            case 3:
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)
         }
         break;

//    Laplace equation
      case LAPLACE:
         switch (_dim) {
  
            case 1:
               CHK_MSG_PDE(Sdm!=FE_P1 && Sdm!=FE_P2,"Approximation of equation not implemented in rita.",pde)

            case 2:
               CHK_MSG_PDE(Sdm!=FE_P1 && !axi,"Approximation of equation not implemented in rita.",spd)

            case 3:
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)
         }
         break;

//    Heat equation
      case HEAT:
         switch (_dim) {

            case 1:
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",pde)

            case 2:
               CHK_MSG_PDE(Sdm!=FE_P1 && Sdm!=FE_P2,"Approximation of equation not implemented in rita.",spd)

            case 3:
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)
         }
         break;

//    Wave equation
      case WAVE:
         switch (_dim) {

            case 1:
               CHK_MSG_PDE(Sdm==FE_P1,"Equation not implemented in rita.",pde)
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)

            case 2:
               CHK_MSG_PDE(Sdm==FE_P1,"Equation not implemented in rita.",pde)
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)

            case 3:
               CHK_MSG_PDE(Sdm==FE_P1,"Equation not implemented in rita.",pde)
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)
         }
         break;

//    Linear transport equation
      case TRANSPORT:
         switch (_dim) {

            case 1:
               CHK_MSG_PDE(Sdm==FE_P1,"Equation not implemented in rita.",pde)
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)

            case 2:
               CHK_MSG_PDE(Sdm==FE_P1,"Equation not implemented in rita.",pde)
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)

            case 3:
               CHK_MSG_PDE(Sdm==FE_P1,"Equation not implemented in rita.",pde)
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)
         }
         break;

      case LINEAR_ELASTICITY:
         switch (_dim) {

            case 1:
               MSG_PDE("Equation not implemented in rita.",pde)

            case 2:
               CHK_MSG_PDE(Sdm!=FE_P1 && Sdm!=FE_Q1,"Approximation of equation not implemented in rita.",pde)

            case 3:
               CHK_MSG_PDE(Sdm!=FE_P1,"Approximation of equation not implemented in rita.",spd)
         }
         break;

      case TRUSS:
         switch (_dim) {

            case 1:
               MSG_PDE("Equation not implemented in rita.",pde)

            case 2:
               CHK_MSG_PDE(Sdm!=FE_P1,"Only 2-D P1 finite element is available for bar equation.",pde)

            case 3:
               MSG_PDE("Equation not implemented in rita.",pde)
         }
         break;
         
      case BEAM:
         switch (_dim) {

            case 1:
               CHK_MSG_PDE(Sdm==FE_P1,"No implementation for 1-D P1 finite element available for beam equation.",pde)

            case 2:
               CHK_MSG_PDE(Sdm==FE_P1,"No implementation for 2-D P1 finite element available for beam equation.",pde)
               CHK_MSG_PDE(Sdm==FE_P2,"No implementation for 2-D P2 finite element available for beam equation.",pde)

            case 3:
               CHK_MSG_PDE(Sdm==FE_P1,"No implementation for 3-D P1 finite element available for beam equation.",pde)
         }
         break;

      case INCOMPRESSIBLE_NAVIER_STOKES:
         switch (_dim) {

            case 1:
               CHK_MSG_PDE(Sdm==FE_P1,"No implementation for 1-D available for"
                                      " incompressible Navier-Stokes equations.",pde)

            case 2:
               CHK_MSG_PDE(Sdm!=FE_P1,"Only P1 finite element is implemented is available for"
                                      " incompressible Navier-Stokes equations.",pde)

            case 3:
               MSG_PDE("Equation not implemented in rita.",pde)
         }
         break;
         
      case COMPRESSIBLE_EULER:
         switch (_dim) {

            case 1:
               CHK_MSG_PDE(Sdm==FE_P1,"Space discretization method not implemented for this PDE.",pde)

            case 2:
               CHK_MSG_PDE(Sdm==FE_P1,"Space discretization method not implemented for this PDE.",pde)
               CHK_MSG_PDE(Sdm==FE_P2,"Space discretization method not implemented for this PDE.",pde)

            case 3:
               CHK_MSG_PDE(Sdm==FE_P1,"Space discretization method not implemented for this PDE.",pde)
         }
         break;

      case INCOMPRESSIBLE_POROUS_1PHASE:
         switch (_dim) {

            case 1:
               CHK_MSG_PDE(Sdm==FE_P1,"Space discretization method not implemented for this PDE.",pde)

            case 2:
               CHK_MSG_PDE(Sdm==FE_P1,"Space discretization method not implemented for this PDE.",pde)
               CHK_MSG_PDE(Sdm==FE_P2,"Space discretization method not implemented for this PDE.",pde)

            case 3:
               CHK_MSG_PDE(Sdm==FE_P1,"Space discretization method not implemented for this PDE.",pde)
         }
         break;

      default:
         log.pde = true;
         _rita->msg(_pr,"Equation not implemented in rita.");
         break;
   }
   return 0;
}


int pde::setCoef()
{
   string pr=_pr+"coef>";
   static const string H = "Command: coef [c00=x] [c10=x] [c01=x] [c20=x] [c02=x]\n"
                           "              [rho=x] [Cp=x] [kappa=x] [Mu=x] [sigma=x] [mu=x] [epsilon=x] [omega=x]\n"
                           "              [beta=x] [v=x] [young=x] [poisson=x]\n\n";
   const static vector<string> kw {"c00","c10","c01","c20","c02","rho","density","Cp","specific-heat","kappa","thermal-conductivity",
                                   "Mu","magnetic-permeability","sigma","electric-conductivity","mu","viscosity","epsilon",
                                   "electric-permittivity","omega","angular-frequency","beta","thermal-dilatation","v",
                                   "velocity","young","poisson"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(pr)
   int nb=0;
   for (size_t i=0; i<nb_args; ++i) {
      int key = _cmd->getArgs(nb);
      NO_VALUE_ARG(pr)
      switch (key) {

         case 100:
            cout << "Available arguments: " << PDE_Coef_help << endl;
            return 0;

         case 101:
            cout << PDE_Coef_Help << endl;
            return 0;

         case  0:
            _c00 = _cmd->string_token(0);
            _c00_set = true;
            break;

         case  1:
            _c10 = _cmd->string_token(0);
            _c10_set = true;
            break;

         case  2:
            _c01 = _cmd->string_token(0);
            _c01_set = true;
            break;

         case  3:
            _c20 = _cmd->string_token(0);
            _c20_set = true;
            break;

         case  4:
            _c02 = _cmd->string_token(0);
            _c02_set = true;
            break;

         case  5:
         case  6:
            _rho_exp = _cmd->string_token(0);
            _rho_set = true;
            break;

         case  7:
         case  8:
            _Cp_exp = _cmd->string_token(0);
            _Cp_set = true;
            break;

         case  9:
         case 10:
            _kappa_exp = _cmd->string_token(0);
            _kappa_set = true;
            break;

         case 11:
         case 12:
            _Mu_exp = _cmd->string_token(0);
            _Mu_set = true;
            break;

         case 13:
         case 14:
            _sigma_exp = _cmd->string_token(0);
            _sigma_set = true;
            break;

         case 15:
         case 16:
            _mu_exp = _cmd->string_token(0);
            _mu_set = true;
            break;

         case 17:
         case 18:
            _epsilon_exp = _cmd->string_token(0);
            _epsilon_set = true;
            break;

         case 19:
         case 20:
            _omega_exp = _cmd->string_token(0);
            _omega_set = true;
            break;

         case 21:
         case 22:
            _beta_exp = _cmd->string_token(0);
            _beta_set = true;
            break;

         case 23:
         case 24:
            _v_exp = _cmd->string_token(0);
            _v_set = true;
            break;

         case 25:
            _young_exp = _cmd->string_token(0);
            _young_set = true;
            break;

         case 26:
            _poisson_exp = _cmd->string_token(0);
            _poisson_set = true;
            break;

         default:
            MSGR(pr,"Unknown argument: "+_cmd->getToken())
      }
   }
   if (nb_args>0) {
      CHK_MSGR_PDE(_c00_set && ieq!=LINEAR_PDE,pr,"This coefficient is for generic linear PDEs' only.")
      CHK_MSGR_PDE(_c10_set && ieq!=LINEAR_PDE,pr,"This coefficient is for generic linear PDEs' only.")
      CHK_MSGR_PDE(_c01_set && ieq!=LINEAR_PDE,pr,"This coefficient is for generic linear PDEs' only.")
      CHK_MSGR_PDE(_c20_set && ieq!=LINEAR_PDE,pr,"This coefficient is for generic linear PDEs' only.")
      CHK_MSGR_PDE(_c02_set && ieq!=LINEAR_PDE,pr,"This coefficient is for generic linear PDEs' only.")
      CHK_MSGR_PDE(_rho_set && ieq!=HEAT && ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES,pr,"This PDE doesn't need density input.")
      CHK_MSGR_PDE(_Cp_set && ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES,pr,"This PDE doesn't need specific heat input.")
      CHK_MSGR_PDE(_kappa_set && ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES,pr,"This PDE doesn't need thermal conductivity input.")
      CHK_MSGR_PDE(_mu_set && ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES,pr,"This PDE doesn't need viscosity input.")
      CHK_MSGR_PDE(_sigma_set && ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ,pr,"This PDE doesn't need electric conductivity input.")
      CHK_MSGR_PDE(_Mu_set && ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ,pr,"This PDE doesn't need viscosity input.")
      CHK_MSGR_PDE(_epsilon_set && ieq!=MAXWELL && ieq!=HELMHOLTZ,pr,"This PDE doesn't need electric permittivity input.")
      CHK_MSGR_PDE(_omega_set && ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ,pr,"This PDE doesn't need angular frequency input.")
      CHK_MSGR_PDE(_beta_set && ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES,pr,"This PDE doesn't need thermal expansion coefficient input.")
      CHK_MSGR_PDE(_v_set && ieq!=WAVE && ieq!=TRANSPORT,pr,"This PDE doesn't need velocity input.")
      CHK_MSGR_PDE(_young_set && ieq!=LINEAR_ELASTICITY && ieq!=BEAM,pr,"This PDE doesn't need Young's modulus input.")
      CHK_MSGR_PDE(_poisson_set && ieq!=LINEAR_ELASTICITY && ieq!=BEAM,pr,"This PDE doesn't need Poisson ratio input")
   }
   /*
   else {
      while (1) {
         if (_cmd->readline(sPrompt+"pde>coef> ")<0)
            continue;
         switch (key=_cmd->getKW(kw)) {

            case  0:
            case  1:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands:\n";
               cout << "rho:      Density\n";
               cout << "Cp:       Specific heat at constant pressure\n";
               cout << "kappa:    Thermal conductivity\n";
               cout << "mu:       Viscosity\n";
               cout << "sigma:    Electric conductivity\n";
               cout << "Mu:       Magnetic permeability\n";
               cout << "epsilon:  Electric permittivity\n";
               cout << "omega:    Angular frequency\n";
               cout << "beta:     Thermal dilatation coefficient\n";
               cout << "v:        Velocity\n";
               cout << "young:    Young modulus\n";
               cout << "poisson:  Poisson ratio\n";
               cout << "end or <: go back to higher level" << endl;
               break;

            case  2:
               _rita->setConfigure();
               break;

            case  3:
            case  4:
               if (ieq!=HEAT && ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need density input" << endl;
                  *_ofl << "In "+sPrompt+"pde>coef>rho>: This PDE doesn't need density input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for density as function of x,y,z,t.")) {
                  *_ofl << "In "+sPrompt"+pde>coef>rho>: Missing regular expression for density" << endl;
                  break;
               }
               ret = _cmd->get(_rho_exp);
               if (!ret) {
                  **_rita->ofh << "    rho " << _rho_exp << endl;
                  _rho_set = true;
	       }
               break;

            case  5:
            case  6:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need specific heat input" << endl;
                  *_ofl << "In "+sPrompt+"pde>coef>Cp>: This PDE doesn't need specific heat input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for specific heat as function of x,y,z,t.")) {
                  *_ofl << "In "+sPrompt+"pde>coef>Cp>: Missing regular expression for specific heat" << endl;
                  break;
               }
               ret = _cmd->get(_Cp_exp);
               if (!ret)
                  **_rita->ofh << "    Cp " << _Cp_exp << endl;
               break;

            case  7:
            case  8:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need thermal conductivity input" << endl;
                  *_ofl << "In rita>pde>coef>kappa>: This PDE doesn't need thermal conductivity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for thermal conductivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>Cp>: Missing regular expression for thermal conductivity" << endl;
                  break;
               }
               ret = _cmd->get(_kappa_exp);
               if (!ret) {
                  **_rita->ofh << "    kappa " << _kappa_exp << endl;
                  _kappa_set = true;
               }
               break;

            case  9:
            case 10:
               if (ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need viscosity input" << endl;
                  *_ofl << "In rita>pde>coef>mu>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for viscosity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>mu>: Missing regular expression for viscosity" << endl;
                  break;
               }
               ret = _cmd->get(_mu_exp);
               if (!ret) {
                  **_rita->ofh << "    mu " << _mu_exp << endl;
                  _mu_set = true;
               }
               break;

            case 11:
            case 12:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need electric conductivity input" << endl;
                  *_ofl << "In rita>pde>coef>sigma>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for electric conductivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>sigma>: Missing regular expression for electric conductivity" << endl;
                  break;
               }
               ret = _cmd->get(_sigma_exp);
               if (!ret) {
                  **_rita->ofh << "    sigma " << _sigma_exp << endl;
                  _sigma_set = true;
               }
               break;

            case 13:
            case 14:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need magnetic permeability input" << endl;
                  *_ofl << "In rita>pde>coef>Mu>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for magnetic permeability as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>Mu>: Missing regular expression for magnetic permeability" << endl;
                  break;
               }
               ret = _cmd->get(_Mu_exp);
               if (!ret) {
                  **_rita->ofh << "    Mu " << _Mu_exp << endl;
                  _Mu_set = true;
               }
               break;

            case 15:
            case 16:
               if (ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need electric permittivity input" << endl;
                  *_ofl << "In rita>pde>coef>epsilon>: This PDE doesn't need electric permittivity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for electric permittivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>epsilon>: Missing regular expression for electric permittivity" << endl;
                  break;
               }
               ret = _cmd->get(_epsilon_exp);
               if (!ret) {
                  **_rita->ofh << "    epsilon " << _epsilon_exp << endl;
                  _epsilon_set = true;
               }
               break;

            case 17:
            case 18:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need angular frequency input" << endl;
                  *_ofl << "In rita>pde>coef>omega>: This PDE doesn't need angular frequency input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give value of angular frequency.")) {
                  *_ofl << "In rita>pde>coef>omega>: Missing value of angular frequency" << endl;
                  break;
               }
               ret = _cmd->get(_omega_exp);
               if (!ret) {
                  **_rita->ofh << "    omega " << _omega_exp << endl;
                  _omega_set = true;
               }
               break;

            case 19:
            case 20:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need thermal expansion coefficient input" << endl;
                  *_ofl << "In rita>pde>coef>beta>: This PDE doesn't need thermal expansion coefficient input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for thermal expansion coefficient as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>beta>: Missing regular expression for thermal expansion coefficient" << endl;
                  break;
               }
               ret = _cmd->get(_beta_exp);
               if (!ret) {
                  **_rita->ofh << "    beta " << _beta_exp << endl;
                  _beta_set = true;
               }
               break;

            case 21:
            case 22:
               if (ieq!=WAVE && ieq!=TRANSPORT) {
                  cout << "Error: This PDE doesn't need velocity input" << endl;
                  *_ofl << "In rita>pde>coef>v>: This PDE doesn't need velocity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for velocity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>v>: Missing regular expression for velocity" << endl;
                  break;
               }
               ret = _cmd->get(_v_exp);
               if (!ret) {
                  **_rita->ofh << "    v " << _v_exp << endl;
                  _v_set = true;
               }
               break;

            case 23:
               if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
                  cout << "Error: This PDE doesn't need Young modulus input" << endl;
                  *_ofl << "In rita>pde>coef>young>: This PDE doesn't need Young's modulus input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for Young's modulus as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>young>: Missing regular expression for Young's modulus" << endl;
                  break;
               }
               ret = _cmd->get(_young_exp);
               if (!ret) {
                  **_rita->ofh << "    young " << _young_exp << endl;
                  _young_set = true;
               }
               break;

            case 24:
               if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
                  cout << "Error: This PDE doesn't need Poisson ratio input" << endl;
                  *_ofl << "In rita>pde>coef>poisson>: This PDE doesn't need Poisson ratio input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for Poisson ratio as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>poisson>: Missing regular expression for Poisson ratio" << endl;
                  break;
               }
               ret = _cmd->get(_poisson_exp);
               if (!ret) {
                  **_rita->ofh << "    poisson " << _poisson_exp << endl;
                  _poisson_set = true;
               }
               break;

            case 25:
            case 26:
               _ret = 0;
               **_rita->ofh << "    end" << endl;
               return;

            case 27:
            case 28:
               _ret = 100;
               return;

            case 29:
               _ret = 200;
               return;

            case -2:
            case -3:
            case -4:
               break;

            default:
               cout << "Unknown Command: " << _cmd->token() << endl;
               cout << "Available commands: rho, Cp, kappa, mu, sigma, Mu, epsilon" << endl;
	       cout << "                    omega, beta, v, young, poisson, end, <" << endl;
               cout << "Global commands:    help, ?, set, quit, exit" << endl;
               *_ofl << "In rita>pde>coef>: Unknown PDE Coefficient " << _cmd->token() << endl;
               break;
         }
      }
      }*/
   return 0;
}


void pde::setPDECoef()
{
   if (_c00_set)
      theEquation->setPDECoef(OFELI::PDECoefType::C00,_c00);
   if (_c10_set)
      theEquation->setPDECoef(PDECoefType::C10,_c10);
   if (_c01_set)
      theEquation->setPDECoef(PDECoefType::C01,_c01);
   if (_c20_set)
      theEquation->setPDECoef(PDECoefType::C20,_c20);
   if (_c02_set)
      theEquation->setPDECoef(PDECoefType::C02,_c02);
}


void pde::set()
{
   switch (ieq) {

      case LINEAR_PDE:
         switch (_dim) {
  
            case 1:
               theEquation = new LinearPDE1D(*theMesh);
               setPDECoef();
               break;

            case 2:
               theEquation = new LinearPDE2D(*theMesh);
               setPDECoef();
               break;

            case 3:
               theEquation = new LinearPDE3D(*theMesh);
               setPDECoef();
               break;

         }
         break;

      case LAPLACE:
         switch (_dim) {
  
            case 1:
               if (Sdm==FE_P1)
                  theEquation = new Laplace1DL2(*theMesh);
               else if (Sdm==FE_P2)
                  theEquation = new Laplace1DL3(*theMesh);
               break;

            case 2:
               if (Sdm==FE_P1)
                  theEquation = new Laplace2DT3(*theMesh);
               break;

            case 3:
               if (Sdm==FE_P1)
                  theEquation = new Laplace3DT4(*theMesh);
               break;
         }
         break;

      case HEAT:
         switch (_dim) {

            case 1:
               if (Sdm==FE_P1) {
                  theEquation = new DC1DL2(*theMesh);
                  theEquation->setTerms(PDE_Terms::LUMPED_CAPACITY);
                  theEquation->setTerms(PDE_Terms::DIFFUSION);
               }
               break;

            case 2:
               if (Sdm==FE_P1 && !axi)
                  theEquation = new DC2DT3(*theMesh);
               if (Sdm==FE_P1 && axi)
                  theEquation = new DC3DAT3(*theMesh);
               else if (Sdm==FE_P2 && !axi)
                  theEquation = new DC2DT6(*theMesh);
               theEquation->setTerms(PDE_Terms::LUMPED_CAPACITY);
               theEquation->setTerms(PDE_Terms::DIFFUSION);
               break;

            case 3:
               if (Sdm==FE_P1) {
                  theEquation = new DC3DT4(*theMesh);
                  theEquation->setTerms(PDE_Terms::LUMPED_CAPACITY);
                  theEquation->setTerms(PDE_Terms::DIFFUSION);
               }
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_Cp_set)
            theEquation->set_Cp(_Cp_exp);
         if (_kappa_set)
            theEquation->set_kappa(_kappa_exp);
         break;


      case LINEAR_ELASTICITY:
         switch (_dim) {

            case 2:
               if (Sdm==FE_P1)
                  theEquation = new Elas2DT3(*theMesh);
               else if (Sdm==FE_Q1)
                  theEquation = new Elas2DQ4(*theMesh);
               break;

            case 3:
               if (Sdm==FE_P1)
                  theEquation = new Elas3DT4(*theMesh);
               else if (Sdm==FE_Q1)
                  theEquation = new Elas3DH8(*theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_young_set)
            theEquation->set_young(_young_exp);
         if (_poisson_set)
            theEquation->set_poisson(_poisson_exp);
         break;

      case TRUSS:
         switch (_dim) {

            case 2:
               if (Sdm==FE_P1)
                  theEquation = new Bar2DL2(*theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_young_set)
            theEquation->set_young(_young_exp);
         break;

      case INCOMPRESSIBLE_NAVIER_STOKES:
         switch (_dim) {

            case 2:
               if (Sdm==FE_P1)
                  theEquation = new TINS2DT3S(*theMesh);
               break;

            case 3:
               if (Sdm==FE_P1)
                  theEquation = new TINS3DT4S(*theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_mu_set)
            theEquation->set_mu(_mu_exp);
         if (_beta_set)
            theEquation->set_beta(_beta_exp);
         break;
   }
   if (theEquation->SolverIsSet()==false)
      ls = CG_SOLVER, prec = DILU_PREC;
}


void pde::print(ostream& s) const
{
   static map<int,string> nm = {{LINEAR_PDE,"Generic Linear PDE"},
                                {LAPLACE,"Laplace equation"},
                                {HEAT,"Heat equation"},
                                {WAVE,"Wave equation"},
                                {TRANSPORT,"Linear transport equation"},
                                {LINEAR_ELASTICITY,"Linear Elasticity equations"},
                                {TRUSS,"2-D Elastic truss equations"},
                                {BEAM,"3-D Elastic beam equations"},
                                {INCOMPRESSIBLE_NAVIER_STOKES,"Incompressible Navier-Stokes equations"},
                                {COMPRESSIBLE_EULER,"Compressible Euler equations"},
                                {INCOMPRESSIBLE_POROUS_1PHASE,"Incompressible One-phase Porous media equation"}};
   static map<int,string> sd = {{FD,"Finite Differences"},
                                {FE_P1,"P1 Finite Elements"},
                                {FE_P2,"P2 Finite Elements"},
                                {FE_Q1,"Q1 Finite Elements"},
                                {FV,"Finite Volumes"},
                                {DG,"Discontinuous Galerkin"}};
   if (log.fail()) {
      s << "PDE improperly or not defined !" << endl;
      return;
   }
   s << "PDE Name: " << name << endl;
   s << "PDE id.: " << eq << endl;
   s << "PDE: " << nm[ieq] << endl;
   s << "PDE unknown vector(s): ";
   for (int i=0; i<nb_vectors-1; ++i) {
      if (_data->dn[fd[i].fn].active==false)
         s << "WARNING: Vector " << fd[i].fn << " has been removed." << endl;
   }
   for (int i=0; i<nb_vectors-1; ++i)
      s << fd[i].fn << ", ";
   s << fd[nb_vectors-1].fn << endl;
   s << "PDE linear solver: " << lsolv << endl;
   if (lsolv!="direct")
      s << "PDE linear preconditioner: " << lprec << endl;
   s << "Space dimension: " << _dim << endl;
   s << "Space discretization: " << sd[Sdm] << endl;
   if (solved==0) {
      s << "PDE unsolved !" << endl;
      return;
   }
   s << "Solution of PDE stored in vector(s): "; 
   s << fd[0].fn;
   for (int i=1; i<nb_vectors; ++i)
      s << ", " << fd[i].fn;
   s << endl;
   if (solved==false)
      cout << "PDE unsolved." << endl;
}


ostream& operator<<(ostream& s, const pde& e)
{
   e.print(s);
   return s;
}

} /* namespace RITA */
