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

                         Implementation of class 'data'

  ==============================================================================*/

#include "data.h"
#include "configure.h"
#include "rita.h"
#include "ls.h"
#include "ae.h"
#include "ode.h"
#include "pde.h"
#include "linear_algebra/Matrix.h"
#include "io/IOField.h"
#include "calc.h"
#include "optim.h"
#include "eigen.h"
#include "helps.h"

using std::cout;
using std::endl;
using std::map;

namespace RITA {

data::data(rita*      r,
           cmd*       command,
           configure* config)
     : _rita(r), _size(0), _verb(1), _configure(config), _cmd(command),
       _theMesh(nullptr), _theTab(nullptr), _theGrid(nullptr), _theFct(nullptr), 
       _theVector(nullptr), _theMatrix(nullptr), _u(nullptr), _hu(nullptr)
{
   nb_vectors = nb_params = nb_fcts = nb_meshes = nb_tabs = nb_grids = nb_matrices = nb_hvectors = 0;
   nb_ls = nb_pde = nb_ode = nb_ae = nb_int = nb_eig = nb_opt = nb_eq = 0;
   theParam.push_back(0.), theVector.push_back(nullptr), VectorTime.push_back(0.), theFct.push_back(nullptr);
   theMatrix.push_back(nullptr), theMesh.push_back(nullptr), theGrid.push_back(nullptr), theHVector.push_back(nullptr);
   theLS.push_back(nullptr), theAE.push_back(nullptr), theODE.push_back(nullptr), thePDE.push_back(nullptr);
   theTab.push_back(nullptr), theEig.push_back(nullptr), theOpt.push_back(nullptr);
   FctName.push_back(""), MeshName.push_back(""), GridName.push_back(""), TabName.push_back(""), VectorName.push_back("");
   MatrixName.push_back(""), HVectorName.push_back(""), OptName.push_back(""), EigName.push_back(""), ParamName.push_back("");
   LSName.push_back(""), AEName.push_back(""), ODEName.push_back(""), PDEName.push_back("");
   iParam = iVector = iHVector = iMatrix = iGrid = iMesh = iTab = iFct = iLS = iAE = iODE = iPDE = iOpt = iEig = 0;
   VectorSizeType.push_back(DataSize::GIVEN_SIZE);
   VectorEquation.push_back(0);
   eq_type.push_back(eqType::AE);
   eq_type.push_back(eqType::LS);
   nb_dof.push_back(0);
   eqq.push_back(0);
   _ret = 0;
   ok = false;
}


int data::setDataExt(int key)
{
   int ret = 0;
   string td = "", fn="", name="", str1="", str2="", old_name="", new_name="";
   static const string H = "parameters, grids, meshes, vectors, hvectors, functions, tabulations, matrices";
   switch (key) {

      case 200:
         ret = setGrid();
         break;

      case 201:
         _rita->_mesh = new mesh(_rita,_cmd,_configure);
         ret = _rita->_mesh->run();
         if (ret) {
            delete _rita->_mesh;
            _rita->_mesh = nullptr;
         }
         break;

      case 202:
         ret = setVector();
         break;

      case 203:
         ret = setTab();
         break;

      case 204:
         ret = setFunction();
         break;

      case 205:
         ret = setMatrix();
         break;

      case 206:
         ret = setSample();
         break;

      case 207:
         ret = save();
         break;

      case 208:
      case 209:
         CHK_MSG1B(_cmd->setNbArg(1,"Data name to delete."),"","Data name to delete.","",1)
         _cmd->get(td);
         ret = remove(td);
         if (!ret)
            *_rita->ofh << "remove " << td << endl;
         break;

      case 210:
         _cmd->get(name);
         _cmd->get(td);
         *_rita->ofh << "desc " << name << " " << td << endl;
         setDesc(name,td);
         break;

      case 211:
         _cmd->get(str1);
         _cmd->get(str2);
         CHK_MSGR(VectorLabel[str1]==0,"","Vector "+str1+" not found.")
         {
            Vect<double> *v = theVector[VectorLabel[str1]];
            ret = setHistory(str2,v->size());
            vect_hist[str1] = str2;
            if (!ret)
               *_rita->ofh << "history " << str1 << " " << str2 << endl;
         }
         break;

      case 212:
         *_rita->ofh << "data" << endl;
         Summary();
         break;

      case 213:
         CHK_MSG1B(_cmd->setNbArg(1,"Type of data to list."),"","Missing data type to list.\nAvailable data: "+H,"",1)
         _cmd->get(td);
         if (td.substr(0,5)=="param") {
            *_rita->ofh << "list param" << endl;
            ListParams(1);
         }
         else if (td.substr(0,4)=="grid") {
            *_rita->ofh << "list grid" << endl;
            ListGrids(1);
         }
         else if (td.substr(0,4)=="mesh") {
            *_rita->ofh << "list mesh" << endl;
            ListMeshes(1);
         }
         else if (td.substr(0,4)=="vect") {
            *_rita->ofh << "list vect" << endl;
            ListVectors(1);
         }
         else if (td.substr(0,5)=="hvect") {
            *_rita->ofh << "list hvect" << endl;
            ListHVectors(1);
         }
         else if (td.substr(0,5)=="funct") {
            *_rita->ofh << "list funct" << endl;
            ListFunctions(1);
         }
         else if (td.substr(0,3)=="tab") {
            *_rita->ofh << "list tab" << endl;
            ListTabs(1);
         }
         else if (td.substr(0,5)=="matri") {
            *_rita->ofh << "list matrices" << endl;
            ListMatrices(1);
         }
         else if (td.substr(0,2)=="ls") {
            *_rita->ofh << "list ls" << endl;
            ListLS(1);
         }
         else if (td.substr(0,2)=="ae" || td.substr(0,9)=="algebraic") {
            *_rita->ofh << "list algebraic" << endl;
            ListAE(1);
         }
         else if (td.substr(0,3)=="ode") {
            *_rita->ofh << "list ode" << endl;
            ListODE(1);
         }
         else if (td.substr(0,3)=="pde") {
            *_rita->ofh << "list pde" << endl;
            ListPDE(1);
         }
         else if (td.substr(0,3)=="opt") {
            *_rita->ofh << "list opt" << endl;
            ListOpt(1);
         }
         else if (td.substr(0,3)=="eig") {
            *_rita->ofh << "list eigen" << endl;
            ListEig(1);
         }
         else
            _rita->msg("list>","Unknown data type "+td);
         ret = 0;
         break;

      case 214:
      case 215:
         CHK_MSG1B(_cmd->setNbArg(1,"Data name to be given.",1),"","Missing data to display.","",1)
         if (!_cmd->get(fn)) {
            print(fn);
            *_rita->ofh << "print " << fn << endl;
         }
         break;

      case 216:
         _cmd->get(old_name);
         _cmd->get(new_name);
         ret = rename(old_name,new_name);
         if (!ret)
            *_rita->ofh << "rename " << old_name << " " << new_name << endl;
         break;

      case 217:
         plot();
         break;

      case 218:
         _rita->_calc->run();
         break;
   }
   return ret;
}


int data::rename(const string& old_name,
                 const string& new_name)
{
   if (dn[new_name].active) {
      cout << "Entity " << new_name << " is already used." << endl;
      return 1;
   }
   if (!dn[old_name].active) {
      cout << "Entity " << old_name << " has been removed." << endl;
      return 1;
   }
   dn[old_name].active = false;
   dn[new_name].active = true;
   DataType d = getType(old_name);
   dn[new_name].dt = d;
   switch (d) {

      case DataType::PARAM:
         ParamName[ParamLabel[old_name]] = new_name;
         ParamLabel[new_name] = dn[new_name].i = dn[old_name].i;
         ParamLabel[old_name] = 0;
         break;

      case DataType::VECTOR:
         VectorName[VectorLabel[old_name]] = new_name;
         VectorLabel[new_name] = dn[new_name].i = dn[old_name].i;
         VectorLabel[old_name] = 0;
         theVector[VectorLabel[new_name]]->setName(new_name);
         break;

      case DataType::HVECTOR:
         HVectorName[HVectorLabel[old_name]] = new_name;
         HVectorLabel[new_name] = dn[new_name].i = dn[old_name].i;
         HVectorLabel[old_name] = 0;
         break;

      case DataType::MATRIX:
         MatrixName[VectorLabel[old_name]] = new_name;
         MatrixLabel[new_name] = dn[new_name].i = dn[old_name].i;
         MatrixLabel[old_name] = 0;
         theMatrix[MatrixLabel[new_name]]->setName(new_name);
         break;

      case DataType::GRID:
         GridName[GridLabel[old_name]] = new_name;
         GridLabel[new_name] = dn[new_name].i = dn[old_name].i;
         GridLabel[old_name] = 0;
         break;

      case DataType::MESH:
         MeshName[MeshLabel[old_name]] = new_name;
         MeshLabel[new_name] = dn[new_name].i = dn[old_name].i;
         MeshLabel[old_name] = 0;
         break;

      case DataType::TAB:
         TabName[TabLabel[old_name]] = new_name;
         TabLabel[new_name] = dn[new_name].i = dn[old_name].i;
         TabLabel[old_name] = 0;
         break;

      case DataType::FCT:
         FctName[FctLabel[old_name]] = new_name;
         FctLabel[new_name] = dn[new_name].i = dn[old_name].i;
         FctLabel[old_name] = 0;
         break;

      case DataType::LS:
         LSName[FctLabel[old_name]] = new_name;
         LSLabel[new_name] = dn[new_name].i = dn[old_name].i;
         LSLabel[old_name] = 0;
         break;

      case DataType::AE:
         AEName[FctLabel[old_name]] = new_name;
         AELabel[new_name] = dn[new_name].i = dn[old_name].i;
         AELabel[old_name] = 0;
         break;

      case DataType::ODE:
         ODEName[FctLabel[old_name]] = new_name;
         dn[new_name].dt = d;
         ODELabel[new_name] = dn[new_name].i = dn[old_name].i;
         break;

      case DataType::PDE:
         PDEName[FctLabel[old_name]] = new_name;
         PDELabel[new_name] = dn[new_name].i = dn[old_name].i;
         PDELabel[old_name] = 0;
         break;

      case DataType::OPTIM:
         OptName[FctLabel[old_name]] = new_name;
         OptLabel[new_name] = dn[new_name].i = dn[old_name].i;
         OptLabel[old_name] = 0;
         break;

      case DataType::EIGEN:
         EigName[FctLabel[old_name]] = new_name;
         EigLabel[new_name] = dn[new_name].i = dn[old_name].i;
         EigLabel[old_name] = 0;
         break;

      default:
         MSGR("","No data named "+old_name+" found")
   }
   Desc[new_name] = Desc[old_name];
   return 0;
}


int data::add2History(const string& s, OFELI::Vect<double>& v, double t)
{
   CHK_MSGR(HVectorLabel[s]==0,"history>","History Vector "+s+" not found.")
   int k = HVectorLabel.at(s);
   _hu = theHVector[k];
   CHK_MSGR(int(v.size()) != _hu->size,"history>",
            "Vector does not have the right size to be stored in history vector "+s+".")
   theHVector[k]->set(v,t);
   return 0;
}


int data::setHistory(const string& s, size_t n)
{
   if (HVectorLabel[s]==0) {
      _hu = new HVect;
      HVectorName.push_back(s);
      nb_hvectors++;
      theHVector.push_back(_hu);
      dn[s].active = true;
      Desc[s] = "History vector containing vectors with size "+to_string(n);
   }
   else
      _hu = theHVector[HVectorLabel.at(s)];
   HVectorLabel[s] = nb_hvectors;
   _hu->setSize(n);
   dn[s].i = nb_hvectors;
   dn[s].dt = DataType::HVECTOR;
   return 0;
}


int data::setSample()
{
   static const vector<string> kw {"fun$ction","nb","x","y","z","grid","mesh","vec$tor"};
   int nb=0, nb_var=0, nx=0, ny=0, nz=0, dim=1;
   double mx=0., my=0., mz=0., Mx=1., My=1., Mz=1.;
   string fct="", gr="", vec="", ms="", name="";
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,"sample>","No argument given for command.")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Sample_help << endl;
            return 0;

         case 101:
            cout << Sample_Help << endl;
            return 0;

         case 0:
            fct = _cmd->string_token(0);
            break;

         case 1:
            DDEF_PAR_R(0,"sample>",nx)
            if (nb>1) {
               DDEF_PAR_R(1,"sample>",ny)
               dim = 2;
            }
            if (nb>2) {
               DDEF_PAR_R(2,"sample>",nz)
               dim = 3;
            }
            break;

         case 2:
            DDEF_PAR_R(0,"sample>",mx)
            DDEF_PAR_R(1,"sample>",Mx)
            break;

         case 3:
            DDEF_PAR_R(0,"sample>",my)
            DDEF_PAR_R(1,"sample>",My)
            break;

         case 4:
            DDEF_PAR_R(0,"sample>",mz)
            DDEF_PAR_R(1,"sample>",Mz)
            break;

         case 5:
            gr = _cmd->string_token(0);
            break;

         case 6:
            ms = _cmd->string_token(0);
            break;

         case 7:
            vec = _cmd->string_token(0);
            break;

         default:
            UNKNOWN_ARG("sample>")
      }
   }

   if (nb_args>0) {
      CHK_MSGR(fct=="","sample>","No function given for sampling.")
      int k = FctLabel[fct];
      CHK_MSGR(k==0,"sample>","Reference to undefined function.")
      CHK_MSGR(nx==0 && gr=="" && ms=="","sample>","A grid, mesh or a number of points must be given.")
      CHK_MSGR(ms!="" && gr!="","sample>","Mesh and grid cannot be given simultaneously.")
      CHK_MSGR(nx!=0 && gr!="","sample>","Grid and number of points cannot be given simultaneously.")
      CHK_MSGR(nx!=0 && ms!="","sample>","Mesh and number of points cannot be given simultaneously.")
      CHK_MSGR(vec=="","sample>","No vector given to contain sampling.")
      _theFct = theFct[k];
      nb_var = _theFct->nb_var;
      CHK_MSGR(nb_var<dim && nx,"sample>","Number of variables must be at least "+to_string(dim))
      *_rita->ofh << "sample function=" << fct;
      int kv = VectorLabel[vec];
      if (kv>0) {
         cout << "Vector " << vec << " is overwritten" << endl;
         _theVector = theVector[kv];
      }
      else {
         _theVector = new OFELI::Vect<double>;
         if (nx)
            addVector(_theVector,vec);
      }
      *_rita->ofh << " vector=" << vec;
      if (nx) {
         if (dim==1) {
            *_rita->ofh << " nb=" << nx;
            *_rita->ofh << " x=" << mx << "," << Mx;
            _theGrid = new OFELI::Grid(mx,Mx,nx-1);
            _theVector->setGrid(*_theGrid);
            for (int i=1; i<=nx; ++i)
               (*_theVector)(i) = (*_theFct)(_theGrid->getX(i));
         }
         if (dim==2) {
            *_rita->ofh << " nb=" << nx << "," << ny;
            *_rita->ofh << " x=" << mx << "," << Mx << " y=" << my << "," << My;
            _theGrid = new OFELI::Grid(mx,Mx,my,My,nx-1,ny-1);
            _theVector->setGrid(*_theGrid);
            for (int i=1; i<=nx; ++i)
               for (int j=1; j<=ny; ++j)
                  (*_theVector)(i,j) = (*_theFct)(_theGrid->getX(i),_theGrid->getY(j));
         }
         if (dim==3) {
            *_rita->ofh << " nb=" << nx << "," << ny << "," << nz;
            *_rita->ofh << " x=" << mx << "," << Mx << " y=" << my << "," << My << " z=" << mz << "," << Mz;
            _theGrid = new OFELI::Grid(mx,Mx,my,My,mz,Mz,nx-1,ny-1,nz-1);
            _theVector->setGrid(*_theGrid);
            for (int i=1; i<=nx; ++i)
               for (int j=1; j<=ny; ++j)
                  for (int k=1; k<=nz; ++k)
                     (*_theVector)(i,j,k) = (*_theFct)(_theGrid->getX(i),_theGrid->getY(j),_theGrid->getZ(k));
         }
         addGrid(_theGrid,name);
      }
      if (gr!="") {
         *_rita->ofh << " grid=" << gr;
         int k = GridLabel[gr];
         CHK_MSGR(k==0,"sample>","Reference to undefined grid.")
         _theGrid = theGrid[k];
         _theVector->setGrid(*_theGrid);
         if (_theGrid->getDim()==1) {
            for (size_t i=1; i<=_theGrid->getNx(); ++i)
               (*_theVector)(i) = (*_theFct)(_theGrid->getX(i));
         }
         else if (_theGrid->getDim()==2) {
            for (size_t i=1; i<=_theGrid->getNx(); ++i)
               for (size_t j=1; j<=_theGrid->getNy(); ++j)
                  (*_theVector)(i,j) = (*_theFct)(_theGrid->getX(i),_theGrid->getY(j));
         }
         else if (_theGrid->getDim()==3) {
            for (size_t i=1; i<=_theGrid->getNx(); ++i)
               for (size_t j=1; j<=_theGrid->getNy(); ++j)
                  for (size_t k=1; k<=_theGrid->getNz(); ++k)
                     (*_theVector)(i,j,k) = (*_theFct)(_theGrid->getX(i),_theGrid->getY(j),_theGrid->getZ(k));
         }
         if (kv==0)
            addVector(_theVector,vec);
         _theVector->setGrid(*_theGrid);
      }
      else if (ms!="") {
         *_rita->ofh << " mesh=" << ms;
         int k = MeshLabel[ms];
         CHK_MSGR(k==0,"sample>","Reference to undefined mesh.")
         _theMesh = theMesh[k];
         for (size_t n=1; n<=_theMesh->getNbNodes(); ++n)
            (*_theVector)(n) = (*_theFct)((*_theMesh)[n]->getCoord());
         _theVector->setMesh(*_theMesh);
      }
      _theVector->setName(vec);
      *_rita->ofh << endl;
   }
   else
      *_rita->ofh << "sample" << endl;
   return 0;
}


int data::setTab2GridVector(OFELI::Tabulation* tab,
                            string&            name)
{
   double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0;
   int nb = tab->getNbVar(1);
   OFELI::fct &f = tab->Funct[0];
   int nx = f.Np[0], ny = f.Np[1], nz = f.Np[2];
   switch (nb)
   {
      case 1:
         xmin = tab->getMinVar(1,1);
         xmax = tab->getMaxVar(1,1);
         _theGrid = new OFELI::Grid(xmin,xmax,nx);
         break;   

      case 2:
         xmin = tab->getMinVar(1,1);
         xmax = tab->getMaxVar(1,1);
         ymin = tab->getMinVar(1,2);
         ymax = tab->getMaxVar(1,2);
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,nx,ny);
         break;   

      case 3:
         xmin = tab->getMinVar(1,1);
         xmax = tab->getMaxVar(1,1);
         ymin = tab->getMinVar(1,2);
         ymax = tab->getMaxVar(1,2);
         zmin = tab->getMinVar(1,3);
         zmax = tab->getMaxVar(1,3);
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
         break;   
   }
   string gname = "";
   int k = addGrid(_theGrid,gname);
   _u = new OFELI::Vect<double>(*_theGrid);
   k = addVector(_u,name);
   return k;
}


void data::setTab2Vector(OFELI::Tabulation* tab,
                         OFELI::Vect<double>& v)
{
   int nb = tab->getNbVar(1);
   OFELI::fct &f = tab->Funct[0];
   int nx = f.Np[0], ny = f.Np[1], nz = f.Np[2];
   string name="";
//   addGridVector(name,tab->getFunctName(1),1);

   switch (nb) {
      case 1:
         v.setSize(nx);
         for (int i=1; i<=nx; ++i)
            v(i) = f.Val(i);
         break;   

      case 2:
         v.setSize(nx,ny);
         for (int i=1; i<=nx; ++i)
            for (int j=1; j<=ny; ++j)
               v(i,j) = f.Val(i,j);
         break;   

      case 3:
         v.setSize(nx,ny,nz);
         for (int i=1; i<=nx; ++i)
            for (int j=1; j<=ny; ++j)
               for (int k=1; k<=nz; ++k)
                  v(i,j,k) = f.Val(i,j,k);
         break;   
   }
}


int data::checkParam(const string& name,
                     double&       value)
{
   if (ParamLabel[name]) {
      int n = ParamLabel[name];
      value = theParam[n];
      return n;
   }
   else
      return -1;
}


int data::checkParam(const string& name,
                     int&          value)
{
   if (ParamLabel[name]) {
      int n = ParamLabel[name];
      value = theParam[n];
      return n;
   }
   else
      return -1;
}


bool data::FindName(const string& s) const
{
   for (auto const& n: Names) {
      if (n==s)
         return true; 
   }
   return false;
}


int data::addLS(ls*     e,
                string& name)
{
   int opt = 0;
   iLS = theLS.size();
   DEFAULT_LS_NAME(name);
   if (checkAlreadyUsed(name,DataType::LS))
      return 0;
   if (LSLabel[name]) {
      CHK_MSGR0(dn[name].active,"","Linear system "+name+" already defined.")
      iLS = LSLabel[name];
      theLS[iLS] = e;
      opt = 1;
   }
   Desc[name] = "";
   e->name = name;
   p2s = name, pt2s = DataType::LS;
   LSLabel[name] = dn[name].i = ++nb_ls;
   dn[name].dt = DataType::LS;
   dn[name].active = true;
   if (opt==0) {
      theLS.push_back(e);
      LSName.push_back(name);
   }
   return iLS;
}


int data::addAE(ae*     e,
                string& name)
{
   int opt = 0;
   iAE = theAE.size();
   DEFAULT_AE_NAME(name);
   if (checkAlreadyUsed(name,DataType::AE))
      return 0;
   if (AELabel[name]) {
      CHK_MSGR0(dn[name].active,"","Algebraic equation "+name+" already defined.")
      iAE = AELabel[name];
      theAE[iAE] = e;
      opt = 1;
   }
   Desc[name] = "";
   e->name = name;
   p2s = name, pt2s = DataType::AE;
   AELabel[name] = dn[name].i = ++nb_ae;
   dn[name].dt = DataType::AE;
   dn[name].active = true;
   if (opt==0) {
      theAE.push_back(e);
      AEName.push_back(name);
   }
   return iAE;
}


int data::addODE(ode*    e,
                 string& name)
{
   int opt = 0;
   iODE = theODE.size();
   DEFAULT_ODE_NAME(name);
   if (checkAlreadyUsed(name,DataType::ODE))
      return 0;
   if (ODELabel[name]) {
      CHK_MSGR0(dn[name].active,"","ODE "+name+" already defined.")
      iODE = ODELabel[name];
      theODE[iODE] = e;
      opt = 1;
   }
   Desc[name] = "";
   e->name = name;
   p2s = name, pt2s = DataType::ODE;
   ODELabel[name] = dn[name].i = ++nb_ode;
   dn[name].dt = DataType::ODE;
   dn[name].active = true;
   if (opt==0) {
      theODE.push_back(e);
      ODEName.push_back(name);
   }
   return iODE;
}


int data::addPDE(pde*    e,
                 string& name)
{
   int opt = 0;
   iPDE = thePDE.size();
   DEFAULT_PDE_NAME(name);
   if (checkAlreadyUsed(name,DataType::PDE))
      return 0;
   if (PDELabel[name]) {
      CHK_MSGR0(dn[name].active,"","PDE "+name+" already defined.")
      iPDE = PDELabel[name];
      thePDE[iPDE] = e;
      opt = 1;
   }
   Desc[name] = "";
   e->name = name;
   p2s = name, pt2s = DataType::PDE;
   PDELabel[name] = dn[name].i = ++nb_pde;
   dn[name].dt = DataType::PDE;
   dn[name].active = true;
   if (opt==0) {
      thePDE.push_back(e);
      PDEName.push_back(name);
   }
   return iPDE;
}


int data::addOpt(optim*  o,
                 string& name)
{
   int opt = 0;
   iOpt = theOpt.size();
   DEFAULT_OPT_NAME(name);
   if (checkAlreadyUsed(name,DataType::OPTIM))
      return 0;
   if (OptLabel[name]) {
      CHK_MSGR0(dn[name].active,"","Optimization problem "+name+" already defined.")
      iOpt = OptLabel[name];
      theOpt[iOpt] = o;
      opt = 1;
   }
   Desc[name] = "";
   o->name = name;
   p2s = name, pt2s = DataType::OPTIM;
   OptLabel[name] = dn[name].i = ++nb_opt;
   dn[name].dt = DataType::OPTIM;
   dn[name].active = true;
   if (opt==0) {
      theOpt.push_back(o);
      OptName.push_back(name);
   }
   return iOpt;
}


int data::addEig(eigen*  e,
                 string& name)
{
   int opt = 0;
   iEig = theEig.size();
   DEFAULT_EIG_NAME(name);
   if (checkAlreadyUsed(name,DataType::EIGEN))
      return 0;
   if (EigLabel[name]) {
      CHK_MSGR0(dn[name].active,"","Eigenvalue problem "+name+" already defined.")
      iEig = EigLabel[name];
      theEig[iEig] = e;
      opt = 1;
   }
   Desc[name] = "";
   e->name = name;
   p2s = name, pt2s = DataType::EIGEN;
   EigLabel[name] = dn[name].i = ++nb_eig;
   dn[name].dt = DataType::EIGEN;
   dn[name].active = true;
   if (opt==0) {
      theEig.push_back(e);
      EigName.push_back(name);
   }
   return iEig;
}


int data::addFunction(OFELI::Fct& f)
{
   iFct = theFct.size();
   string name = f.getName();
   DEFAULT_FCT_NAME(name);
   if (checkAlreadyUsed(name,DataType::FCT))
      return 0;
   _theFct = new funct(_rita,name);
   _theFct->set(f);
   Desc[name] = "";
   if (FctLabel[name]) {
      REDEFINED("Function ")
      iFct = FctLabel[name];
      dn[name].active  = true;
      theFct[iFct] = _theFct;
      return iFct;
   }
   FctName.push_back(name);
   theFct.push_back(_theFct);
   dn[name].active = true;
   FctName[iFct] = name;
   FctLabel[name] = dn[name].i = ++nb_fcts;
   dn[name].dt = DataType::FCT;
   return iFct;
}


int data::addFunction(string& name)
{
   iFct = theFct.size();
   DEFAULT_FCT_NAME(name);
   if (checkAlreadyUsed(name,DataType::FCT))
      return 0;
   _theFct = new funct(_rita,name);
   Desc[name] = "";
   if (FctLabel[name]) {
      REDEFINED("Function ")
      iFct = FctLabel[name];
      dn[name].active  = true;
      theFct[iFct] = _theFct;
      return iFct;
   }
   FctName.push_back(name);
   theFct.push_back(_theFct);
   dn[name].active = true;
   FctName[iFct] = name;
   FctLabel[name] = dn[name].i = ++nb_fcts;
   dn[name].dt = DataType::FCT;
   return iFct;
}


int data::addFunction(string&       name,
                      const string& var,
                      const string& exp)
{
   vector<string> v;
   v.push_back(var);
   return addFunction(name,v,exp);
}


int data::addFunction(string&               name,
                      const vector<string>& var,
                      const string&         exp,
                      int                   n)
{
   iFct = theFct.size();
   DEFAULT_FCT_NAME(name);
   if (checkAlreadyUsed(name,DataType::FCT))
      return 0;
   _theFct = new funct(_rita,name);
   if (var[0]=="t") {
      _theFct->setVar(var[0]);
      if (var.size()>1) {
         _theFct->setVar(var[1],n);
         for (size_t i=2; i<var.size(); ++i)
            _theFct->setVar(var[i],n);
      }
   }
   else {
      _theFct->setVar(var[0],n);
      for (size_t i=1; i<var.size(); ++i)
         _theFct->setVar(var[i],n);
   }
   int ret = _theFct->setExpr(exp);
   if (ret) {
      delete _theFct;
      return -1;
   }
   Desc[name] = "";
   if (FctLabel[name]) {
      REDEFINED("Function ")
      iFct = FctLabel[name];
      dn[name].active  = true;
      theFct[iFct] = _theFct;
      return iFct;
   }
   FctName.push_back(name);
   theFct.push_back(_theFct);
   dn[name].active = true;
   FctName[iFct] = name;
   FctLabel[name] = dn[name].i = ++nb_fcts;
   dn[name].dt = DataType::FCT;
   return iFct;
}


int data::addGrid(OFELI::Grid* g,
                  string&      name)
{
   iGrid = theGrid.size();
   DEFAULT_GRID_NAME(name);
   if (checkAlreadyUsed(name,DataType::GRID))
      return 0;
   Desc[name] = "";
   if (GridLabel[name]) {
      REDEFINED("Grid ")
      iGrid = GridLabel[name];
      dn[name].active = true;
      theGrid[iGrid] = g;
      return iGrid;
   }
   theGrid.push_back(g);
   GridName.push_back(name);
   dn[name].active = true;
   GridLabel[name] = dn[name].i = ++nb_grids;
   dn[name].dt = DataType::GRID;
   return iGrid;
}


int data::addMesh(OFELI::Mesh* ms,
                  string&      name)
{
   iMesh = theMesh.size();
   DEFAULT_MESH_NAME(name);
   if (checkAlreadyUsed(name,DataType::MESH))
      return 0;
   Desc[name] = "";
   if (MeshLabel[name]) {
      REDEFINED("Mesh ")
      iMesh = MeshLabel[name];
      dn[name].active = true;
      theMesh[iMesh] = ms;
      return iMesh;
   }
   theMesh.push_back(ms);
   MeshName.push_back(name);
   dn[name].active = true;
   MeshLabel[name] = dn[name].i = ++nb_meshes;
   dn[name].dt = DataType::MESH;
   return iMesh;
}


int data::addParam(string name,
                   double value,
                   bool   opt)
{
   if (checkAlreadyUsed(name,DataType::PARAM))
      return 0;
   Desc[name] = "";
   if (!opt)
      _rita->_calc->setParam(name,value);
   if (ParamLabel[name]) {
      iParam = ParamLabel[name];
      dn[name].active = true;
      theParam[iParam] = value;
      ParamLabel[name] = dn[name].i = iParam;
      dn[name].dt = DataType::PARAM;
      return iParam;
   }
   iParam = theParam.size();
   theParam.push_back(value);
   ParamName.push_back(name);
   dn[name].active = true;
   ParamLabel[name] = dn[name].i = ++nb_params;
   dn[name].dt = DataType::PARAM;
   return iParam;
}


int data::addVector(string name,
                    double t,
                    int    n,
                    string file,
                    int    opt)
{
   int nbd = 1;
   iVector = theVector.size();
   DEFAULT_VECTOR_NAME(name);
   if (checkAlreadyUsed(name,DataType::VECTOR))
      return 0;
   Desc[name] = "";
   if (file!="")
      n = 1;
   _u = new OFELI::Vect<double>(n);
   if (file!="") {
      OFELI::XMLParser xml(file,EType::VECTOR);
      xml.get(*_u);
      nbd = _u->getNbDOF();
   }
   _u->setName(name);
   vect_hist[name] = "%$§&";
   dn[name].dt = DataType::VECTOR;
   if (VectorLabel[name]) {
      if (!opt)
         _rita->_calc->setVectorValue(name,_u);
      iVector = VectorLabel[name];
      if (dn[name].active==false)
         nb_vectors++;
      if (file!="") {
         OFELI::XMLParser xml(file,EType::VECTOR);
         xml.get(*_u);
      }
      *theVector[iVector] = *_u;
      VectorTime[iVector] = t;
      VectorName[iVector] = name;
      nb_dof[iVector] = nbd;
      dn[name].active = true;
      return iVector;
   }
   if (!opt)
      _rita->_calc->setVector(_u);
   VectorName.push_back(name);
   theVector.push_back(_u);
   VectorTime.push_back(t);
   dn[name].active = true;
   VectorSizeType.push_back(DataSize::GIVEN_SIZE);
   VectorLabel[name] = dn[name].i = ++nb_vectors;
   nb_dof.push_back(nbd);
   return iVector;
}


void data::setTime(string s,
                   double t)
{
   int it = 0;
   if ((it=VectorLabel[s]))
      VectorTime[it] = t;
}


DataType data::getType(string s)
{
   DataType dt = dn[s].dt;
   if (dt<DataType::PARAM || dt>DataType::EIGEN)
      return DataType::NOTHING;
   else
      return dt;
}


int data::remove(const string& name)
{
   DataType d = getType(name);
   switch (d) {

      case DataType::PARAM:
         if (ParamLabel[name]) {
            dn[name].active = false;
            _rita->_calc->theParser.RemoveVar(name);
            nb_params--;
         }
         break;

      case DataType::VECTOR:
         if (VectorLabel[name]) {
            dn[name].active = false;
            _rita->_calc->theParser.RemoveVar(name);
            nb_vectors--;
         }
         break;

      case DataType::HVECTOR:
         if (HVectorLabel[name]) {
            dn[name].active = false;
            nb_hvectors--;
         }
         break;

      case DataType::MATRIX:
         if (MatrixLabel[name]) {
            theMatrix[MatrixLabel[name]]->setSize(0);
            dn[name].active = false;
            _rita->_calc->theParser.RemoveVar(name);
            nb_matrices--;
         }
         break;

      case DataType::GRID:
         if (GridLabel[name]) {
            delete theGrid[GridLabel[name]];
            theGrid[GridLabel[name]] = nullptr;
            dn[name].active = false;
            nb_grids--;
         }
         break;

      case DataType::MESH:
         if (MeshLabel[name]) {
            delete theMesh[MeshLabel[name]];
            theMesh[MeshLabel[name]] = nullptr;
            dn[name].active = false;
            nb_meshes--;
         }
         break;

      case DataType::TAB:
         if (TabLabel[name]) {
            delete theTab[TabLabel[name]];
            theTab[TabLabel[name]] = nullptr;
            dn[name].active = false;
            nb_tabs--;
         }
         break;

      case DataType::FCT:
         if (FctLabel[name]) {
            delete theFct[FctLabel[name]];
            theFct[FctLabel[name]] = nullptr;
            dn[name].active = false;
            nb_fcts--;
         }
         break;

      case DataType::LS:
         if (LSLabel[name]) {
            delete theLS[LSLabel[name]];
            theLS[LSLabel[name]] = nullptr;
            dn[name].active = false;
            nb_ls--;
         }
         break;

      case DataType::AE:
         if (AELabel[name]) {
            delete theAE[AELabel[name]];
            theAE[AELabel[name]] = nullptr;
            dn[name].active = false;
            nb_ae--;
         }
         break;

      case DataType::ODE:
         if (ODELabel[name]) {
            delete theODE[ODELabel[name]];
            theODE[ODELabel[name]] = nullptr;
            dn[name].active = false;
            nb_ode--;
         }
         break;

      case DataType::PDE:
         if (PDELabel[name]) {
            delete thePDE[PDELabel[name]];
            thePDE[PDELabel[name]] = nullptr;
            dn[name].active = false;
            nb_pde--;
         }
         break;

      case DataType::OPTIM:
         if (OptLabel[name]) {
            delete theOpt[OptLabel[name]];
            theOpt[OptLabel[name]] = nullptr;
            dn[name].active = false;
            nb_opt--;
         }
         break;

      case DataType::EIGEN:
         if (EigLabel[name]) {
            delete theEig[EigLabel[name]];
            theEig[EigLabel[name]] = nullptr;
            dn[name].active = false;
            nb_eig--;
         }
         break;

      default:
         MSGR("","No data named "+name+" found")
   }
   return 0;
}


int data::addMatrix(string name,
                    int    nr,
                    int    nc,
                    string file,
                    string s,
                    bool   opt)
{
   iMatrix = theMatrix.size();
   string nm = name;
   if (checkAlreadyUsed(name,DataType::MATRIX))
      return 0;
   Desc[name] = "";
   if (file=="") {
      if (s=="dense")
         _theMatrix = new DMatrix<double>(nr,nc);
      else if (s=="dense-symmetric")
         _theMatrix = new DSMatrix<double>(nr);
      else if (s=="sparse")
         _theMatrix = new SpMatrix<double>(nr,nc);
      else if (s=="tridiagonal")
         _theMatrix = new TrMatrix<double>(nr);
   }
   else {
      ifstream f(file.c_str());
      CHK_MSGR(!f.good(),"","File "+file+" not found")
      OFELI::XMLParser xml(file,EType::MATRIX);
      xml.get(_theMatrix);
      MatrixSize msize = xml.MSize();
      if (nm=="")
         name = msize.name;
      if (msize.mt==DENSE)
         _theMatrix = new DMatrix<double>(msize.nb_rows,msize.nb_cols);
      else if (msize.mt==TRIDIAGONAL)
         _theMatrix = new TrMatrix<double>(msize.size);
      else if (msize.mt==BAND)
         _theMatrix = new BMatrix<double>(msize.size,msize.ld,msize.ud);
      else if (msize.mt==SKYLINE)
         _theMatrix = new SkMatrix<double>(msize.ch);
      else if (msize.mt==SPARSE)
         _theMatrix = new SpMatrix<double>(msize.IJ);
      _theMatrix->setName(name);
      xml.get(_theMatrix);
   }
   DEFAULT_MATRIX_NAME(name);
   dn[name].dt = DataType::MATRIX;
   if (MatrixLabel[name]) {
      iMatrix = MatrixLabel[name];
      if (dn[name].active==false)
         nb_matrices++;
      dn[name].active = true;
      theMatrix[iMatrix] = _theMatrix;
      MatrixLabel[name] = dn[name].i = nb_matrices;
      if (!opt)
         _rita->_calc->setMatrixValue(name,_theMatrix);
      return iMatrix;
   }
   MatrixName.push_back(name);
   theMatrix.push_back(_theMatrix);
   dn[name].active = true;
   MatrixLabel[name] = dn[name].i = ++nb_matrices;
   if (!opt)
      _rita->_calc->setMatrix(_theMatrix);
   return iMatrix;
}


int data::addTab(OFELI::Tabulation* tab,
                 string&            name)
{
   iTab = theTab.size();
   DEFAULT_TAB_NAME(name);
   if (checkAlreadyUsed(name,DataType::TAB))
      return 0;
   Desc[name] = "";
   if (TabLabel[name]) {
      REDEFINED("Tabulation ")
      iTab = TabLabel[name];
      dn[name].active = true;
      theTab[iTab] = tab;
      return iTab;
   }
   theTab.push_back(tab);
   TabName.push_back(name);
   dn[name].active = true;
   TabLabel[name] = dn[name].i = ++nb_tabs;
   dn[name].dt = DataType::TAB;
   return iTab;
}


int data::addMeshVector(const string& nm1,
                        string&       nm2,
                        DataSize      s,
                        int           ndof,
                        bool          opt)
{
   iVector = theVector.size();
   DEFAULT_VECTOR_NAME(nm2);
   if (checkAlreadyUsed(nm2,DataType::VECTOR))
      return 0;
   CHK_MSGR0(MeshLabel[nm1]==0,"","Mesh "+nm1+" not found")
   _theMesh = theMesh[MeshLabel[nm1]];
   CHK_MSGR(s==DataSize::NODES&&_theMesh->getNbNodes()==0,"","Mesh has no nodes")
   CHK_MSGR(s==DataSize::ELEMENTS&&_theMesh->getNbElements()==0,"","Mesh has no elements.")
   CHK_MSGR(s==DataSize::SIDES&&_theMesh->getNbSides()==0,"","Mesh has no sides.")
   CHK_MSGR(s==DataSize::EDGES&&_theMesh->getNbEdges()==0,"","Mesh has no edges.")
   _u = new OFELI::Vect<double>;
   _u->setName(nm2);
   if (s==DataSize::NODES) {
      _u->setMesh(*_theMesh,NODE_DOF,ndof);
      VectorSizeType.push_back(DataSize::NODES);
   }
   else if (s==DataSize::ELEMENTS) {
      _u->setMesh(*_theMesh,ELEMENT_DOF,ndof);
      VectorSizeType.push_back(DataSize::ELEMENTS);
   }
   else if (s==DataSize::SIDES) {
      _u->setMesh(*_theMesh,SIDE_DOF,ndof);
      VectorSizeType.push_back(DataSize::SIDES);
   }
   else if (s==DataSize::EDGES) {
      _u->setMesh(*_theMesh,EDGE_DOF,ndof);
      VectorSizeType.push_back(DataSize::EDGES);
   }
   Desc[nm2] = "Vector associated to mesh " + nm1;
   dn[nm2].dt = DataType::VECTOR;
   if (VectorLabel[nm2]) {
      iVector = VectorLabel[nm2];
      dn[nm2].active = true;
      theVector[iVector] = _u;
      VectorLabel[nm2] = dn[nm2].i = ++nb_vectors;
      if (!opt)
         _rita->_calc->setVectorValue(nm2,_u);
      return iVector;
   }
   VectorName.push_back(nm2);
   theVector.push_back(_u);
   nb_dof.push_back(ndof);
   dn[nm2].active = true;
   VectorLabel[nm2] = dn[nm2].i = ++nb_vectors;
   if (!opt)
      _rita->_calc->setVector(_u);
   return iVector;
}


int data::addGridVector(const string& nm1,
                        string&       nm2,
                        int           ndof,
                        bool          opt)
{
   iVector = theVector.size();
   DEFAULT_VECTOR_NAME(nm2);
   if (checkAlreadyUsed(nm2,DataType::VECTOR))
      return 0;
   CHK_MSGR0(GridLabel[nm1]==0,"","Grid "+nm1+" not found")
   _theGrid = theGrid[GridLabel[nm1]];
   _u = new OFELI::Vect<double>(*_theGrid);
   _u->setName(nm2);
   _theGrid->setNbDOF(ndof);
   dn[nm2].dt = DataType::VECTOR;
   Desc[nm2] = "Vector associated to grid " + nm1;
   if (VectorLabel[nm2]) {
      iVector = VectorLabel[nm2];
      dn[nm2].active = true;
      theVector[iVector] = _u;
      _rita->_calc->setVectorValue(nm2,_u);
      return iVector;
   }
   VectorName.push_back(nm2);
   theVector.push_back(_u);
   dn[nm2].active = true;
   nb_dof.push_back(ndof);
   VectorLabel[nm2] = dn[nm2].i = ++nb_vectors;
   VectorSizeType.push_back(DataSize::GRID);
   if (!opt)
      _rita->_calc->setVector(_u);
   return iVector;
}


int data::setVectorValue(int                        k,
                         const OFELI::Vect<double>* v)
{
   CHK_MSGR(k<=0,"","Vector not found.")
   *theVector[k] = *v;
   _rita->_calc->setVectorValue(VectorName[k],v);
   return 0;
}


int data::addVector(OFELI::Vect<double>* v,
                    string&              name,
                    bool                 opt)
{
   iVector = theVector.size();
   DEFAULT_VECTOR_NAME(name);
   v->setName(name);
   Desc[name] = "";
   if (!opt)
      _rita->_calc->setVector(v);
   if (VectorLabel[name]) {
      REDEFINED("Vector ")
      iVector = VectorLabel[name];
      dn[name].active = true;
      theVector[iVector] = _u;
      return iVector;
   }
   iVector = theVector.size();
   theVector.push_back(v);
   VectorName.push_back(name);
   dn[name].active = true;
   VectorLabel[name] = dn[name].i = ++nb_vectors;
   VectorSizeType.push_back(DataSize::GIVEN_SIZE);
   nb_dof.push_back(1);
   dn[name].dt = DataType::VECTOR;
   return iVector;
}


int data::print(const string& s)
{
   if (dn[s].i==0) {
      cout << "Entity " << s << " not found." << endl;
      return 1;
   }
   if (!dn[s].active) {
      cout << "Entity " << s << " has been removed." << endl;
      return 1;
   }
   DataType d = getType(s);
   switch (d) {

      case DataType::PARAM:
         cout << s << " = " << theParam[ParamLabel[s]] << endl;
         break;

      case DataType::VECTOR:
         cout << s << " =\n" << *theVector[VectorLabel[s]] << endl;
         break;

      case DataType::HVECTOR:
         cout << s << " = " << *theHVector[HVectorLabel[s]] << endl;
         break;

      case DataType::MATRIX:
         {
            cout << "Matrix " << s << endl;
            int k = MatrixLabel[s];
            for (size_t i=1; i<=theMatrix[k]->getNbRows(); ++i) {
               cout << "Row " << i << ": ";
               for (size_t j=1; j<=theMatrix[k]->getNbColumns(); ++j)
                  cout << theMatrix[k]->at(i,j) << "  ";
               cout << endl;
            }
         }
         break;

      case DataType::MESH:
         cout << *theMesh[MeshLabel[s]];
         break;

      case DataType::GRID:
         cout << *theGrid[GridLabel[s]];
         break;

      case DataType::FCT:
         cout << *theFct[FctLabel[s]];
         break;

      case DataType::TAB:
         cout << *theTab[TabLabel[s]];
         break;

      case DataType::LS:
         cout << *theLS[LSLabel[s]];
         break;

      case DataType::AE:
         cout << *theAE[AELabel[s]];
         break;

      case DataType::ODE:
         cout << *theODE[ODELabel[s]];
         break;

      case DataType::PDE:
         cout << *thePDE[PDELabel[s]];
         break;

      case DataType::OPTIM:
         cout << *theOpt[OptLabel[s]];
         break;

      case DataType::EIGEN:
         cout << *theEig[EigLabel[s]];
         break;

      default:
         MSGR("","No entity named "+s+" found")
   }
   return 0;
}


int data::setConfigure()
{
   _configure->setVerbose(_verb);
   int ret = _configure->run();
   _verb = _configure->getVerbose();
   return ret;
}


int data::setGrid()
{
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, dim=0, dd=0, nn=0, nx=10, ny=10, nz=10;
   string name="";
   static const vector<string> kw {"name","x","y","z","ne"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,"grid>","No argument given for command.")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Grid_help << endl;
            return 0;

         case 101:
            cout << Grid_Help << endl;
            return 0;

         case   0:
            name = _cmd->string_token(0);
            break;

         case   1:
            CHK_MSGR(nb!=2,"grid>","Illegal number of arguments")
            dim = 1, nn++;
            DDEF_PAR_R(0,"grid>",xmin)
            DDEF_PAR_R(1,"grid>",xmax)
            break;

         case   2:
            CHK_MSGR(nb!=2,"grid>","Illegal number of arguments")
            dim = 2;
            DDEF_PAR_R(0,"grid>",ymin)
            DDEF_PAR_R(1,"grid>",ymax)
            break;

         case   3:
            CHK_MSGR(nb!=2,"grid>","Illegal number of arguments")
            dim = 3;
            DDEF_PAR_R(0,"grid>",zmin)
            DDEF_PAR_R(1,"grid>",zmax)
            break;

         case   4:
            CHK_MSGR(nb==0 || nb>3,"grid>","Illegal number of arguments")
            dd = nb;
            DDEF_PAR_R(0,"grid>",nx)
            if (nb>1)
               DDEF_PAR_R(1,"grid>",ny)
            if (nb>2)
               DDEF_PAR_R(2,"grid>",nz)
            break;

         default:
            UNKNOWN_ARG("grid>")
      }
   }
   CHK_MSGR(dd!=dim,"grid>","Dimensions do not match.")
   CHK_MSGR(xmin>=xmax || ymin>=ymax || zmin>=zmax,"grid>","Domain definition is incorrect.")
   CHK_MSGR(nx<=0,"grid>","Number of x-subdivisions is incorrect.")
   CHK_MSGR(dim==2 && ny<=0,"grid>","Number of y-subdivisions is incorrect.")
   CHK_MSGR(dim==3 && nz<=0,"grid>","Number of z-subdivisions is incorrect.")
   *_rita->ofh << "grid";
   if (name!="")
      *_rita->ofh << " name=" << name;
   if (dim==1) {
      _theGrid = new OFELI::Grid(xmin,xmax,nx);
      *_rita->ofh << " x=" << xmin << "," << xmax << " ne=" << nx;
   }
   else if (dim==2) {
      _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,nx,ny);
      *_rita->ofh << " x=" << xmin << "," << xmax << " y=" << ymin << "," << ymax << " ne=" << nx << "," << ny;
   }
   else if (dim==3) {
      _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
      *_rita->ofh << " x=" << xmin << "," << xmax << " y=" << ymin << "," << ymax 
                  << " z=" << zmin << "," << zmax << " ne=" << nx << "," << ny << "," << nz;
   }
   *_rita->ofh << endl;
   addGrid(_theGrid,name);
   return 0;
}


int data::setVector()
{
   OFELI::Vect<double> *v;
   int nb=0, size=0, ret=0;
   string name="V"+to_string(nbAllVector()+1);
   string name_grid="", name_mesh="", file="";
   static const vector<string> kw {"name","size","grid","mesh","file","def$ine","set"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,"vector>","No argument given for command.")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Vect_help << endl;
            return 0;

         case 101:
            cout << Vect_Help << endl;
            return 0;

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
            DDEF_PAR_R(0,"vector>",size)
            break;

         case 2:
            name_grid = _cmd->string_token(0);
            break;

         case 3:
            name_mesh = _cmd->string_token(0);
            break;

         case 4:
            file = _cmd->string_token(0);
            break;

         case 5:
            break;

         case 6:
            break;

         default:
            UNKNOWN_ARG("vector>")
      }
   }
   if (nb_args>0) {
      CHK_MSGR(name_grid!="" && size>0,"vector>","Grid data and vector size cannot be given simultaneously")
      CHK_MSGR(name_mesh!="" && size>0,"vector>","Mesh data and vector size cannot be given simultaneously")
      CHK_MSGR(name_grid!="" && name_mesh!="","vector>","Grid and Mesh data cannot be given simultaneously")
      CHK_MSGR(file!="" && size>0,"vector>","file and vector size cannot be given simultaneously")
      CHK_MSGR(name_grid=="" && name_mesh=="" && file=="" && size==0,"vector>","No size given for vector")
      *_rita->ofh << "vector";
      if (name_grid!="") {
         int k = GridLabel[name];
         CHK_MSGR(k==0,"vector>","Reference to undefined grid.")
         _theGrid = theGrid[k];
         v = new OFELI::Vect<double>(*_theGrid);
         addVector(v,name);
         *_rita->ofh << " grid=" << name_grid;
      }
      if (name_mesh!="") {
         int k = MeshLabel[name];
         CHK_MSGR(k<=0,"vector>","Reference to undefined mesh.")
         _theMesh = theMesh[k];
         v = new OFELI::Vect<double>(*_theMesh);
         addVector(v,name);
         *_rita->ofh << " mesh=" << name_mesh;
      }
      if (size>0) {
         addVector(name,0.,size,"",false);
         *_rita->ofh << " size=" << size;
      }
      *_rita->ofh << endl;
   }
   return ret;
}


int data::setMatrix()
{
   string name="", file="", storage="dense";
   int nr=0, nc=0, nb=0;
   static const vector<string> kw {"name","file","storage","nr","nc","def$ine","set"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,"matrix>","No agument given for command.")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Matrix_help << endl;
            return 0;

         case 101:
            cout << Matrix_Help << endl;
            return 0;

         case   0:
            name = _cmd->string_token(0);
            break;

         case   1:
            file = _cmd->string_token(0);
            break;

         case   2:
            storage = _cmd->string_token(0);
            break;

         case   3:
            DDEF_PAR_R(0,"matrix>",nr)
            if (nc==0) nc = nr;
            break;

         case   4:
            DDEF_PAR_R(0,"matrix>",nc)
            if (nr==0) nr = nc;
            break;

         default:
            UNKNOWN_ARG("matrix>")
      }
   }
   CHK_MSGR(file=="" && nr==0,"matrix>","Matrix must be given either by its size or read from a file.")
   CHK_MSGR(storage!="dense","matrix>","Only dense storage is allowed in the current release.")
   *_rita->ofh << "matrix";
   if (file!="")
      *_rita->ofh << " file=" << file;
   if (nr)
      *_rita->ofh << " size=" << nr;
   *_rita->ofh << endl;
   addMatrix(name,nr,nc,file,storage,false);
   return 0;
}


int data::setFunction()
{
   _ret = 0;
   _pr = "function>";
   bool var_ok=false, def_ok=false;
   string vv="", def="", name="";
   int nb=0, nbv=1;
   vector<string> var;
   static const vector<string> kw {"name","var$iable","vec$tor","nb","def$inition"};

   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,"function>","No argument for command")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Fct_help << endl;
            return 0;

         case 101:
            cout << Fct_Help << endl;
            return 0;

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
         case 2:
            for (int i=0; i<nb; ++i) {
               vv = _cmd->string_token(i);
               var.push_back(vv);
            }
            var_ok = true;
            break;

         case 3:
            DDEF_PAR_R(0,"",nbv)
            break;

         case 4:
            def = _cmd->string_token(0);
            def_ok = true;
            break;

         default:
            UNKNOWN_ARG("function>")
      }
   }
   CHK_MSGR(nb_args==0,"function>","No command argument given.")
   CHK_MSGR(!var_ok,"function>","No variable defined.")
   CHK_MSGR(!def_ok,"function>","No function definition given.")
   int k = addFunction(name,var,def);
   FCT_NOT_DEFINED("",name)
   FCT_ALREADY_DEFINED("",name)
   if (k>0) {
      *_rita->ofh << "function name=" << name;
      if (nbv>1)
         *_rita->ofh << " nb=" << nbv;
      *_rita->ofh << " variable=" << var[0];
      *_rita->ofh << " definition=" << def << endl;
   }
   return 0;
}


int data::setNbDOF()
{
   CHK_MSGR(_cmd->setNbArg(1,"Give number of degrees of freedom."),"nbdof>","Missing value of nbdof.")
   int ndof = 1;
   _ret = _cmd->get(ndof);
   if (!_ret)
      *_rita->ofh << "  nbdof " << ndof;
   return _ret;
}


int data::setTab()
{
   int dim1=0, dim2=0, dim3=0;
   string file="", name="Tab-"+to_string(nbAllTab()+1);
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, nx=10, ny=10, nz=10, grid_ok=0, file_ok=0, vector_ok=0;
   static const vector<string> kw {"name","file","min","max","ne","vector"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args==0,"tabulation>","No argument for command")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " + Tab_help << endl;
            return 0;

         case 101:
            cout << Tab_Help << endl;
            return 0;

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
            file = _cmd->string_token(0);
            file_ok = 1;
            break;

         case 2:
            dim1 = nb;
            DDEF_PAR_R(0,"tabulation>",xmin)
            if (nb>1)
               DDEF_PAR_R(1,"tabulation>",ymin)
            if (nb>2)
               DDEF_PAR_R(2,"tabulation>",zmin)
            grid_ok += 1;
            break;

         case 3:
            dim2 = nb;
            DDEF_PAR_R(0,"tabulation>",xmax)
            if (nb>1)
               DDEF_PAR_R(1,"tabulation>",ymax)
            if (nb>2)
               DDEF_PAR_R(2,"tabulation>",zmax)
            grid_ok += 10;
            break;

         case 4:
            dim3 = nb;
            DDEF_PAR_R(0,"tabulation>",nx)
            if (nb>1)
               DDEF_PAR_R(1,"tabulation>",ny)
            if (nb>2)
               DDEF_PAR_R(2,"tabulation>",nz)
            grid_ok += 100;
            break;

         case 5:
	   //            fd = _cmd->string_token();
            vector_ok = true;
            break;

         default:
            UNKNOWN_ARG("tabulation>")
      }
   }
   CHK_MSGR(!file_ok && grid_ok<111,"tabulation>","No grid data given.")
   CHK_MSGR(!vector_ok && !file_ok,"tabulation>","No associated vector given.")
   CHK_MSGR(!file_ok && (dim1!=dim2 || dim1!=dim3 || dim2!=dim3),"tabulation>",
            "Incompatible space dimensions as given by grid data.")
   *_rita->ofh << "tabulation name=" << name;
   if (file_ok) {
      ifstream ip(file);
      CHK_MSGR(!ip.is_open(),"tabulation>","Unable to open file: "+file)
      ip.close();
      _theTab = new OFELI::Tabulation(file);
      *_rita->ofh << " file=" << file;
   }
   else {
      if (dim1==1)
         _theGrid = new OFELI::Grid(xmin,xmax,nx);
      else if (dim1==2)
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,nx,ny);
      if (dim1==3)
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
   }
   addTab(_theTab,name);
   *_rita->ofh << endl;
   return 0;
}


int data::getPar(int n, const string& msg, int& v)
{
   string s="";
   if (n<0)
      _cmd->get(s);
   else
      s = _cmd->string_token(n);
   if (_cmd->isNumeric(s,v))
      return 0;
   CHK_MSGR(ParamLabel[s]==0,msg,"Undefined parameter "+s)
   v = _rita->_calc->get_int(s);
   return 0;
}


int data::getPar(int n, const string& msg, double& v)
{
   string s="";
   if (n<0) {
      if (_cmd->get(s))
         return 1;
   }
   else
      s = _cmd->string_token(n);
   if (_cmd->isNumeric(s,v))
      return 0;
   CHK_MSGR(ParamLabel[s]==0,msg,"Undefined parameter "+s)
   v = _rita->_calc->get_double(s);
   return 0;
}


int data::setDesc(const string& name,
                  const string& desc)
{
   if (getType(name)==DataType::NOTHING)
      return 1;
   Desc[name] = desc;
   return 0;
}


int data::save()
{
   string name="", file="", format="ofeli", t="";
   double every=1, e=1, ret=0;
   bool name_ok=false, file_ok=false;
   static const vector<string> kw {"name","file","format","every"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   int nb=0;
   CHK_MSGR(nb_args==0,"save>","No argument given.")
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Save_help << endl;
            return 0;

         case 101:
            cout << Save_Help << endl;
            return 0;

         case 0:
            if ((t=_cmd->string_token(0))!="")
               name = t;
            name_ok = true;
            break;

         case 1:
            if ((t=_cmd->string_token(0))!="")
               file = t;
            file_ok = true;
            break;

         case 2:
            if ((t=_cmd->string_token(0))!="")
               format = t;
            break;

         case 3:
            if ((e=_cmd->int_token(0))!=0)
               every = e;
            break;

         default:
            UNKNOWN_ARG("save>")
      }
   }
   CHK_MSGR(name_ok==false,"save>","No data to save.")
   CHK_MSGR(file_ok==false,"save>","No file given.")
   CHK_MSGR(dn[name].i==0,"save>","Data "+name+" not found. Nothing to save!")
   *_rita->ofh << "save name=" << name;

// Case of a parameter
   CHK_MSGR(dn[name].dt==DataType::PARAM,"save>","Parameters cannot be saved.")

// Save a vector
   if (dn[name].dt==DataType::VECTOR && dn[name].active) {
      int k=VectorLabel[name];
      if (format=="ofeli") {
         *_rita->ofh << " format=ofeli file=" << file;
         OFELI::IOField ffo(file,OFELI::IOField::OUT);
         ffo.put(*theVector[k]);
         return 0;
      }
      else if (format=="gmsh") {
         *_rita->ofh << " format=gmsh file=" << file;
         int ret = saveGmsh(file,*theVector[k]);
         CHK_MSGR(ret==1,"save>","Cannot save in gmsh format a vector without mesh association.")
      }
      else if (format=="vtk") { 
         *_rita->ofh << " format=vtk file=" << file;
         int ret = saveVTK(file,*theVector[k]);
         CHK_MSGR(ret==1,"save>","Cannot save in gmsh format a vector without mesh association.")
      }
      else if (format=="tecplot") { 
         *_rita->ofh << " format=tecplot file=" << file;
         int ret = saveTecplot(file,*theVector[k]);
         CHK_MSGR(ret==1,"save>","Cannot save in gmsh format a vector without mesh association.")
      }
      else if (format=="gnuplot") {
         _theMesh = &(theVector[k]->getMesh());
         _theGrid = &(theVector[k]->getGrid());
         CHK_MSGR(_theMesh==nullptr && _theGrid==nullptr,"save>",
                  "Vector "+VectorName[k]+" cannot be plotted in gnuplot format: No associated mesh or grid detected.")
         *_rita->ofh << " format=gnuplot file=" << file;
         saveGnuplot(file,*theVector[k]);
         return 0;
      }
      else
         CHK_MSGR(ret==1,"save>","Format "+format+" is not available for saving vector.")
   }

// Save a history vector
   else if (dn[name].dt==DataType::HVECTOR && dn[name].active) {
      CHK_MSGR(theHVector[HVectorLabel[name]]->nt==0,"save>","History vector "+name+" empty.")
      if (format=="ofeli") {
         int ret = theHVector[HVectorLabel[name]]->saveOFELI(file,every);
         *_rita->ofh << " format=ofeli file=" << file << " every=" << every;
         return ret;
      }
      else if (format=="gmsh") {
         int ret = theHVector[HVectorLabel[name]]->saveGmsh(file,every);
         CHK_MSGR(ret==1,"save>","Cannot save in gmsh format a vector without mesh association.")
         *_rita->ofh << " format=gmsh file=" << file << " every=" << every;
      }
      else if (format=="vtk") {
         int ret = theHVector[HVectorLabel[name]]->saveVTK(file,every);
         CHK_MSGR(ret==1,"save>","Cannot save in vtk format a vector without mesh association.")
         *_rita->ofh << " format=vtk file=" << file << " every=" << every;
      }
      else if (format=="tecplot") {
         int ret = theHVector[HVectorLabel[name]]->saveTecplot(file,every);
         CHK_MSGR(ret==1,"save>","Cannot save in vtk format a vector without mesh association.")
         *_rita->ofh << " format=vtk file=" << file << " every=" << every;
      }
      else if (format=="gnuplot") {
         int ret = theHVector[HVectorLabel[name]]->saveGnuplot(file,every);
         *_rita->ofh << " format=gnuplot file=" << file << " every=" << every;
         return ret;
      }
      else
         CHK_MSGR(ret==1,"save>","Format "+format+" is not available for saving history vector.")
   }

// Save a matrix
   else if (dn[name].dt==DataType::MATRIX && dn[name].active) {
      OFELI::saveMatrix(theMatrix[MatrixLabel[name]],file);
      *_rita->ofh << " file=" << file;
      return 0;
   }

// Save a mesh
   else if (dn[name].dt==DataType::MESH && dn[name].active) {
      if (format=="ofeli") {
         saveMesh(file,*theMesh[MeshLabel[name]],OFELI_FF);
         *_rita->ofh << " file=" << file << " format=ofeli";
         return 0;
      }
      else if (format=="gmsh") {
         saveMesh(file,*theMesh[MeshLabel[name]],GMSH);
         *_rita->ofh << " file=" << file << " format=gmsh";
         return 0;
      }
      else if (format=="vtk") {
         saveMesh(file,*theMesh[MeshLabel[name]],VTK);
         *_rita->ofh << " file=" << file << " format=vtk";
         return 0;
      }
      else if (format=="gnuplot") {
         saveMesh(file,*theMesh[MeshLabel[name]],GNUPLOT);
         *_rita->ofh << " file=" << file << " format=vtk";
         return 0;
      }
      else
         CHK_MSGR(ret==1,"save>","Format "+format+" is not available for saving mesh.")
   }

// Save a grid
   CHK_MSGR(dn[name].dt==DataType::GRID && dn[name].active,"save>","This option is not yet implemented.")

// Case of a function
   CHK_MSGR(dn[name].dt==DataType::FCT && dn[name].active,"save>","Functions cannot be saved.")

// Save a tabulation
   CHK_MSGR(dn[name].dt==DataType::TAB && dn[name].active,"save>","This type of data cannot be saved in file.")
   *_rita->ofh << endl;
   return ret;
}


void data::ListParams(int opt)
{
   if (opt && !nb_params) {
      cout << "No defined parameters." << endl;
      return;
   }
   cout << "Number of parameters: " << nb_params << endl;
   for (size_t i=1; i<theParam.size(); ++i) {
      string s = ParamName[i];
      if (dn[s].active) {
         cout << s << " = " << theParam[i] << endl;
         if (Desc[s]!="")
            cout << "Description: " << Desc[ParamName[i]] << endl;
      }
   }
}


void data::ListVectors(int opt)
{
   if (opt && !nb_vectors) {
      cout << "No defined vectors." << endl;
      return;
   }
   cout << "Number of vectors: " << nb_vectors << endl;
   for (size_t i=1; i<theVector.size(); ++i) {
      string s = VectorName[i];
      if (dn[s].active) {
         if (VectorSizeType[i]==DataSize::GIVEN_SIZE)
            cout << "Vector: " << s << ", Size: " << theVector[i]->size() << endl;
         else
            cout << "Vector: " << s << ", Number of degrees of freedom: " << nb_dof[i] << endl;
      }
   }
}


void data::ListHVectors(int opt)
{
   if (opt && !nb_hvectors) {
      cout << "No defined history vectors." << endl;
      return;
   }
   cout << "Number of history vectors: " << nb_hvectors << endl;
   for (size_t i=1; i<theHVector.size(); ++i) {
      string s = HVectorName[i];
      if (dn[s].active)
         cout << "History Vector: " << s << ", Size: " << theHVector[i]->nt << endl;
   }
}


void data::ListFunctions(int opt)
{
   if (opt && !nb_fcts) {
      cout << "No defined functions." << endl;
      return;
   }
   cout << "Number of functions: " << nb_fcts << endl;
   for (size_t i=1; i<theFct.size(); ++i) {
      string s = FctName[i];
      if (dn[s].active)
         cout << *theFct[i];
   }
}


void data::ListTabs(int opt)
{
   if (opt && !nb_tabs) {
      cout << "No defined tabulations." << endl;
      return;
   }
   cout << "Number of tabulations: " << nb_tabs << endl;
   for (size_t i=1; i<theTab.size(); ++i) {
      string s = TabName[i];
      if (dn[s].active)
         cout << "Tabulation: " << s << ", Nb. of variables: " 
              << theTab[i]->getNbVar(1) << ", Size: " << theTab[i]->getSize(1,1) << endl;
   }
}


void data::ListMatrices(int opt)
{
   if (opt && !nb_matrices) {
      cout << "No defined matrices." << endl;
      return;
   }
   cout << "Number of matrices: " << nb_matrices << endl;
   for (size_t i=1; i<theMatrix.size(); ++i) {
      string s = MatrixName[i];
      if (dn[s].active)
         cout << "Matrix: " << s << ", Size: " << theMatrix[i]->getNbRows()
              << " x " << theMatrix[i]->getNbColumns() << endl;
   }
}


void data::ListGrids(int opt)
{
   if (opt && !nb_grids) {
      cout << "No defined grids." << endl;
      return;
   }
   cout << "Number of grids: " << nb_grids << endl;
   for (size_t i=1; i<theGrid.size(); ++i) {
      string s = GridName[i];
      if (dn[s].active) {
         OFELI::Grid *g = theGrid[i];
         cout << "Grid name:           " << s << endl;
         cout << "Space dimension:     " << g->getDim() << endl;
         if (g->getDim()==1) {
            cout << "Domain:              (" << g->getX(1) << ","
                 << g->getX(theGrid[i]->getNx()+1) << ")" << endl;
            cout << "Number of intervals: " << g->getNx() << endl;
         }
         else if (g->getDim()==2) {
            cout << "Domain:                    (" << g->getX(1) << ","
                 << g->getX(g->getNx()+1) << ")x(" << g->getY(1) << ","
                 << g->getY(g->getNy()+1) << ")" << endl;
            cout << "Number of intervals:    " << g->getNx() << " x " << g->getNy() << endl;
         }
         else if (g->getDim()==3) {
            cout << "Domain:                    (" << g->getX(1) << ","
                 << g->getX(g->getNx()+1) << ")x(" << g->getY(1) << ","
                 << g->getY(g->getNy()+1) << ")x(" << g->getZ(1) << ","
                 << g->getZ(g->getNz()+1) << ")" << endl;
            cout << "Number of intervals:    " << g->getNx() << " x " << g->getNy() << " x " << g->getNz() << endl;
         }
      }
   }
}


void data::ListMeshes(int opt)
{
   if (opt && !nb_meshes) {
      cout << "No defined meshes." << endl;
      return;
   }
   cout << "Number of meshes: " << nb_meshes << endl;
   for (size_t i=1; i<theMesh.size(); ++i) {
      string s = MeshName[i];
      if (dn[s].active) {
         OFELI::Mesh *m = theMesh[i];
         cout << "Mesh name:          " << s << endl;
         cout << "Number of nodes:    " << m->getNbNodes() << endl;
         cout << "Number of elements: " << m->getNbElements() << endl;
         cout << "Number of sides:    " << m->getNbSides() << endl;
      }
   }
}


void data::ListODE(int opt)
{
   if (opt && !nb_ode) {
      cout << "No defined ordinary differential equations." << endl;
      return;
   }
   cout << "Number of ordinary differential equations: " << nb_ode << endl;
   for (size_t i=1; i<theODE.size(); ++i) {
      if (dn[ODEName[i]].active)
         cout << *theODE[i];
   }
}


void data::ListPDE(int opt)
{
   if (opt && !nb_pde) {
      cout << "No defined partial differential equations." << endl;
      return;
   }
   cout << "Number of partial differential equations: " << nb_pde << endl;
   for (size_t i=1; i<thePDE.size(); ++i) {
      string s = PDEName[i];
      if (dn[s].active) {
         pde *e = thePDE[i];
         cout << "PDE name: " << PDEName[i] << endl;
         cout << "PDE id: " << e->eq << endl;
         cout << "PDE unknown vector(s): ";
         for (int j=0; j<e->nb_vectors-1; ++j)
            cout << e->fd[j].fn << ", ";
         cout << e->fd[e->nb_vectors-1].fn << endl;
      }
   }
}


void data::ListLS(int opt)
{
   if (opt && !nb_ls) {
      cout << "No defined linear systems." << endl;
      return;
   }
   cout << "Number of linear systems: " << nb_ls << endl;
   for (size_t i=1; i<theLS.size(); ++i) {
      string s = LSName[i];
      if (dn[s].active) {
         ls *e = theLS[i];
         cout << "Linear System name:     " << s << endl;
         cout << "Matrix:                 " << MatrixName[e->iMat] << endl;
         cout << "Right-hand side vector: " << VectorName[e->iRHS] << endl;
         cout << "Solution vector:        " << VectorName[e->iSol] << endl;
      }
   }
}


void data::ListAE(int opt)
{
   if (opt && !nb_ae) {
      cout << "No defined algebraic equations." << endl;
      return;
   }
   cout << "Number of algebraic equations: " << nb_ae << endl;
   for (size_t i=1; i<theAE.size(); ++i)
      if (dn[AEName[i]].active)
         cout << *theAE[i];
}


void data::ListOpt(int opt)
{
   if (opt && !nb_opt) {
      cout << "No defined optimization problems." << endl;
      return;
   }
   cout << "Number of optimization problems:  " << nb_opt << endl;
   for (size_t i=1; i<theOpt.size(); ++i) {
      string s = OptName[i];
      if (dn[s].active)
         cout << *theOpt[i];
   }
}


void data::ListEig(int opt)
{
   if (opt && !nb_eig) {
      cout << "No defined eigen problems." << endl;
      return;
   }
   cout << "Number of eigen problems: " << nb_eig << endl;
   for (size_t i=1; i<theEig.size(); ++i) {
      string s = EigName[i];
      if (dn[s].active) {
         eigen *e = theEig[i];
         cout << "Eigen problem name:                         " << e->name << endl;
         cout << "Matrix size:                                " << e->M->getNbRows() << "x"
              << e->M->getNbColumns() << endl;
         cout << "Vector of eigenvalues (real and imag):      " << e->evR << ", " << e->evI << endl;
         if (e->eig_vec)
            cout << "Vectors of eigenvectors (real and imag): " << e->evectR[0] << ", "
                 << e->evectI[0] << " .. " << e->evectR[e->nb_eigv-1] << ", "
                 << e->evectI[e->nb_eigv-1] << endl;
      }
   }
}


void data::Summary()
{
   if (nb_params+nb_vectors+nb_matrices+nb_fcts+nb_tabs+nb_grids+nb_meshes+
       nb_ls+nb_ae+nb_ode+nb_pde==0) {
      cout << "No defined data or entities." << endl;
      return;
   }
   cout << "\nSUMMARY OF DATA AND ENTITIES:" << endl;
   cout << "---------------------------------------------------------------" << endl;
   if (nb_params) {
      ListParams(0);
      cout << endl;
   }
   if (nb_vectors) {
      ListVectors(0);
      cout << endl;
   }
   if (nb_matrices) {
      ListMatrices(0);
      cout << endl;
   }
   if (nb_hvectors) {
      ListHVectors(0);
      cout << endl;
   }
   if (nb_fcts) {
      ListFunctions(0);
      cout << endl;
   }
   if (nb_tabs) {
      ListTabs(0);
      cout << endl;
   }
   if (nb_grids) {
      ListGrids(0);
      cout << endl;
   }
   if (nb_meshes) {
      ListMeshes(0);
      cout << endl;
   }
   if (nb_ls) {
      ListLS(0);
      cout << endl;
   }
   if (nb_ae) {
      ListAE(0);
      cout << endl;
   }
   if (nb_ode) {
      ListODE(0);
      cout << endl;
   }
   if (nb_pde) {
      ListPDE(0);
      cout << endl;
   }
   if (nb_opt) {
      ListOpt(0);
      cout << endl;
   }
   if (nb_eig)
      ListEig(0);
   cout << "---------------------------------------------------------------" << endl;
}


int data::saveGmsh(const string &file, const Vect<double>& v)
{
   int nb_en=0;
   ofstream ff(file);

   if (v.WithMesh()==false)
      return 1;
   OFELI::Mesh &ms = v.getMesh();
   size_t nb_dof = v.getNbDOF();
   ff << "View \"" << v.getName() << "\" {" << endl;
   switch (ms.getDim()) {

      case 1:
         for (size_t ne=1; ne<=ms.getNbElements(); ++ne) {
            Element &el = *ms(ne);
            ff << "SL(";
            ff << el(1)->getX() <<  ", 0., 0., " << el(2)->getX() <<  ", 0., 0. ) {" << endl;
            ff << v(el(1)->n(),1) << "," << v(el(2)->n(),1) << endl;
         }
         ff << "};" << endl;
         break;

      case 2:
         for (size_t ne=1; ne<=ms.getNbElements(); ++ne) {
            Element &el = *ms(ne);
            if (nb_dof==1)
               ff << 'S';
            else
               ff << 'V';
            size_t nb_en = el.getNbNodes();
            if (nb_en==3)
               ff << "T(";
            else if (nb_en==4)
               ff << "Q(";
            for (size_t k=1; k<nb_en; ++k)
               ff << el(k)->getX() << "," << el(k)->getY() << ",0.,";
            ff << el(nb_en)->getX() << "," << el(nb_en)->getY() << ",0.) {" << endl;
            for (size_t k=1; k<=nb_en; ++k) {
               ff << v(el(k)->n(),1);
               if (nb_dof > 1)
                  ff << "," << v(el(k)->n(),2) << ",0.0";
               if (k<nb_en)
                  ff << ",";
            }
            ff << "};" << endl;
         }
         ff << "};" << endl;
            break;

         case 3:
            for (size_t ne=1; ne<=ms.getNbElements(); ++ne) {
               Element &el = *ms(ne);
               if (nb_dof==1)
                  ff << 'S';
               else
                  ff << 'V';
               if ((nb_en=el.getNbNodes())==4)
                  ff << "S(";
               else if (nb_en==8)
                  ff << "H(";
               else if (nb_en==6)
                  ff << "I(";
               else if (nb_en==5)
                  ff << "Y(";
               for (int k=1; k<=nb_en-1; ++k)
                  ff << el(k)->getX() << ","
                     << el(k)->getY() << ","
                     << el(k)->getZ() << ",";
               ff << el(nb_en)->getX() << ","
                  << el(nb_en)->getY() << ","
                  << el(nb_en)->getZ() << ") {" << endl;
               for (int k=1; k<=nb_en; ++k) {
                  ff << v(el(k)->n(),1);
                  if (nb_dof > 1)
                     ff << "," << v(el(k)->n(),1) << ","
                        << v(el(k)->n(),3) << endl;
               }
               ff << "};" << endl;
            }
            ff << "};" << endl;
            break;
      }
      return 0;
}


int data::saveVTK(const string &file, const Vect<double>& v)
{
   static map<int,int> shCode = {{LINE,2},{TRIANGLE,3},{QUADRILATERAL,4},{TETRAHEDRON,4},
                                 {HEXAHEDRON,8},{PENTAHEDRON,6}};
   static map<int,int> ShCode = {{LINE,3},{TRIANGLE,5},{QUADRILATERAL,9},{TETRAHEDRON,10},
                                 {HEXAHEDRON,12},{PENTAHEDRON,13}};

   if (v.WithMesh()==false)
      return 1;
   size_t nb_dof = v.getNbDOF();
   OFELI::Mesh &ms = v.getMesh();
   int sz=0;
   for (size_t ne=1; ne<=ms.getNbElements(); ++ne) {
      Element &el = *ms(ne);
      sz += shCode[el.getShape()] + 1;
   }
   ofstream ff(file.c_str());
   ff << setprecision(16) << std::scientific;
   ff << "# vtk DataFile Version 2.0\n# Imported from OFELI files\nASCII" << endl;
   ff << "DATASET UNSTRUCTURED_GRID\nPOINTS " << ms.getNbNodes() << " double" << endl;
   for (auto const& n: ms.theNodes)
      ff << n->getX() << "  " << n->getY() << "  " << n->getZ() << endl;
   ff << "\nCELLS " << ms.getNbElements() << setw(10) << sz << endl;
    for (auto const& e: ms.theElements) {
      ff << setw(8) << shCode[e->getShape()];
      for (int i=1; i<=shCode[e->getShape()]; ++i)
         ff << setw(10) << (*e)(i)->n()-1;
      ff << endl;
   }
   ff << "\nCELL_TYPES  " << ms.getNbElements() << endl;
   int k=0;
   for (auto const& e: ms.theElements) {
      ff << setw(4) << ShCode[e->getShape()];
      if (++k%30 == 0)
         ff << endl;
   }
   ff << "\nPOINT_DATA  " << ms.getNbNodes() << endl;

   if (nb_dof==1)
      ff << "SCALARS  " << v.getName()<< "  double  1\nLOOKUP_TABLE  default" << endl;
   else
      ff << "VECTORS  " << v.getName() << "  double" << endl;

   for (size_t n=1; n<=ms.getNbNodes(); ++n) {
      ff << v(n,1) << " ";
      if (nb_dof>1) {
         ff << v(n,2) << " ";
         if (nb_dof > 2)
            ff << v(n,3) << " ";
         else
            ff << 0. << " ";
      }
      ff << endl;
   }
   return 0;
}


int data::saveGnuplot(const string& file, const Vect<double>& v)
{
   ofstream ff(file.c_str());
   if (_theMesh!=nullptr) {
      switch (_theMesh->getDim()) {

         case 1:
            for (size_t n=1; n<=_theMesh->getNbNodes(); ++n)
               ff << (*_theMesh)[n]->getCoord().x << "  " << v(n) << endl;
            break;

         case 2:
            break;

         case 3:
            break;

      }
      ff.close();
   }
   else if (_theGrid!=nullptr) {
      switch (_theGrid->getDim()) {

         case 1:
            for (size_t n=1; n<=_theGrid->getNbNodes(); ++n)
               ff << _theGrid->getX(n) << "  " << v(n) << endl;
            break;

         case 2:
            break;

         case 3:
            break;

      }
      ff.close();
   }
   return 0;
}


int data::saveTecplot(const string& file, const Vect<double>& v)
{
   using namespace OFELI;
   static map<int,string> shape = {{LINE,"LINESEG"},{QUADRILATERAL,"QUADRILATERAL"},{TRIANGLE,"TRIANGLE"},
                                   {TETRAHEDRON,"TETRAHEDRON"},{HEXAHEDRON,"HEXAHEDRON"},
                                   {PENTAHEDRON,"HEXAHEDRON"}};
   if (v.WithMesh()==false)
      return 1;
   size_t nb_dof = v.getNbDOF();
   Mesh &ms = v.getMesh();
   ofstream ff(file.c_str());
   ff.setf(ios::right|ios::scientific);
   ff << "TITLE = \" \"\n\nVARIABLES = \"X\", \"Y\"";
   if (ms.getDim()==3)
      ff << ", \"Z\"";
   if (nb_dof==1)
      ff << ", \"T\"";
   else if (nb_dof==2)
      ff << ", \"UX\", \"UY\"";
   else if (nb_dof==3)
      ff << ", \"UX\", \"UY\", \"UZ\"";
   ff << "\n\nZONE T=\"" << "step-1" << "\", N=" << ms.getNbNodes() << ", E="
      << ms.getNbElements() << ", F=FEPOINT, ET=" << shape[ms.getShape()]
      << ", SOLUTIONTIME=" << v.getTime();
   ff << ", D=(1,";
   if (ms.getDim()>1)
      ff << "2,";
   if (ms.getDim()==3)
      ff << "3,";
   ff << "FECONNECT)";
   ff << endl;
   for (size_t n=1; n<=ms.getNbNodes(); ++n) {
      Node &nd=*ms[n];
      for (size_t i=1; i<=ms.getDim(); ++i)
         ff << "  " << nd.getCoord(i);
      for (size_t j=0; j<nb_dof; j++)
         ff << "  " << v[nb_dof*(nd.n()-1)+j];
      ff << endl;
   }
   for (size_t ne=1; ne<=ms.getNbElements(); ++ne) {
      Element &el=*ms(ne);
     for (size_t i=1; i<=el.getNbNodes(); ++i)
         ff << setw(10) << el(i)->n();
      ff << endl;
   }
   ff.close();
   return 0;
}


int data::check_variable_name(const string& n)
{
   for (size_t i=0; i<n.length(); ++i)
      CHK_MSGR(!isalpha(n[i]) && !isdigit(n[i]) && n[i]!='_',"calc>","Illegal variable name: "+n)
   return 0;
}


int data::checkAlreadyUsed(const string& s, const DataType& dt)
{
   if (dn[s].i==0)
      return 0;
   DataType d = dn[s].dt;
   if ((dt!=DataType::MESH  && s.substr(0,3)==DEF_MESH_NAME) ||
       (dt!=DataType::GRID  && s.substr(0,3)==DEF_GRID_NAME) ||
       (dt!=DataType::LS    && s.substr(0,3)==DEF_LS_NAME)   ||
       (dt!=DataType::AE    && s.substr(0,3)==DEF_AE_NAME)   ||
       (dt!=DataType::ODE   && s.substr(0,3)==DEF_ODE_NAME)  ||
       (dt!=DataType::PDE   && s.substr(0,3)==DEF_PDE_NAME)  ||
       (dt!=DataType::OPTIM && s.substr(0,3)==DEF_OPT_NAME)  || 
       (dt!=DataType::EIGEN && s.substr(0,3)==DEF_EIG_NAME)  ||
       (dt!=DataType::FCT   && s.substr(0,3)==DEF_FCT_NAME))
      ILLEGAL_PREFIX

   if (d!=dt) {
  
      switch (d) {

         case DataType::PARAM:
            ALREADY_USED(ParamLabel,"a parameter")
            break;

         case DataType::VECTOR:
            ALREADY_USED(VectorLabel,"a vector")
            break;

         case DataType::HVECTOR:
            ALREADY_USED(HVectorLabel,"a history vector")
            break;

         case DataType::MATRIX:
            ALREADY_USED(MatrixLabel,"a matrix")
            break;

         case DataType::GRID:
            ALREADY_USED(GridLabel,"a grid")
            break;

         case DataType::MESH:
            ALREADY_USED(MeshLabel,"a mesh")
            break;

         case DataType::TAB:
            ALREADY_USED(TabLabel,"a tabulation")
            break;

         case DataType::FCT:
            ALREADY_USED(FctLabel,"a function")
            break;

         case DataType::LS:
            ALREADY_USED(LSLabel,"a linear system")
            break;

         case DataType::AE:
            ALREADY_USED(AELabel,"a parameter")
            break;

         case DataType::ODE:
            ALREADY_USED(ODELabel,"an ordinary differential equation")
            break;

         case DataType::PDE:
            ALREADY_USED(PDELabel,"a partial differential equation")
            break;

         case DataType::EIGEN:
            ALREADY_USED(EigLabel,"an eigenvalue problem")
            break;

         case DataType::OPTIM:
            ALREADY_USED(OptLabel,"an optimization problem")
            break;

         default:
            break;
      }
   }
   return 0;
}


void data::remove_temp()
{
   for (auto const& f: temp_file)
      remove(f.c_str());
}

} /* namespace RITA */