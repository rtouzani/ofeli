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

                        Definition of class 'data'

  ==============================================================================*/

#pragma once

#include <string>
#include <iostream>
using std::string;

#include "mesh.h"
#include "HVect.h"
#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Matrix.h"
#include "funct.h"
#include "io/Tabulation.h"
#include "OFELI.h"
#include "gnuplot.h"

namespace RITA {

class rita;
class configure;
class cmd;
class ls;
class ae;
class ode;
class pde;
class optim;
class eigen;


#define DEF_LS_NAME     "sys"
#define DEF_AE_NAME     "alg"
#define DEF_ODE_NAME    "ode"
#define DEF_PDE_NAME    "pde"
#define DEF_OPT_NAME    "opt"
#define DEF_EIG_NAME    "eig"
#define DEF_FCT_NAME    "fun"
#define DEF_GRID_NAME   "grd"
#define DEF_MESH_NAME   "msh"
#define DEF_VECTOR_NAME "v"
#define DEF_MATRIX_NAME "M"
#define DEF_TAB_NAME    "t"

#define DEFAULT_LS_NAME(s)     if (s=="") s = DEF_LS_NAME+to_string(iLS)
#define DEFAULT_AE_NAME(s)     if (s=="") s = DEF_AE_NAME+to_string(iAE)
#define DEFAULT_ODE_NAME(s)    if (s=="") s = DEF_ODE_NAME+to_string(iODE)
#define DEFAULT_PDE_NAME(s)    if (s=="") s = DEF_PDE_NAME+to_string(iPDE)
#define DEFAULT_OPT_NAME(s)    if (s=="") s = DEF_OPT_NAME+to_string(iOpt)
#define DEFAULT_EIG_NAME(s)    if (s=="") s = DEF_EIG_NAME+to_string(iEig)
#define DEFAULT_FCT_NAME(s)    if (s=="") s = DEF_FCT_NAME+to_string(iFct)
#define DEFAULT_GRID_NAME(s)   if (s=="") s = DEF_GRID_NAME+to_string(iGrid)
#define DEFAULT_MESH_NAME(s)   if (s=="") s = DEF_MESH_NAME+to_string(iMesh)
#define DEFAULT_VECTOR_NAME(s) if (s=="") s = DEF_VECTOR_NAME+to_string(iVector)
#define DEFAULT_MATRIX_NAME(s) if (s=="") s = DEF_MATRIX_NAME+to_string(iMatrix)
#define DEFAULT_TAB_NAME(s)    if (s=="") s = DEF_TAB_NAME+to_string(iTab)

#define CHECK_VECT(s)          if (VectorLabel.count(s)==0) { _rita->msg("","Vector "+s+" undefined."); return 1; }
#define CHECK_ACTIVE(s)        if (!dn[s].active) { _rita->msg("","Reference to non-active data "+s); return 2; }
#define CHECK_HVECT(s)         if (HVectorLabel.count(s)==0) { _rita->msg("","History vector "+s+" undefined."); return 1; }
#define CHECK_MATRIX(s)        if (MatrixLabel.count(s)==0) { _rita->msg("","Matrix "+s+" undefined."); return 1; }
#define CHECK_TAB(s)           if (TabLabel.count(s)==0) { _rita->msg("","Tabulation "+s+" undefined."); return 1; }
#define ILLEGAL_PREFIX         { _rita->msg("","Illegal name prefix for entity "+s); return 1; }

enum class DataType { NOTHING, PARAM, VECTOR, HVECTOR, MATRIX, GRID, MESH, TAB, FCT,
                      LS, AE, ODE, PDE, OPTIM, EIGEN };

struct pls {
   double mx, Mx, my, My;
   int lx, ly, component, contour;
   string title, xaxist, yaxist, mark, label, soft;
};

class data
{

 public:

    enum class eqType { LS, AE, ODE, PDE, OPT, EIGEN, INTEGR };
    enum class DataSize { GIVEN_SIZE, GRID, NODES, ELEMENTS, SIDES, EDGES };
    enum class Storage { DENSE, SPARSE, SKYLINE, BAND, TRIDIAGONAL, DIAGONAL };
    struct Dat { int i; DataType dt; bool active;
                 Dat() { dt=DataType::NOTHING; active=false; }
               };

    data(rita *r, cmd *command, configure *config);
    ~data() { }
    int remove(const string& name);
    int rename(const string& data1, const string& data2);
    int addVector(string name, double t=0., int n=1, string file="", int opt=1);
    int addVector(OFELI::Vect<double>* v, string& name, bool opt=true);
    int addMeshVector(const string& nm1, string& nm2, DataSize s, int ndof=1, bool opt=true);
    int addGridVector(const string& nm1, string& nm2, int ndof=1, bool opt=true);
    int addParam(string name, double value, bool opt=true);
    int addLS(ls *ls, string& name);
    int addAE(ae *e, string& name);
    int addODE(ode *e, string& name);
    int addPDE(pde *e, string& name);
    int addTab(OFELI::Tabulation *tab, string& name);
    int addOpt(optim *o, string& name);
    int addEig(eigen *e, string& name);
    void setTime(string s, double t);
    DataType getType(string s);
    int checkParam(const string& name, double& value);
    int checkParam(const string& name, int& value);
    int setVectorValue(int k, const OFELI::Vect<double>* v);
    int save();
    void setVerbose(int verb) { _verb = verb; }
    void setSave(int s) { _sr = s; }
    int getVerbose() const { return _verb; }
    void set(cmd* command) { _cmd = command; }
    int ret() const { return _ret; }
    int setDataExt(int key);
    int addFunction(string& name);
    int addFunction(string& name, const string& var, const string& exp);
    int addFunction(string& name, const vector<string>& var, const string& exp, const vector<size_t>& n);
    int addFunction(string& name, const vector<string>& var, const string& exp);
    int addFunction(OFELI::Fct& f);
    int addMatrix(string name, int nr, int nc, string file="", string s="dense", bool opt=true);
    int addGrid(OFELI::Grid* g, string& name);
    int addMesh(OFELI::Mesh* ms, string& name);
    int getPar(int n, const string& msg, double& v);
    int getPar(int n, const string& msg, int& v);
    int print(const string &s);
    int getNbEq() const { return nb_ae+nb_ode+nb_pde; }
    int setTab2GridVector(OFELI::Tabulation* tab, string& name);
    void setTab2Vector(OFELI::Tabulation* tab, OFELI::Vect<double>& v);
    int add2History(const string& s, OFELI::Vect<double>& v, double t=0.0);
    int setHistory(const string& s, size_t n);
    int setPhaseHistory(const string& name);
    int checkAlreadyUsed(const string& s, const DataType& dt);
    void Summary();

  
  /*
   *  Each entity E:
   *  - iE is the index in the global array
   *  - theE[iE] is the pointer to entity number iE
   *  - NameE[iE] is the name of entity number iE
   *  - EName[name] is the index of the entity named name in the global array
   * 
   *  Entities are:
   *  Vector, Matrix, Mesh, Grid, Tab, Param, Fct, AE, ODE, PDE, ls
   */
    string p2s;
    DataType pt2s;
    int nb_vectors, nb_hvectors, nb_fcts, nb_tabs, nb_meshes, nb_grids, nb_params, nb_matrices;
    int nb_pde, nb_ode, nb_ae, nb_ls, nb_opt, nb_int, nb_eig, nb_eq;
    vector<int> nb_dof;
    int iParam, iVector, iHVector, iMatrix, iGrid, iMesh, iTab, iFct, iLS, iAE, iODE, iPDE, iOpt, iEig;
    vector<OFELI::Vect<double> *> theVector;
    vector<HVect *> theHVector;
    map<string,string> vect_hist;
    vector<double> theParam, VectorTime;
    map<string,int> VectorLabel, HVectorLabel, ParamLabel, MatrixLabel, TabLabel, GridLabel, MeshLabel;
    map<string,int> FctLabel, LSLabel, AELabel, ODELabel, PDELabel, OptLabel, EigLabel;
    map<string,string> Desc;
    map<string,Dat> dn;
    vector<DataSize> VectorSizeType;
    vector<int> VectorEquation;
    vector<OFELI::Grid *> theGrid;
    vector<OFELI::Mesh *> theMesh;
    vector<OFELI::Tabulation *> theTab;
    vector<funct *> theFct;
    vector<OFELI::Matrix<double> *> theMatrix;
    vector<ae *> theAE;
    vector<ode *> theODE;
    vector<ls *> theLS;
    vector<pde *> thePDE;
    vector<optim *> theOpt;
    vector<eigen *> theEig;
    vector<string> LSName, AEName, ODEName, PDEName, GridName, TabName, MeshName, FctName, VectorName;
    vector<string> HVectorName, ParamName, MatrixName, OptName, EigName, Names, temp_file;
    double obj, integral;
    bool ok;
    vector<eqType> eq_type;
    vector<int> eqq;
    bool FindName(const string& s) const;
    int nbAllParam() const { return theParam.size()-1; }
    int nbAllVector() const { return theVector.size()-1; }
    int nbAllHVector() const { return theHVector.size()-1; }
    int nbAllMatrix() const { return theMatrix.size()-1; }
    int nbAllTab() const { return theTab.size()-1; }
    int nbAllGrid() const { return theGrid.size()-1; }
    int nbAllMesh() const { return theMesh.size()-1; }
    int nbAllFct() const { return theFct.size()-1; }
    int nbAllLS() const { return theLS.size()-1; }
    int nbAllAE() const { return theAE.size()-1; }
    int nbAllODE() const { return theODE.size()-1; }
    int nbAllPDE() const { return thePDE.size()-1; }
    int nbAllOpt() const { return theOpt.size()-1; }
    int nbAllEig() const { return theEig.size()-1; }
    int setDesc(const string& name, const string& desc);
//    void setVar(const string& v, vector<string>& w, int n, int opt=0);
//    void setVar(const vector<string>& v, vector<string>& w, int n, int opt=0);
    int setVector();
    int setMatrix();
    int setGrid();
    int setTab();
    int setFunction();
    int setDerivative();
    int setSample();
    void ListParams(int opt);
    void ListGrids(int opt);
    void ListMeshes(int opt);
    void ListVectors(int opt);
    void ListHVectors(int opt);
    void ListFunctions(int opt);
    void ListTabs(int opt);
    void ListMatrices(int opt);
    void ListLS(int opt);
    void ListAE(int opt);
    void ListODE(int opt);
    void ListPDE(int opt);
    void ListOpt(int opt);
    void ListEig(int opt);

   int saveGnuplot(const string& file, const Vect<double>& v);
   int saveGmsh(const string &file, const Vect<double>& v);
   int saveVTK(const string &file, const Vect<double>& v);
   int saveTecplot(const string &file, const Vect<double>& v);
   void remove_temp();

    map<DataType,string> type_st = {{DataType::NOTHING,""},
                                    {DataType::PARAM,"parameter"},
                                    {DataType::VECTOR,"vector"},
                                    {DataType::HVECTOR,"history vector"},
                                    {DataType::MATRIX,"matrix"},
                                    {DataType::GRID,"grid"},
                                    {DataType::MESH,"mesh"},
                                    {DataType::TAB,"tabulation"},
                                    {DataType::FCT,"function"},
                                    {DataType::LS,"linear system"},
                                    {DataType::AE,"algrbraic equation"},
                                    {DataType::ODE,"ode"},
                                    {DataType::PDE,"pde"},
                                    {DataType::OPTIM,"optimization problem"},
                                    {DataType::EIGEN,"eigenvalue problem"}};

    map<string,OFELI::MatrixType> st_map = {{"dense",OFELI::DENSE},
                                            {"band",OFELI::BAND},
                                            {"skyline",OFELI::SKYLINE},
                                            {"sparse",OFELI::SPARSE}};

    map<string,int> ff = {{"gmsh",GMSH},
                          {"gnuplot",GNUPLOT},
                          {"vtk",VTK},
                          {"tecplot",TECPLOT},
                          {"matlab",MATLAB},
                          {"ofeli",OFELI_FF}};

    map<string,OFELI::Iteration> Ls = {{"direct",OFELI::DIRECT_SOLVER},
                                       {"cg",OFELI::CG_SOLVER},
                                       {"cgs",OFELI::CGS_SOLVER},
                                       {"bicg",OFELI::BICG_SOLVER},
                                       {"bicg-stab",OFELI::BICG_STAB_SOLVER},
                                       {"gmres",OFELI::GMRES_SOLVER}};

    map<string,OFELI::Preconditioner> Pr = {{"ident",OFELI::IDENT_PREC},
                                            {"diag",OFELI::DIAG_PREC},
                                            {"dilu",OFELI::DILU_PREC},
                                            {"ilu",OFELI::ILU_PREC},
                                            {"ssor",OFELI::SSOR_PREC}};

    map<string,NonLinearIter> NLs = {{"bisection",BISECTION},
                                     {"regula-falsi",REGULA_FALSI},
                                     {"picard",PICARD},
                                     {"secant",SECANT},
                                     {"newton",NEWTON}};

   const string Fct_help = "name, variable, nb, definition";
   const string Fct_Help = "function [name=f] [variable=x] [nb=n] [definition=d]\n"
                           "f: Name to give to function\n"
                           "x: Name of variable\n"
                           "n: Number of variables. If many variables and variable name is x, actual variables are x1, x2, ...\n"
                           "d: Regular expression defining function\n";
                   

template<class T1_,class T2_>
bool contains(const map<T1_,T2_>& m, const T1_& key)
{
   for (auto x: m) {
      if (x.first==key)
         return true;
   }
   return false;
}


 private:

    rita *_rita;
    int _size, _nb_args, _verb, _ret, _nb, _sr;
    configure *_configure;
    cmd *_cmd;
    vector<string> _pl_cmd;
    string _pr;

    OFELI::Mesh *_theMesh;
    OFELI::Tabulation *_theTab;
    OFELI::Grid *_theGrid;
    funct *_theFct;
    OFELI::Vect<double> *_theVector;
    HVect *_theHVector;
    OFELI::Matrix<double> *_theMatrix;
    Vect<double> *_u;
    HVect *_hu;
    pls _pl;
    int setNbDOF();
    int plot();
    int plot_fct();
    int plot_hist();
    int plot_vect();
    int plot_tab();
    int plot_mesh();
    void modify_exp(const string& s1, string& s2);
    int setConfigure();
    vector<double> _xv;
    int check_variable_name(const string& n);
};

} /* namespace RITA */
