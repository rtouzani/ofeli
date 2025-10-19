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

             Implementation of class 'calc' and other related classes

  ==============================================================================*/


#define _USE_MATH_DEFINES 
#undef __STRICT_ANSI__

#if defined(_WIN32)
// Memory leak dumping
  #if defined(_DEBUG)
    #define _CRTDBG_MAP_ALLOC
    #include <stdlib.h>
    #include <crtdbg.h>
    #define CREATE_LEAKAGE_REPORT
  #endif

// Needed for windows console UTF-8 support
  #include <fcntl.h>
  #include <io.h>
#endif

#include "calc.h"
#include "cmd.h"
#include "data.h"
#include "helps.h"

namespace RITA {

void FctFun::Eval(ptr_val_type&       ret,
                  const ptr_val_type* a_pArg,
                  int                 argc)
{
   if (argc==0) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   _arg.clear();
   for (int i=0; i<argc; ++i)
      _arg.push_back(a_pArg[i]->GetFloat());
   *ret = (*_fun)(_arg);
}


void FctMatrix::Eval(ptr_val_type&       ret,
                     const ptr_val_type* a_pArg,
                     int                 argc)
{
   if (argc<1 || argc>2) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type m = a_pArg[0]->GetInteger(),
            n = (argc==1) ? m : a_pArg[1]->GetInteger();
   if (m==n && n==1)
      *ret = 0.0;
   else
      *ret = matrix_type(m, n, 0.0);
}


void FctMatrixNorm2::Eval(ptr_val_type&       ret,
                          const ptr_val_type* a_pArg,
                          int                 argc)
{
   if (argc!=1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type nr = a_pArg[0]->GetRows(), nc = a_pArg[0]->GetCols();
   float_type x = 0.0, z = 0.0;
   for (int i=0; i<nr; ++i) {
      for (int j=0; j<nc; ++j) {
         z = a_pArg[0]->At(i,j).GetFloat();
         x += z*z;
      }
   }
   *ret = sqrt(x);
}


void FctMatrixNormMax::Eval(ptr_val_type&       ret,
                            const ptr_val_type* a_pArg,
                            int                 argc)
{
   if (argc!=1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type nr = a_pArg[0]->GetRows(), nc = a_pArg[0]->GetCols();
   float_type x = 0.0, z = 0.0;
   for (int i=0; i<nr; ++i) {
      for (int j=0; j<nc; ++j) {
         z = a_pArg[0]->At(i,j).GetFloat();
         x = std::max(x,std::abs(z));
      }
   }
   *ret = x;
}


void FctVector::Eval(ptr_val_type&       ret,
                     const ptr_val_type* a_pArg,
                     int                 argc)
{
   if (argc != 1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type n = a_pArg[0]->GetInteger();
   if (n==1)
      *ret = 0.0;
   else
      *ret = matrix_type(n, 1, 0.0);
}


void FctRowVector::Eval(ptr_val_type&       ret,
                        const ptr_val_type* a_pArg,
                        int                 argc)
{
   if (argc != 1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type n = a_pArg[0]->GetInteger();
   if (n==1)
      *ret = 0.0;
   else
      *ret = matrix_type(1, n, 0.0);
}


void FctRowVectorNorm1::Eval(ptr_val_type&       ret,
                             const ptr_val_type* a_pArg,
                             int                 argc)
{
   if (argc!=1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type nr = a_pArg[0]->GetRows(), nc = a_pArg[0]->GetCols();
   float_type x = 0.0, z = 0.0;
   for (int i=0; i<nr; ++i) {
      for (int j=0; j<nc; ++j) {
         z = a_pArg[0]->At(i,j).GetFloat();
         x += std::abs(z);
      }
   }
   *ret = x;
}


void FctRowVectorNorm2::Eval(ptr_val_type&       ret,
                             const ptr_val_type* a_pArg,
                             int                 argc)
{
   if (argc!=1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type nr = a_pArg[0]->GetRows(), nc = a_pArg[0]->GetCols();
   float_type x = 0.0, z = 0.0;
   for (int i=0; i<nr; ++i) {
      for (int j=0; j<nc; ++j) {
         z = a_pArg[0]->At(i,j).GetFloat();
         x += z*z;
      }
   }
   *ret = sqrt(x);
}


void FctRowVectorNormMax::Eval(ptr_val_type&       ret,
                               const ptr_val_type* a_pArg,
                               int                 argc)
{
   if (argc!=1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type nr = a_pArg[0]->GetRows(), nc = a_pArg[0]->GetCols();
   float_type x = 0.0, z = 0.0;
   for (int i=0; i<nr; ++i) {
      for (int j=0; j<nc; ++j) {
         z = a_pArg[0]->At(i,j).GetFloat();
         x = std::max(x,std::abs(z));
      }
   }
   *ret = x;
}


void FctRowVectorCanonical::Eval(ptr_val_type&       ret,
                                 const ptr_val_type* a_pArg,
                                 int                 argc)
{
   if (argc != 2) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type n = a_pArg[0]->GetInteger(),
            i = a_pArg[1]->GetInteger();
   if (i>n) {
      ErrorContext err;
      err.Errc = ecARRAY_SIZE_MISMATCH;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   matrix_type M(1,int(n),0.0);
   M.At(i-1) = 1.0;
   *ret = M;
}


void FctVectorNorm1::Eval(ptr_val_type&       ret,
                          const ptr_val_type* a_pArg,
                          int                 argc)
{
   if (argc!=1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type nr = a_pArg[0]->GetRows(), nc = a_pArg[0]->GetCols();
   float_type x = 0.0, z = 0.0;
   for (int i=0; i<nr; ++i) {
      for (int j=0; j<nc; ++j) {
         z = a_pArg[0]->At(i,j).GetFloat();
         x += std::abs(z);
      }
   }
   *ret = x;
}


void FctVectorNorm2::Eval(ptr_val_type&       ret,
                          const ptr_val_type* a_pArg,
                          int                 argc)
{
   if (argc!=1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type nr = a_pArg[0]->GetRows(), nc = a_pArg[0]->GetCols();
   float_type x = 0.0, z = 0.0;
   for (int i=0; i<nr; ++i) {
      for (int j=0; j<nc; ++j) {
         z = a_pArg[0]->At(i,j).GetFloat();
         x += z*z;
      }
   }
   *ret = sqrt(x);
}


void FctVectorNormMax::Eval(ptr_val_type&       ret,
                            const ptr_val_type* a_pArg,
                            int                 argc)
{
   if (argc!=1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type nr = a_pArg[0]->GetRows(), nc = a_pArg[0]->GetCols();
   float_type x = 0.0, z = 0.0;
   for (int i=0; i<nr; ++i) {
      for (int j=0; j<nc; ++j) {
         z = a_pArg[0]->At(i,j).GetFloat();
         x = std::max(x,std::abs(z));
      }
   }
   *ret = x;
}


void FctVectorCanonical::Eval(ptr_val_type&       ret,
                              const ptr_val_type* a_pArg,
                              int                 argc)
{
   if (argc != 2) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type n = a_pArg[0]->GetInteger(),
            i = a_pArg[1]->GetInteger();
   if (i>n) {
      ErrorContext err;
      err.Errc = ecARRAY_SIZE_MISMATCH;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   matrix_type M(int(n),1,0.0);
   M.At(i-1) = 1.0;
   *ret = M;
}


void FctMatrixLaplace1D::Eval(ptr_val_type&       ret,
                              const ptr_val_type* a_pArg,
                              int                 argc)
{
   if (argc<1 || argc>2) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type n = a_pArg[0]->GetInteger();
   double c = 1.0;
   if (argc==2)
      c = a_pArg[1]->GetFloat();
   matrix_type M(n,n,0.0);
   for (int i=0; i<n; ++i)
      M.At(i,i) = 2.0*c;
   for (int i=0; i<n-1; ++i)
      M.At(i,i+1) = -c;
   for (int i=1; i<n; ++i)
      M.At(i,i-1) = -c;
   *ret = M;
}


void FctMatrixDiag::Eval(ptr_val_type&       ret,
                         const ptr_val_type* a_pArg,
                         int                 argc)
{
   if (argc<1 || argc>2) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type n = a_pArg[0]->GetInteger();
   double d = 1.0;
   if (argc==2)
      d = a_pArg[1]->GetFloat();
   matrix_type M(n,n,0.0);
   for (int i=0; i<n; ++i)
      M.At(i,i) = d;
   *ret = M;
}

/* ----------------------------------------------------------------------------------------------- */
/*                             Implementation of calc member functions                             */
/* ----------------------------------------------------------------------------------------------- */

calc::calc(data*     d,
           cmd*      command,
           ofstream* ofh, 
           ofstream* ofl)
     : _data(d), _cmd(command), _ofh(ofh), _ofl(ofl), _nb_def_funs(0)
{
   theParser = new ParserX(pckALL_NON_COMPLEX);
   theParser->EnableAutoCreateVar(true);
   theParser->DefineFun(new FctMatrix);
   theParser->DefineFun(new FctMatrixNorm2);
   theParser->DefineFun(new FctMatrixNormMax);
   theParser->DefineFun(new FctMatrixLaplace1D);
   theParser->DefineFun(new FctMatrixDiag);
   theParser->DefineFun(new FctVector);
   theParser->DefineFun(new FctVectorCanonical);
   theParser->DefineFun(new FctVectorNorm1);
   theParser->DefineFun(new FctVectorNorm2);
   theParser->DefineFun(new FctVectorNormMax);
   theParser->DefineFun(new FctRowVector);
   theParser->DefineFun(new FctRowVectorCanonical);
   theParser->DefineFun(new FctRowVectorNorm1);
   theParser->DefineFun(new FctRowVectorNorm2);
   theParser->DefineFun(new FctRowVectorNormMax);
}


void calc::ListExprVar()
{
   cout << "\nVariables found in : \"" << theParser->GetExpr() << "\"\n";
   cout << "-----------------------------\n";
   var_maptype vmap = theParser->GetExprVar();
   if (!vmap.size())
      cout << "Expression does not contain variables" << endl;
	else {
      for (const auto& v: vmap)
         cout << "  " << v.first << " = " << (Variable&)(*(v.second)) << endl;
	}
}


void calc::ListFun()
{
   if (_nb_def_funs==0) {
      cout << "No defined functions found.\n";
      return;
   }
   cout << "List of defined functions:" << endl;
   for (const auto& f: _funs)
      cout << f << endl;
}


void calc::ListVar()
{
   var_maptype var = theParser->GetVar();
   if (var.size()==0) {
      cout << "No variables defined.\n";
      return;
   }
   if (var.size()==1)
      cout << "One defined variable:\n";
   else
      cout << "\n" << var.size() << " defined variables:\n";
   for (const auto& v: var)
      cout << v.first << " = " << *v.second << std::endl;
}


void calc::ListConst()
{
   cout << "Parser constants:\n";
   var_maptype cmap = theParser->GetConst();
   if (!cmap.size())
      cout << "Expression does not contain constants\n";
	else {
      for (const auto& v: cmap)
         cout << "  " << v.first << " = " << *v.second << std::endl;
   }
}


int calc::setParam(const string& name,
                   double        v)
{
   Value *V = new Value(v);
   theParser->DefineVar(name,Variable(V));
   return 0;
}


int calc::setFun(funct*                f,
                 const string&         exp,
                 const vector<string>& var)
{
   theParser->SetExpr(exp);
   var_maptype vmap = theParser->GetExprVar();
   if (vmap.size()==0) {
      cout << "Expression does not contain variables" << endl;
      return 1;
   }
   f->setVar(var);
   f->par_value.clear();
   for (const auto& p: vmap) {
      for (size_t i=1; i<_data->theParam.size(); ++i) {
         if (p.first==_data->ParamName[i]) {
            f->setPar(p.first);
            Variable &w = (Variable&)(*(p.second));
            f->par_value.push_back(w.GetFloat());
         }
      }
   }
   f->setExpr(exp);
   if (std::find(_funs.begin(),_funs.end(),f->name) != _funs.end())
      ;
   else {
      theParser->DefineFun(new FctFun(f));
      _funs.push_back(f->name);
      _nb_def_funs++;
   }
   return 0;
}


int calc::setFun(funct*                f,
                 const string&         exp,
                 const vector<string>& var,
                 const vector<size_t>& n)
{
   theParser->SetExpr(exp);
   var_maptype vmap = theParser->GetExprVar();
   if (vmap.size()==0) {
      cout << "Expression does not contain variables" << endl;
      return 1;
   }
   f->setVar(var,n);
   f->par_value.clear();
   for (const auto& p: vmap) {
      for (size_t i=1; i<_data->theParam.size(); ++i) {
         if (p.first==_data->ParamName[i]) {
            f->setPar(p.first);
            Variable &w = (Variable&)(*(p.second));
            f->par_value.push_back(w.GetFloat());
         }
      }
   }
   f->setExpr(exp);
   if (std::find(_funs.begin(),_funs.end(),f->name) != _funs.end())
      ;
   else {
      theParser->DefineFun(new FctFun(f));
      _funs.push_back(f->name);
      _nb_def_funs++;
   }
   return 0;
}


int calc::setFun(funct* f)
{
   theParser->SetExpr(f->getExpression());
   var_maptype vmap = theParser->GetExprVar();
   if (vmap.size()==0) {
      cout << "Expression does not contain variables" << endl;
      return 1;
   }
	else {
      size_t np = _data->theParam.size() - 1;
      for (const auto& p: vmap) {
         for (size_t i=1; i<np; ++i) {
            if (p.first==_data->ParamName[i]) {
               f->nb_par++;
               f->par.push_back(p.first);
            }
         }
      }
	}

   if (std::find(_funs.begin(),_funs.end(),f->name) != _funs.end())
      ;
   else {
      theParser->DefineFun(new FctFun(f));
      _funs.push_back(f->name);
      _nb_def_funs++;
   }
   return 0;
}


int calc::setVector(OFELI::Vect<double>* u)
{
   size_t n = u->size();
   Value *v = new Value(n,1);
   for (size_t i=0; i<n; ++i)
      v->At(i,0) = (*u)[i];
   theParser->DefineVar(u->getName(),Variable(v));
   return 0;
}


int calc::setVectorValue(const string&              name,
                         const OFELI::Vect<double>* u)
{
   var_maptype var = theParser->GetVar();
   if (var.count(name)==0)
      return 1;
   else {
      Variable &w = (Variable&)(*(var.at(name)));
      for (size_t i=0; i<u->size(); ++i)
         w.At(i,0) = (*u)[i];
      return 0;
   }
}


int calc::setMatrixValue(const string&          name,
                         OFELI::Matrix<double>* M)
{
   var_maptype var = theParser->GetVar();
   if (var.count(name)==0)
      return 1;
   else {
      Variable &w = (Variable&)(*(var.at(name)));
      for (size_t i=0; i<M->getNbRows(); ++i)
         for (size_t j=0; j<M->getNbColumns(); ++j)
            w.At(i,j) = M->at(i+1,j+1);
      return 0;
   }
}


int calc::setMatrix(OFELI::Matrix<double>* M)
{
   int nr=M->getNbRows(), nc=M->getNbColumns();
   Value *v = new Value(nr,nc,0.);
   for (int i=1; i<=nr; ++i)
      for (int j=1; j<=nc; ++j)
         v->At(i-1,j-1) = M->at(i,j);
   theParser->DefineVar(M->getName(),Variable(v));
   return 0;
}


int calc::CheckKeywords(const string& sLine)
{
   int ret = 0;
   if (sLine=="data") {
      ListConst();
      ListVar();
      ListFun();
      ret = 1;
   }
   else if (sLine=="list_const")
      ListConst(), ret = 1;
   else if (sLine=="list_var")
      ListVar(), ret = 1;
   else if (sLine=="list_expr_var")
      ListExprVar(), ret = 1;
   else if (sLine=="list_func")
      ListFun(), ret = 1;
   else if (sLine=="help")
      cout << Calc_Help << endl, ret = 1;
   return ret;
}


int calc::run()
{
   try {
      string sLine = _cmd->buffer();
      *_ofh << sLine << endl;
      if (sLine[sLine.size()-1]==';')
         sLine.pop_back();
      switch (CheckKeywords(sLine))
      {
         case  0: break;
         case  1: return 1;
         case -1: return -1;
      }
      theParser->SetExpr(sLine);
      Value val = theParser->Eval();
      const var_maptype& expr_var = theParser->GetExprVar();
      const var_maptype& var = theParser->GetVar();
      if (sLine.find('=') == string::npos) {
         if (val.GetType()=='v') {
            cout << "Unknown variable !" << endl;
            return 1;
         }
         cout << " = " << val << endl;
         return 0;
      }

      for (auto const& v: expr_var) {
         Variable &w = (Variable&)(*(v.second));
         string name = v.first;
         var_maptype::const_iterator it = var.find(name);
         switch (w.GetType()) {

            case 'i':
            case 'f':
               if (it==var.end())
                  _data->setParam(name,w.GetFloat());
               else
                  _data->addParam(name,w.GetFloat());
//               cout << name << " = " << val << endl;
               break;

            case 'm':
               {
                  int nr=w.GetRows(), nc=w.GetCols(), n=std::max(nr,nc);
                  if (nr==1 || nc==1) {
                     if (it==var.end()) {
                        OFELI::Vect<double> *v = new OFELI::Vect<double>(n);
                        for (int i=0; i<n; ++i)
                           (*v)[i] = w.GetArray().At(0,i).GetFloat();
                        _data->setVector(name,v);
                     }
                     else {
                        int k = _data->addVector(name,0.,n,"");
                        if (k==0)
                           break;
                        OFELI::Vect<double> *v = _data->theVector[k];
                        for (int i=0; i<n; ++i)
                           (*v)[i] = w.GetArray().At(0,i).GetFloat();
                     }
                  }
                  else {
                     if (it==var.end()) {
                        OFELI::DMatrix<double> *M = new OFELI::DMatrix<double>(nr,nc);
                        for (int i=0; i<nr; ++i)
                           for (int j=0; j<nc; ++j)
                              (*M)(i+1,j+1) = w.GetArray().At(i,j).GetFloat();
                        _data->setMatrix(name,M);
                     }
                     else {
                        int k = _data->addMatrix(name,nr,nc,"","dense");
                        OFELI::Matrix<double> *M = _data->theMatrix[k];
                        for (int i=0; i<nr; ++i)
                           for (int j=0; j<nc; ++j)
                              (*M)(i+1,j+1) = w.GetArray().At(i,j).GetFloat();
                     }
                  }
//                  cout << name << " = " << *it->second.Get() << endl;
                  break;
               }

            case 'v':
               theParser->RemoveVar(name);
               if (_data->dn[name].active==-1) {
                  cout << "No data or entity " << name << " found." << endl;
                  return 1;
               }
               cout << name << " = " << val << endl;
               break;

            default:
               break;
         }
      }
   }
   catch(ParserError &e) {
      if (e.GetPos()!=-1) {
         string_type sMarker;
         sMarker.insert(0,e.GetPos()+10,' ');
         sMarker += "^\n";
         cout << sMarker;
      }
      cout << "Error: " << e.GetMsg() << std::dec << endl;
      *_ofl << "In calc: " << e.GetMsg() << std::dec << endl;
   }
   return 0;
}


int calc::get_int(const string& s)
{
   theParser->SetExpr(s);
   return theParser->Eval().GetInteger();
}


double calc::get_double(const string& s)
{
   theParser->SetExpr(s);
   return theParser->Eval().GetFloat();
}

} /* namespace RITA */