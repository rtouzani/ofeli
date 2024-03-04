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
#include "rita.h"
#include "cmd.h"
#include "data.h"

namespace RITA {

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


void FctMatrixNormI::Eval(ptr_val_type&       ret,
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


void calc::ListExprVar()
{
   cout << "\nVariables found in : \"" << theParser.GetExpr() << "\"\n";
   cout << "-----------------------------\n";
   var_maptype vmap = theParser.GetExprVar();
   if (!vmap.size())
      cout << "Expression does not contain variables" << endl;
	else {
      for (const auto& v: vmap)
         cout << "  " << v.first << " = " << (Variable&)(*(v.second)) << endl;
	}
}


void calc::ListVar()
{
   var_maptype variables = theParser.GetVar();
   if (!variables.size())
      return;
   cout << "\nParser variables:\n";
   cout << "-----------------\n";
   cout << "Number: " << variables.size() << "\n";
   for (const auto& v: variables)
      cout << "Name: " << v.first << " = " << *v.second << std::endl;
}


void calc::ListConst()
{
   cout << "\nParser constants:\n";
   cout << "-----------------\n";
   var_maptype cmap = theParser.GetConst();
   if (!cmap.size())
      cout << "Expression does not contain constants\n";
	else {
      for (const auto& v: cmap)
         cout << "  " << v.first << " = " << *v.second << std::endl;
   }
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


void FctVectorNormI::Eval(ptr_val_type&       ret,
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


void FctMatrixLaplace1D::Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int argc)
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


void FctMatrixDiag::Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int argc)
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

calc::calc(rita* r,
           cmd*  command)
     : _rita(r), _cmd(command)
{
   theParser.EnableAutoCreateVar(true);
   theParser.DefineFun(new FctMatrix);
   theParser.DefineFun(new FctMatrixNorm2);
   theParser.DefineFun(new FctMatrixNormI);
   theParser.DefineFun(new FctMatrixLaplace1D);
   theParser.DefineFun(new FctMatrixDiag);
   theParser.DefineFun(new FctVector);
   theParser.DefineFun(new FctVectorCanonical);
   theParser.DefineFun(new FctVectorNorm1);
   theParser.DefineFun(new FctVectorNorm2);
   theParser.DefineFun(new FctVectorNormI);
}


int calc::setParam(const string& name,
                   double        v)
{
   _v = new Value(v);
   _theV.push_back(_v);
   theParser.DefineVar(name,Variable(_v));
   return 0;
}


int calc::setVector(OFELI::Vect<double>* u)
{
   int n = u->size();
   _v = new Value(n,1,0.);
   _theV.push_back(_v);
   for (int i=0; i<n; ++i)
      _v->At(i) = (*u)[i];
   theParser.DefineVar(u->getName(),Variable(_v));
   return 0;
}


int calc::setVectorValue(const string&              name,
                         const OFELI::Vect<double>* u)
{
   var_maptype var = theParser.GetVar();
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
   var_maptype var = theParser.GetVar();
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
   _v = new Value(nr,nc,0.0);
   _theV.push_back(_v);
   theParser.DefineVar(M->getName(),Variable(_v));
   for (int i=1; i<=nr; ++i)
      for (int j=1; j<=nc; ++j)
         _v->At(i-1,j-1) = M->at(i,j);
   return 0;
}


int calc::CheckKeywords()
{
   if (_sLine.substr(0,5)=="print") {
      _rita->msg("","Command print is not allowed in this mode","");
      return 0;
   }
   else if (_sLine=="list") {
      ListConst();
      ListVar();
      ListExprVar();
      return 1;
   }
   return 0;
}


int calc::run()
{
   _data = _rita->_data;
   try {
      _sLine = _cmd->buffer();
      *_rita->ofh << _sLine << endl;
      if (_sLine[0]=='%')
         _sLine.erase(0,1);
      if (_sLine[_sLine.size()-1]==';')
         _sLine.pop_back();
      switch (CheckKeywords())
      {
         case  0: break;
         case  1: break;
         case -1: return -1;
      }
      parse();
   }
   catch(ParserError &e) {
      if (e.GetPos()!=-1) {
         string_type sMarker;
         sMarker.insert(0,e.GetPos()+10,' ');
         sMarker += "^\n";
         cout << sMarker;
      }
      cout << "Error: " << e.GetMsg() << std::dec << endl;
      *_rita->ofl << "In " << sPrompt << e.GetMsg() << std::dec << endl;
   }
   return 0;
}


int calc::get_int(const string& s)
{
   theParser.SetExpr(s);
   return theParser.Eval().GetInteger();
}


double calc::get_double(const string& s)
{
   theParser.SetExpr(s);
   return theParser.Eval().GetFloat();
}


int calc::run(const string& buffer)
{
   _data = _rita->_data;
   try {
      _sLine = buffer;
      switch (CheckKeywords())
      {
         case  0: break;
         case  1: break;
         case -1: return -1;
      }
      parse();
   }
   catch(ParserError &e) {
      if (e.GetPos()!=-1) {
         string_type sMarker;
         sMarker.insert(0,e.GetPos()+10,' ');
         sMarker += "^\n";
         cout << sMarker;
      }
      cout << "Error: " << e.GetMsg() << std::dec << endl;
      *_rita->ofl << "In " << sPrompt << e.GetMsg() << std::dec << endl;
   }
   return 0;
}


int calc::parse()
{
   int ret = 0;
   theParser.SetExpr(_sLine);
   cout << "calc> " << _sLine << endl;
   theParser.Eval();
   const var_maptype &expr_var = theParser.GetExprVar();
   for (auto const& v: expr_var) {
      Variable &w = (Variable&)(*(v.second));
      switch (w.GetType()) {

         case 'i':
         case 'f':
            ret = _data->addParam(v.first,w.GetFloat(),true);
            if (!ret)
               return 1;
         break;

         case 'm':
            {
               int nr=w.GetRows(), nc=w.GetCols();
               if (nr==1 || nc==1) {
                  int k = _data->addVector(v.first,0.,std::max(nr,nc),"",true);
                  if (k==0) {
                     ret = 1;
                     break;
                  }
                  OFELI::Vect<double> *u = _data->theVector[k];
                  for (int i=0; i<std::max(nr,nc); ++i)
                     (*u)[i] = w.GetArray().At(0,i).GetFloat();
               }
               else {
                  int k = _data->addMatrix(v.first,nr,nc,"","dense",true);
                  OFELI::Matrix<double> *M = _data->theMatrix[k];
                  for (int i=0; i<nr; ++i)
                     for (int j=0; j<nc; ++j)
                        (*M)(i+1,j+1) = w.GetArray().At(i,j).GetFloat();
               }
               ret = 0;
               break;
            }

         default:
            ret = 1;
            break;
      }
   }
   return ret;
}

} /* namespace RITA */
