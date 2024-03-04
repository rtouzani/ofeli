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

                            Definition of class 'calc'

  ==============================================================================*/

#pragma once

#include "../muparserx/mpParser.h"
#include "../muparserx/mpDefines.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Matrix.h"

namespace RITA {

using namespace mup;

class rita;
class cmd;
class data;

class FctMatrix : public ICallback
{
 public:
   FctMatrix() : ICallback(cmFUNC, _T("matrix"), -1) {}
   virtual ~FctMatrix() {}
   virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int argc) override;
   virtual const char_type* GetDesc() const override { return _T("matrix(x [, y]) - Returns a matrix whose elements are all 0."); }
   virtual IToken* Clone() const override { return new FctMatrix(*this); }
};


class FctMatrixNorm2 : public ICallback
{
 public:
    FctMatrixNorm2() : ICallback(cmFUNC, _T("Norm2"), 1) {}
    ~FctMatrixNorm2() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("Norm2(M) - Returns 2-norm (euclidean) of matrix."); }
    virtual IToken* Clone() const { return new FctMatrixNorm2(*this); }
};


class FctMatrixNormI : public ICallback
{
 public:
    FctMatrixNormI() : ICallback(cmFUNC, _T("NormI"), 1) {}
    ~FctMatrixNormI() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("normI(M) - Returns infinity-norm of matrix."); }
    virtual IToken* Clone() const { return new FctMatrixNormI(*this); }
};


class FctVector : public ICallback
{
 public:
    FctVector() : ICallback(cmFUNC, _T("vector"), -1) {}
    virtual ~FctVector() {}
    virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int argc) override;
    virtual const char_type* GetDesc() const override { return _T("vector(x) - Returns a vector whose elements are all 0."); }
    virtual IToken* Clone() const override { return new FctVector(*this); }
};


class FctVectorNorm1 : public ICallback
{
 public:
    FctVectorNorm1() : ICallback(cmFUNC, _T("norm1"), 1) {}
    ~FctVectorNorm1() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("norm1(v) - Returns 1-norm of vector."); }
    virtual IToken* Clone() const { return new FctVectorNorm1(*this); }
};


class FctVectorNorm2 : public ICallback
{
 public:
    FctVectorNorm2() : ICallback(cmFUNC, _T("norm2"), 1) {}
    ~FctVectorNorm2() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("norm2(v) - Returns 2-norm (euclidean) of vector."); }
    virtual IToken* Clone() const { return new FctVectorNorm2(*this); }
};


class FctVectorNormI : public ICallback
{
 public:
    FctVectorNormI() : ICallback(cmFUNC, _T("normI"), 1) {}
    ~FctVectorNormI() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("normI(v) - Returns infinity-norm of vector."); }
    virtual IToken* Clone() const { return new FctVectorNormI(*this); }
};


class FctVectorCanonical : public ICallback
{
 public:
    FctVectorCanonical() : ICallback(cmFUNC, _T("canonical"), -1) {}
    ~FctVectorCanonical() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("canonical(n,i) - Returns i-th vector of canonical basis of Rn."); }
    virtual IToken* Clone() const { return new FctVectorCanonical(*this); }
};


class FctMatrixDiag : public ICallback
{
  public:
    FctMatrixDiag() : ICallback(cmFUNC, _T("diag"), -1) {}
    ~FctMatrixDiag() {}
    virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int a_iArgc) override;
    virtual const char_type* GetDesc() const override { return _T("diag(n,d) - returns diagonal matrix with entries d"); }
    virtual IToken* Clone() const override { return new FctMatrixDiag(*this); }
};


class FctMatrixLaplace1D : public ICallback
{
  public:
    FctMatrixLaplace1D() : ICallback(cmFUNC, _T("laplace1d"), -1) {}
    ~FctMatrixLaplace1D() {}
    virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int a_iArgc) override;
    virtual const char_type* GetDesc() const override { return _T("laplace1d(n,c) - returns matrix of discretization"
                                                                  " of -u'' with 3-stencil finite difference scheme"); }
    virtual IToken* Clone() const override { return new FctMatrixLaplace1D(*this); }
};


class calc {

 public:
    calc(rita* r, cmd* command);
    ~calc() { }
    int run();
    int run(const string& buffer);
    int CheckKeywords();
    int getVar(string_type& s);
    int setParam(const string& name, double v);
    int setVector(OFELI::Vect<double> *u);
    int setMatrix(OFELI::Matrix<double> *M);
    int setVectorValue(const string& name, const OFELI::Vect<double>* u);
    int setMatrixValue(const string& name, OFELI::Matrix<double>* M);
    int get_int(const string& s);
    double get_double(const string& s);
    ParserX theParser;

 private:
    rita *_rita;
    data *_data;
    string _sLine;
    cmd *_cmd;
    var_maptype _vmap;
    vector<Value *> _theV;
    Value *_v;

    void ListVar();
    void ListConst();
    void ListExprVar();
    void setData();
    int parse();
    void addVar();

};

} /* namespace RITA */