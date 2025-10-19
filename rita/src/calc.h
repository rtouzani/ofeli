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

                            Definition of class 'calc'
                            and other related classes

  ==============================================================================*/

#pragma once

#include "../muparserx/mpParser.h"
#include "../muparserx/mpDefines.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Matrix.h"
#include "funct.h"

namespace RITA {

using namespace mup;

class rita;
class cmd;
class data;

class FctFun : public ICallback
{
 public:

    FctFun(funct *f) : ICallback(cmFUNC, _T(f->name.c_str()),-1), _fun(f) { }
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int argc) override;
    virtual const char_type* GetDesc() const override { return "Returns value of function at given point"; }
    virtual IToken* Clone() const override { return new FctFun(*this); }

 private:
    funct *_fun;
    vector<double> _arg;
};


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


class FctMatrixNormMax : public ICallback
{
 public:
    FctMatrixNormMax() : ICallback(cmFUNC, _T("NormMax"), 1) {}
    ~FctMatrixNormMax() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("normMax(M) - Returns infinity-norm of matrix."); }
    virtual IToken* Clone() const { return new FctMatrixNormMax(*this); }
};


class FctRowVector : public ICallback
{
 public:
    FctRowVector() : ICallback(cmFUNC, _T("rvector"), -1) {}
    virtual ~FctRowVector() {}
    virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int argc) override;
    virtual const char_type* GetDesc() const override { return _T("rvector(x) - Returns a row vector whose elements are all 0."); }
    virtual IToken* Clone() const override { return new FctRowVector(*this); }
};


class FctRowVectorNorm1 : public ICallback
{
 public:
    FctRowVectorNorm1() : ICallback(cmFUNC, _T("rnorm1"), 1) {}
    ~FctRowVectorNorm1() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("rnorm1(v) - Returns 1-norm of row vector."); }
    virtual IToken* Clone() const { return new FctRowVectorNorm1(*this); }
};


class FctRowVectorNorm2 : public ICallback
{
 public:
    FctRowVectorNorm2() : ICallback(cmFUNC, _T("rnorm2"), 1) {}
    ~FctRowVectorNorm2() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("rnorm2(v) - Returns 2-norm (euclidean) of row vector."); }
    virtual IToken* Clone() const { return new FctRowVectorNorm2(*this); }
};


class FctRowVectorNormMax : public ICallback
{
 public:
    FctRowVectorNormMax() : ICallback(cmFUNC, _T("rnormMax"), 1) {}
    ~FctRowVectorNormMax() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("rnormMax(v) - Returns infinity-norm of row vector."); }
    virtual IToken* Clone() const { return new FctRowVectorNormMax(*this); }
};


class FctRowVectorCanonical : public ICallback
{
 public:
    FctRowVectorCanonical() : ICallback(cmFUNC, _T("rcanonical"), -1) {}
    ~FctRowVectorCanonical() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("rcanonical(n,i) - Returns i-th vector of canonical basis of Rn."); }
    virtual IToken* Clone() const { return new FctRowVectorCanonical(*this); }
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


class FctVectorNormMax : public ICallback
{
 public:
    FctVectorNormMax() : ICallback(cmFUNC, _T("normMax"), 1) {}
    ~FctVectorNormMax() {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const { return _T("normMax(v) - Returns infinity-norm of vector."); }
    virtual IToken* Clone() const { return new FctVectorNormMax(*this); }
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
    calc(data* d, cmd* command, ofstream *ofh, ofstream *ofl);
    ~calc() { delete theParser; }
    int run();
    int setParam(const string& name, double v);
    int setFun(funct *f);
    int setFun(funct *f, const string& exp, const vector<string>& var);
    int setFun(funct *f, const string& exp, const vector<string>& var, const vector<size_t>& n);
    int setVector(OFELI::Vect<double> *u);
    int setMatrix(OFELI::Matrix<double> *M);
    int setVectorValue(const string& name, const OFELI::Vect<double>* u);
    int setMatrixValue(const string& name, OFELI::Matrix<double>* M);
    void ListFun();
    void ListVar();
    void ListConst();
    void ListExprVar();
    int get_int(const string& s);
    double get_double(const string& s);
    ParserX *theParser;

 private:
    data *_data;
    cmd *_cmd;
    ofstream *_ofh, *_ofl;
    int _nb_def_funs;
    vector<string> _funs, _params, _vects, _matrices;
    int CheckKeywords(const string& sLine);
};

} /* namespace RITA */