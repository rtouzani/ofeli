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

                     Implementation of class 'approximation'

  ==============================================================================*/

#include "approximation.h"
#include "data.h"
#include "calc.h"
#include "defs.h"
#include "helps.h"

namespace RITA {

approximation::approximation(rita*      r,
                             cmd*       command,
                             configure* config)
              : _rita(r), _configure(config), _cmd(command)
{
   _data = _rita->_data;
   _theTab = nullptr;
   _tab_alloc = false;
   _verb = _configure->getVerbose();
}


approximation::~approximation()
{
   if (_tab_alloc)
      delete _theTab;
}


int approximation::run()
{
   _pr = _PR;
   _rita->_analysis_type = analysis_type::APPROXIMATION;
   string file="", name="", tab="", ff="";
   vector<funct> fct;
   _method = NONE;
   int ret=0, nb=0, key=0;
   size_t nb_args = _cmd->getNbArgs();
   CHK_MSGR(nb_args>1,_pr,"Illegal number of arguments")

   if (nb_args==1) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available sub-commands: " << Approx_help << endl;
            return 0;

         case 101:
            cout << Approx_Help << endl;
            return 0;

         case 104:
         case 105:
            *_rita->ofh << " end" << endl;
            return 0;

         default:
            UNKNOWN_ARG(_pr)
      }
   }

   static const vector<string> kw {"lagrange","fit$ting","bspline","bezier","nurbs"};
   *_rita->ofh << "approximation" << endl;
   while (1) {
      READLINE
      switch (key) {

         case 0:
            _method = LAGRANGE;
            setLagrange();
            break;

         case 1:
            _method = FITTING;
            setFitting();
            break;

         case 2:
            _method = BSPLINE;
            setBSpline();
            break;

         case 3:
            _method = BEZIER;
            setBezier();
            break;

         case 4:
            _method = NURBS;
            setNurbs();
            break;

         case 100:
             cout << "Available commands: " << Approx_help << endl;
             break;

         case 101:
            cout << Approx_Help << endl;
            break;

         case 102:
            _rita->getLicense();
            break;

         case 104:
         case 105:
            *_rita->ofh << "  end" << endl;
            return ret;

         case 106:
            _rita->setEcho(nb_args);
            break;

         default:
            DEFAULT_KW(_rita)
      }
   }
   return 0;
}


void approximation::setXY()
{
   size_t np = _theTab->getSize(1,1);
   _x.resize(np);
   double h = (_theTab->getMaxVar(1,1)-_theTab->getMinVar(1,1))/(np-1);
   _x[0] = _theTab->getMinVar(1,1);
   (*_y)[0] = _theTab->Funct[0].Val(1);
   for (size_t i=1; i<np; ++i) {
      _x[i] = _x[i-1] + h;
      (*_y)[i] = _theTab->Funct[0].Val(i+1);
   }
}


int approximation::setLagrange()
{
   string pr=_pr+"lagrange>";
   int nb=0, k=0;
   _iTab = 0;
   _file_count = _tab_count = 0;
   _tab = _file = _fft = "";
   static const vector<string> kw {"file","tab$ulation","func$tion"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(pr)
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Approx_Lagrange_help << endl;
            return 0;

         case 101:
            cout << Approx_Lagrange_Help << endl;
            return 0;

         case   0:
            NO_VALUE_ARG(pr)
            _file = _cmd->string_token(0);
            _file_count++;
            break;

         case   1:
            NO_VALUE_ARG(pr)
            _tab = _cmd->string_token(0);
            _tab_count++;
            break;

         case   2:
            NO_VALUE_ARG(pr)
            _fft = _cmd->string_token(0);
            break;

         default:
            MSGR(pr,"Unknown argument: "+_cmd->getToken())
      }
   }

   CHK_MSGR(_file_count==0 && _tab_count==0,pr,"A data file or a tabulation name must be given.")
   CHK_MSGR(_file_count>0 && _tab_count>0,pr,"File and tabulation cannot be given simultaneously.")
   *_rita->ofh << "  lagrange";
   if (_tab_count>0 && _tab!="") {
      k = _data->TabLabel[_tab];
      CHK_MSGR(k<=0,pr,"Reference to undefined tabulation"+_tab)
      _theTab = _data->theTab[k];
      *_rita->ofh << " tabulation=" << _tab;
   }
   if (_file_count && _file!="") {
      _theTab = new OFELI::Tabulation(_file);
      _tab_alloc = true;
      _iTab = _data->addTab(_theTab,_tab);
      *_rita->ofh << " file=" << _file;
   }
   if (_fft!="")
      *_rita->ofh << " function=" << _fft;
   *_rita->ofh << endl;
   CHK_MSGR(_theTab->getNbVar(1)>1,pr,"This approximation method is available for one-variable cases only.")
   _y = new Vect<double>(_theTab->getSize(1,1));
   setXY();
   _degree = _x.size() - 1;
   k = _data->addFunction(_fft);
   CHK_MSGR(k<=0,pr,".")
   _data->theFct[k]->setVar("x",1);
   _data->theFct[k]->setNoPar();
   _ff = _data->theFct[k]->getFct();
   OFELI::FuncApprox fa;
   fa.setLagrange(_degree,_x,*_y,*_ff);
   fa.run();
   _rita->_calc->setFun(_data->theFct[k]);
   if (_verb)
      cout << "Lagrange interpolation created. Lagrange polynomial is: " << _ff->getName() 
           << "(" << _ff->getVar(1) << ")" << endl;
   if (_tab_count==0)
      _data->remove(_data->TabName[_iTab]);
   return 0;
}


int approximation::setFitting()
{
   string pr=_pr+"fitting>";
   static const vector<string> kw {"file","tab$ulation","func$tion","degree","basis"};
   _fct_names.clear();
   _tab = _file = _fft = "";
   _basis_count = _degree = 0;
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(pr)
   int nb=0, k=0;
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Approx_Fitting_help << endl;
            return 0;

         case 101:
            cout << Approx_Fitting_Help << endl;
            return 0;

         case   0:
            NO_VALUE_ARG(pr)
            _file = _cmd->string_token(0);
            _file_count++;
            break;

         case   1:
            NO_VALUE_ARG(pr)
            _tab = _cmd->string_token(0);
            _tab_count++;
            break;

         case   2:
            NO_VALUE_ARG(pr)
            _fft = _cmd->string_token(0);
            break;

         case   3:
            NO_VALUE_ARG(pr)
            _degree = _cmd->int_token(0);
            break;

         case   4:
            NO_VALUE_ARG(pr)
            for (int i=0; i<nb; ++i) {
               _fct_names.push_back(_cmd->string_token(i));
               _basis_count++;
            }
            break;

         default:
            MSGR(pr,"Unknown argument: "+_cmd->getToken())
      }
   }

   CHK_MSGR(_file_count==0 && _tab_count==0,pr,"A data file or a tabulation name must be given.")
   CHK_MSGR(_file_count>0 && _tab_count>0,pr,"File and tabulation cannot be given simultaneously.")
   CHK_MSGR(_degree>0 && _basis_count>0,pr,"Polynomial degree and basis functions cannot be given simultaneously.")
   if (_degree==0 && _basis_count==0)
      _degree = 1;
   *_rita->ofh << "  fitting";
   if (_tab_count>0 && _tab!="") {
      k = _data->TabLabel[_tab];
      CHK_MSGR(k<=0,pr,"Reference to undefined tabulation"+_tab)
      _theTab = _data->theTab[k];
      *_rita->ofh << " tabulation=" << _tab;
   }
   if (_file_count && _file!="") {
      _theTab = new OFELI::Tabulation(_file);
      _tab_alloc = true;
      _iTab = _data->addTab(_theTab,_tab);
      *_rita->ofh << " file=" << _file;
   }
   CHK_MSGR(_theTab->getNbVar(1)>1,pr,"This approximation method is available for one-variable cases only.")
   if (_degree)
      *_rita->ofh << " degree=" << _degree;
   _nb_fit = _basis_count;
   if (_nb_fit) {
      *_rita->ofh << " basis=" << _fct_names[0];
      for (int i=1; i<_nb_fit; ++i)
         *_rita->ofh << "," << _fct_names[i];
   }
   if (_fft!="")
      *_rita->ofh << " function=" << _fft;
   *_rita->ofh << endl;
   int np = _theTab->getSize(1,1);
   CHK_MSGR(_degree>0 && np<_degree+1,pr,"Degree of fitting polynomial is too large.")
   _y = new Vect<double>(np);
   setXY();
   FuncApprox fa;
   Vect<double> a;

// Case of polynomial fitting
   if (_degree) {
      a.resize(_degree+1);
      fa.setLeastSquare(_x,*_y,_degree,a);
      fa.run();
   }

// Case of arbitrary basis fitting
   if (_nb_fit) {
      a.resize(_nb_fit);
      for (auto const& f: _fct_names) {
         Fct *ff = new Fct(f,"x");
         _fct.push_back(ff);
      }
      fa.setLeastSquare(_fct,_x,*_y,a);
      fa.run();
   }
   k = _data->addFunction(_fft);
   CHK_MSGR(k<=0,pr,".")
   fa.getLeastSquare(*_data->theFct[k]->getFct());
   _data->theFct[k]->setFromFct();
   _rita->_calc->setFun(_data->theFct[k]);
   if (_verb)
      cout << "Curve fitting created.\nFitting function: " << _fft << endl;
   if (_tab_count==0)
      _data->remove(_data->TabName[_iTab]);
   for (auto const& f: _fct)
      delete f;
   return 0;
}


int approximation::setBSpline()
{
   string pr=_pr+"bspline>";
   static const vector<string> kw {"file","vert$ices","order","nc","np","vect$or"};
   vector<double> x;
   _file = "";
   string vname = "";
   _degree = 2;
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   int npu=0, npv=1, ncu=0, ncv=1;
   NO_ARG(pr)
   int nb=0, k=0;
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Approx_BSpline_help << endl;
            return 0;

         case 101:
            cout << Approx_BSpline_Help << endl;
            return 0;

         case   0:
            NO_VALUE_ARG(pr)
            _file = _cmd->string_token(0);
            _file_count++;
            break;

         case   1:
            NO_VALUE_ARG(pr)
            k = 0;
            for (int i=0; i<nb; ++i)
               x.push_back(_cmd->double_token(k++));
            break;

         case   2:
            NO_VALUE_ARG(pr)
            _degree = _cmd->int_token(0);
            break;

         case   3:
            NO_VALUE_ARG(pr)
            CHK_MSGR(nb>2,pr,"Two values of nc can be given at most.")
            ncu = _cmd->int_token(0);
            if (nb>1)
               ncv = _cmd->int_token(1);
            break;

         case   4:
            NO_VALUE_ARG(pr)
            CHK_MSGR(nb>2,pr,"Two values of np can be given at most.")
            npu = _cmd->int_token(0);
            if (nb>1)
               npv = _cmd->int_token(1);
            break;

         case   5:
            NO_VALUE_ARG(pr)
            vname = _cmd->string_token(0);
            break;

         default:
            MSGR(pr,"Unknown argument: "+_cmd->getToken())
      }
   }

   CHK_MSGR(_file_count==0 && ncu==0,pr,"A data file or a list of vertices must be given.")
   CHK_MSGR(npu==0,pr,"Number of points on resulting curve or surface must be given.")
   *_rita->ofh << "  bspline";
   if (_file_count && _file!="")
      *_rita->ofh << " file=" << _file;
   *_rita->ofh << " order=" << _degree;
   if (ncu>0) {
      *_rita->ofh << " nc=" << ncu;
      if (ncv>1)
         *_rita->ofh << "," << ncv;
      _x.setSize(ncu,ncv,3);
      k = 0;
      for (int i=1; i<=ncu; ++i) {
         for (int j=1; j<=ncv; ++j) {
            _x(i,j,1) = x[k++];
            _x(i,j,2) = x[k++];
            _x(i,j,3) = 0.;
            if (ncv>1)
               _x(i,j,3) = x[k++];
         }
      }
   }
   else {
      XMLParser xml(_file,EType::VECTOR);
      xml.get(_x);
      ncu = _x.getNx(), ncv = _x.getNy();
   }
   *_rita->ofh << " np=" << npu;
   if (npv>1)
      *_rita->ofh << "," << npv;
   if (vname!="")
      *_rita->ofh << " vector=" << vname;
   if (x.size()) {
      *_rita->ofh << " vertices=";
      for (size_t i=0; i<x.size()-1; ++i)
         *_rita->ofh << x[i] << ",";
      *_rita->ofh << x[x.size()-1];
   }
   *_rita->ofh << endl;
   _y = new Vect<double>(npu,npv,3);
   k = _data->addVector(_y,vname,SetCalc::SET);
   FuncApprox fa;
   if (ncv==1)
      fa.setBSpline(ncu,_degree,npu,_x,*_y);
   else
      fa.setBSplineSurface(ncu,ncv,_degree,_degree,npu,npv,_x,*_y);
   fa.run();
   if (_verb)
      cout << "Resulting curve is stored in vector: " << _data->VectorName[k] << endl;
   return 0;
}


int approximation::setBezier()
{
   string pr=_pr+"bezier>";
   static const vector<string> kw {"file","vert$ices","nc","np","vect$or"};
   vector<double> x;
   _file = "";
   string vname = "";
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   int npu=0, npv=1, ncu=0, ncv=1;
   NO_ARG(pr)
   int nb=0, k=0;
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Approx_Bezier_help << endl;
            return 0;

         case 101:
            cout << Approx_Bezier_Help << endl;
            return 0;

         case   0:
            NO_VALUE_ARG(pr)
            _file = _cmd->string_token(0);
            _file_count++;
            break;

         case   1:
            NO_VALUE_ARG(pr)
            k = 0;
            for (int i=0; i<nb; ++i)
               x.push_back(_cmd->double_token(k++));
            break;

         case   2:
            NO_VALUE_ARG(pr)
            CHK_MSGR(nb>2,pr,"Two values of nc can be given at most.")
            ncu = _cmd->int_token(0);
            if (nb>1)
               ncv = _cmd->int_token(1);
            break;

         case   3:
            NO_VALUE_ARG(pr)
            CHK_MSGR(nb>2,pr,"Two values of np can be given at most.")
            npu = _cmd->int_token(0);
            if (nb>1)
               npv = _cmd->int_token(1);
            break;

         case   4:
            NO_VALUE_ARG(pr)
            vname = _cmd->string_token(0);
            break;

         default:
            MSGR(pr,"Unknown argument: "+_cmd->getToken())
      }
   }

   CHK_MSGR(_file_count==0 && ncu==0,pr,"A data file or a list of vertices must be given.")
   CHK_MSGR(npu==0,pr,"Number of points on resulting curve or surface must be given.")
   *_rita->ofh << "  bezier";
   if (_file_count && _file!="")
      *_rita->ofh << " file=" << _file;
   if (ncu>0) {
      *_rita->ofh << " nc=" << ncu;
      if (ncv>1)
         *_rita->ofh << "," << ncv;
      _x.setSize(ncu,ncv,3);
      k = 0;
      for (int i=1; i<=ncu; ++i) {
         for (int j=1; j<=ncv; ++j) {
            _x(i,j,1) = x[k++];
            _x(i,j,2) = x[k++];
            _x(i,j,3) = 0.;
            if (ncv>1)
               _x(i,j,3) = x[k++];
         }
      }
   }
   else {
      XMLParser xml(_file,EType::VECTOR);
      xml.get(_x,"");
      ncu = _x.getNx(), ncv = _x.getNy();
   }
   *_rita->ofh << " np=" << npu;
   if (npv>1)
      *_rita->ofh << "," << npv;
   if (vname!="")
      *_rita->ofh << " vector=" << vname;
   *_rita->ofh << endl;
   _y = new Vect<double>(npu,npv,3);
   k = _data->addVector(_y,vname,SetCalc::SET);
   FuncApprox fa;
   if (ncv==1)
      fa.setBezier(ncu,npu,_x,*_y);
   else
      fa.setBezierSurface(ncu,ncv,npu,npv,_x,*_y);
   fa.run();
   if (_verb)
      cout << "Resulting curve is stored in vector: " << _data->VectorName[k] << endl;
   return 0;
}


int approximation::setNurbs()
{
   string pr=_pr+"nurbs>";
   static const vector<string> kw {"file","vert$ices","nc","np","order","vect$or"};
   vector<double> x;
   _file = "";
   string vname = "";
   _degree = 2;
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   int npu=0, npv=1, ncu=0, ncv=1;
   NO_ARG(pr)
   int nb=0, k=0;
   for (size_t i=0; i<nb_args; ++i) {
      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Approx_Nurbs_help << endl;
            return 0;

         case 101:
            cout << Approx_Nurbs_Help << endl;
            return 0;

         case   0:
            NO_VALUE_ARG(pr)
            _file = _cmd->string_token(0);
            _file_count++;
            break;

         case   1:
            NO_VALUE_ARG(pr)
            k = 0;
            for (int i=0; i<nb; ++i)
               x.push_back(_cmd->double_token(k++));
            break;

         case   2:
            NO_VALUE_ARG(pr)
            CHK_MSGR(nb>2,pr,"Two values of nc can be given at most.")
            ncu = _cmd->int_token(0);
            if (nb>1)
               ncv = _cmd->int_token(1);
            break;

         case   3:
            NO_VALUE_ARG(pr)
            CHK_MSGR(nb>2,pr,"Two values of np can be given at most.")
            npu = _cmd->int_token(0);
            if (nb>1)
               npv = _cmd->int_token(1);
            break;

         case   4:
            NO_VALUE_ARG(pr)
            _degree = _cmd->int_token(0);

         case   5:
            NO_VALUE_ARG(pr)
            vname = _cmd->string_token(0);
            break;

         default:
            MSGR(pr,"Unknown argument: "+_cmd->getToken())
      }
   }

   CHK_MSGR(_file_count==0 && ncu==0,pr,"A data file or a list of vertices must be given.")
   CHK_MSGR(npu==0,pr,"Number of points on resulting curve or surface must be given.")
   *_rita->ofh << "  bezier";
   if (_file_count && _file!="")
      *_rita->ofh << " file=" << _file;
   if (ncu>0) {
      *_rita->ofh << " nc=" << ncu;
      if (ncv>1)
         *_rita->ofh << "," << ncv;
      _x.setSize(ncu,ncv,3);
      k = 0;
      for (int i=1; i<=ncu; ++i) {
         for (int j=1; j<=ncv; ++j) {
            _x(i,j,1) = x[k++];
            _x(i,j,2) = x[k++];
            _x(i,j,3) = 0.;
            if (ncv>1)
               _x(i,j,3) = x[k++];
         }
      }
   }
   else {
      XMLParser xml(_file,EType::VECTOR);
      xml.get(_x,"");
      ncu = _x.getNx(), ncv = _x.getNy();
   }
   *_rita->ofh << " order=" << _degree;
   *_rita->ofh << " np=" << npu;
   if (npv>1)
      *_rita->ofh << "," << npv;
   if (vname!="")
      *_rita->ofh << " vector=" << vname;
   *_rita->ofh << endl;
   _y = new Vect<double>(npu,npv,3);
   k = _data->addVector(_y,vname,SetCalc::SET);
   Vect<double> h(ncu);
   h = 1.;
   FuncApprox fa;
   if (ncv==1)
      fa.setNurbs(ncu,_degree,npu,h,_x,*_y);
   else
      fa.setNurbsSurface(ncu,ncv,_degree,_degree,npu,npv,_x,*_y);
   fa.run();
   if (_verb)
      cout << "Resulting curve is stored in vector: " << _data->VectorName[k] << endl;
   return 0;
}

} /* namespace RITA */
