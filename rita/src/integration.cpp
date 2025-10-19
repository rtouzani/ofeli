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

                      Implementation of class 'integration'

  ==============================================================================*/


#include "integration.h"
#include "data.h"
#include "defs.h"
#include "helps.h"

namespace RITA {

integration::integration(rita*      r,
                         cmd*       command,
                         configure* config)
            : _rita(r), _configure(config), _cmd(command)
{
}


integration::~integration()
{
}


int integration::run()
{
   int ret1=0, ret2=0, ret3=0, ret=0;
   _rita->_analysis_type = analysis_type::INTEGRATION;
   string fct="", def="", var_name="", form="trapezoidal", name="";
   int count_fct=0, count_def=0, count_vector=0;
   int nb=0;
   result="";
   nx = ny = nz = 10;
   xmin = ymin = zmin = 0.;
   xmax = ymax = zmax = 1.;
   dim = 1;
   unif = 1;
   ng = 1;
   var.clear();
   _data = _rita->_data;
   static const vector<string> kw {"func$tion","def$inition","x","var$iable","vect$or","ne","form$ula","res$ult"};
   _cmd->set(kw,_rita->_gkw);
   size_t nb_args = _cmd->getNbArgs();
   NO_ARG(_pr)
   for (size_t k=0; k<nb_args; ++k) {

      switch (_cmd->getArgs(nb)) {

         case 100:
            cout << "Available arguments: " << Int_help << endl;
            return 1;

         case 101:
            cout << Int_Help << endl;
            return 1;

         case   0:
            NO_VALUE_ARG(_pr)
            fct = _cmd->string_token(0);
            count_fct++;
            break;

         case   1:
            NO_VALUE_ARG(_pr)
            def = _cmd->string_token(0);
            count_def++;
            break;

         case   2:
            NO_VALUE_ARG(_pr)
            ret1  = _data->getPar(0,_pr,xmin);
            if (nb==2)
               ret1 += _data->getPar(1,_pr,xmax);
            break;

         case   3:
         case   4:
            NO_VALUE_ARG(_pr)
            var_name = _cmd->string_token(0);
            count_vector++;
            break;

         case   5:
            NO_VALUE_ARG(_pr)
            ret3  = _data->getPar(0,_pr,nx);
            break;

         case   6:
            NO_VALUE_ARG(_pr)
            form = _cmd->string_token(0);
            if (nb>1)
               ret = _data->getPar(1,_pr,ng);
            break;

         case   7:
            result  = _cmd->string_token(0);
            break;

         default:
            UNKNOWN_ARG(_pr)
      }
   }

   if (nb_args>0) {
      CHK_MSGR(count_fct==0 && count_def==0,_pr,"No defined function to integrate.")
      CHK_MSGR(count_fct && count_vector,_pr,"Function already defined.")
      CHK_MSGR(count_fct && count_def,_pr,"Function already defined.")
      CHK_MSGR(count_fct>1 || count_def>1,_pr,"Too many functions defined.")
      CHK_MSGR(count_def && !count_vector,_pr,"Missing variable name.")
      CHK_MSGR(ret1||ret2||ret3,_pr,"Error in data.")
      CHK_MSGR(xmin>=xmax,_pr,"xmin must be smaller than xmax.")
      CHK_MSGR(ret2,_pr,"Error in ne data.")
      *_rita->ofh << "integration";
      if (dim==1)
         *_rita->ofh << " x=" << xmin << "," << xmax << " ne=" << nx;
      if (count_fct) {
         iFct = _data->FctLabel[fct];
         CHK_MSGR(iFct==-1,_pr,"Non defined function "+fct)
         *_rita->ofh << " function=" << fct;
      }
      else {
         *_rita->ofh << " var=" << var_name << " definition=" << def;
         if (dim==1)
            var.push_back(var_name);
         else {
            for (int i=0; i<dim; ++i)
               var.push_back(var_name+to_string(i+1));
         }
         iFct = _data->addFunction(name,var[0],def);
         if (iFct<0)
            return 1;
      }
      nim = Nint[form];
      *_rita->ofh << " formula=" << form;
      if (nim==GAUSS_LEGENDRE || nim==GAUSS_LOBATTO)
         *_rita->ofh << "," << ng;
      if (result!="")
         *_rita->ofh << " result=" << result;
      *_rita->ofh << endl;
   }
   return ret;
}


int integration::go()
{
   res = 0.;
   double dx=0.;
   if (unif==1)
      dx = (xmax-xmin)/nx;
   static const vector<double> xg {0.,-0.5773502691896257,0.5773502691896257,0.,-0.7745966692414834,
                                   0.7745966692414834,-0.3399810435848563,0.3399810435848563,
                                   -0.8611363115940526,0.8611363115940526,0.,-0.5384693101056831,
                                   0.5384693101056831,-0.9061798459386640,0.9061798459386640,
                                   -0.6612093864662645,0.6612093864662645,-0.2386191860831969,
                                   0.2386191860831969,-0.9324695142031521,0.9324695142031521};
   static const vector<double> wg {2.,1.,1.,0.8888888888888889,0.5555555555555556,0.5555555555555556,
                                   0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538,
                                   0.5688888888888889,0.4786286704993665,0.4786286704993665,0.2369268850561891,
                                   0.2369268850561891,0.3607615730481386,0.3607615730481386,0.4679139345726910,
                                   0.4679139345726910,0.1713244923791704,0.1713244923791704};
   static const vector<double> xl {0.,-1.0,1.0,-0.447213595499958,0.447213595499958,-1.0,1.0,0.,-0.475481063909541,
                                   0.475481063909541,-1.0,1.0};
   static const vector<double> wl {1.333333333333333,1.333333333333333,0.333333333333333,0.833333333333333,
                                   0.833333333333333,0.166666666666667,0.166666666666667,0.711111111111111,
                                   0.544444444444444,0.544444444444444,0.1,0.1};

   _x.resize(nx+1);
   _x[0] = xmin;
   for (int i=1; i<=nx; ++i)
      _x[i] = _x[i-1] + dx;
   funct &f = *_data->theFct[iFct];
   for (int ii=0; ii<nx; ++ii) {
      double x1=_x[ii], x2=_x[ii+1], x12=0.5*(_x[ii]+_x[ii+1]); 
      if (nim==LRECTANGLE)
         res += dx*f(x1);
      else if (nim==RRECTANGLE)
         res += dx*f(x2);
      else if (nim==MIDPOINT)
         res += dx*f(x12);
      else if (nim==TRAPEZOIDAL)
         res += 0.5*dx*(f(x1)+f(x2));
      else if (nim==SIMPSON)
         res += (f(x1)+4*f(x12)+f(x2))*dx/6.0;
      else if (nim==GAUSS_LEGENDRE) {
         CHK_MSGR(ng<1 || ng>5,_pr,"For Gauss-Legendre formula, number of points must be between 1 and 6.")
         for (int i=0; i<ng; ++i) {
            int j = ng*(ng-1)/2 + i;
            res += 0.5*dx*wg[j]*f(x12+0.5*dx*xg[j]);
         }
      }
      else if (nim==GAUSS_LOBATTO) {
         CHK_MSGR(ng<3 || ng>5,_pr,"For Gauss-Lobatto formula, number of points must be between 3 and 5.")
         for (int i=0; i<ng; ++i) {
            int j = ng*(ng-1)/2 + i - 3;
            res += 0.5*dx*wl[j]*f(x12+0.5*dx*xl[j]);
         }
      }
   }
   cout << "Approximate Integral: " << res << endl;
   if (result=="")
      result = "I" + to_string(_data->iParam+1);
   _data->addParam(result,res,SetCalc::SET);
   return 0;
}

} /* namespace RITA */
