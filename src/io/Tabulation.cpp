/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                       Implementation of class 'Tabulation'

  ==============================================================================*/


#include "io/Tabulation.h"
#include "io/XMLParser.h"

namespace OFELI {

Tabulation::Tabulation()
{
   _nb_funct = 0;
   _funct.resize(10);
}


Tabulation::Tabulation(string file)
{
   _nb_funct = 0;
   _funct.resize(10);
   setFile(file);
}


Tabulation::~Tabulation() { }


void Tabulation::setFile(string file)
{
   XMLParser p(file,XMLParser::FUNCTION);
   p.get(*this);
}


void Tabulation::setFunction(string label)
{
   _funct[_nb_funct].Label = label;
   _funct[_nb_funct].NbVar = 0;
   _funct_id[label] = ++_nb_funct;
}


void Tabulation::setVariable(string label)
{
   fct &f = _funct[_nb_funct-1];
   size_t n = f.NbVar;
   f.VarLabel[n] = label;
   f.Id[label] = n + 1;
   f.NbVar++;
}


void Tabulation::setSizes()
{
   fct &f = _funct[_nb_funct-1];
   if (f.NbVar==1)
      f.Val.setSize(f.Np[0]);
   else if (f.NbVar==2)
      f.Val.setSize(f.Np[0],f.Np[1]);
   else if (f.NbVar==3)
      f.Val.setSize(f.Np[0],f.Np[1],f.Np[2]);
}


real_t Tabulation::getValue(string funct,
                            real_t v)
{
   fct &f = _funct[_funct_id[funct]-1];
   real_t h = (f.Max[0]-f.Min[0])/(f.Np[0]-1);
   size_t i = size_t((v-f.Min[0])/h);
   real_t s = (v-f.Min[0])/h - i;
   if (v<f.Min[0])
      i=0, s=0;
   else if (v>f.Max[0])
      i=f.Np[0]-1, s=1;
   return (1-s)*f.Val(i+1) + s*f.Val(i+2);
}


real_t Tabulation::getValue(string funct,
                            real_t v1,
                            real_t v2)
{
   fct &f = _funct[_funct_id[funct]-1];
   size_t N1=f.Np[0], N2=f.Np[1];
   real_t h1 = (f.Max[0]-f.Min[0])/(N1-1),
          h2 = (f.Max[1]-f.Min[1])/(N2-1);
   size_t i=size_t((v1-f.Min[0])/h1), j=size_t((v2-f.Min[1])/h2);
   real_t s = (v1-f.Min[0])/h1 - i,
          t = (v2-f.Min[1])/h2 - j;
   if (v1<f.Min[0])
      i=0, s=0;
   if (v1>f.Max[0])
      i=N1-1, s=1;
   if (v2<f.Min[1])
      j=0, t=0;
   if (v2>f.Max[1])
      j=N2-1, t=1;
   return   (1-s)*(1-t)*f.Val(i+1,j+1) + s*(1-t)*f.Val(i+2,j+1)
          + s*t*f.Val(i+2,j+2) + (1-s)*t*f.Val(i+1,j+2);
}


real_t Tabulation::getValue(string funct,
                            real_t v1,
                            real_t v2,
                            real_t v3)
{
   fct &f = _funct[_funct_id[funct]-1];
   size_t N1=f.Np[0], N2=f.Np[1], N3=f.Np[2];
   real_t h1 = (f.Max[0]-f.Min[0])/(N1-1),
          h2 = (f.Max[1]-f.Min[1])/(N2-1),
          h3 = (f.Max[2]-f.Min[2])/(N3-1);
   size_t i=size_t((v1-f.Min[0])/h1),
          j=size_t((v2-f.Min[1])/h2),
          k=size_t((v3-f.Min[2])/h3);
   real_t r = (v1-f.Min[0])/h1 - i,
          s = (v2-f.Min[1])/h2 - j,
          t = (v3-f.Min[2])/h3 - k;
   if (v1<f.Min[0])
      i=0, r=0;
   if (v1>f.Max[0])
      i=N1-1, r=1;
   if (v2<f.Min[1])
      j=0, s=0;
   if (v2>f.Max[1])
      j=N2-1, s=1;
   if (v3<f.Min[2])
      k=0, t=0;
   if (v3>f.Max[2])
      k=N3-1, t=1;
   real_t w1 = (1-r)*(1-s)*(1-t)*f.Val(i+1,j+1,k+1)
             + (1-r)*s*(1-t)*f.Val(i+1,j+2,k+1)
             + (1-r)*(1-s)*t*f.Val(i+1,j+1,k+2);
   real_t w2 = (1-r)*s*t*f.Val(i+1,j+2,k+2)
             + r*(1-s)*(1-t)*f.Val(i+2,j+1,k+1)
             + r*s*(1-t)*f.Val(i+2,j+2,k+1)
             + r*s*t*f.Val(i+2,j+2,k+2);
   return w1+w2;
}


real_t Tabulation::getDerivative(string funct,
                                 real_t v)
{
   fct &f = _funct[_funct_id[funct]-1];
   size_t N = f.Np[0];
   real_t h = (f.Max[0]-f.Min[0])/(N-1);
   size_t i = size_t((v-f.Min[0])/h);
   if (v<f.Min[0])
      i=0;
   else if (v>f.Max[0])
      i=N-1;
   return f.Val(i+2) - f.Val(i+1);
}


ostream& operator<<(ostream&          s,
                    const Tabulation& p)
{
   size_t n = p._nb_funct;
   s << "\n===============================================================================\n";
   s << "Summary of Tabulated function data:\n";
   s << "-------------------------------------------------------------------------------\n";
   for (size_t i=1; i<=n; i++) {
      const fct &f = p._funct[i-1];
      s << " Function: " << f.Label << ", Number of variables: " << f.NbVar << endl;
      for (size_t j=0; j<f.NbVar; j++) {
         s << "    Variable: " << f.VarLabel[j] << ", Min. = " << f.Min[j];
         s << ", Max = " << f.Max[j] << ", Number of points = " << f.Np[j] << endl;
      }
   s << "-------------------------------------------------------------------------------" << endl;
   }
   return s;
}

} /* namespace OFELI */
