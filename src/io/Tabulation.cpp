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
   Funct.resize(10);
}


Tabulation::Tabulation(string file)
{
   _nb_funct = 0;
   Funct.resize(10);
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
   Funct[_nb_funct].Label = label;
   Funct[_nb_funct].NbVar = 0;
   _funct_id[label] = ++_nb_funct;
}


void Tabulation::setVariable(string label)
{
   fct &f = Funct[_nb_funct-1];
   size_t n = f.NbVar;
   f.VarLabel[n] = label;
   f.Id[label] = n + 1;
   f.NbVar++;
}


void Tabulation::setSizes()
{
   fct &f = Funct[_nb_funct-1];
   if (f.NbVar==1)
      f.Val.setSize(f.Np[0]);
   else if (f.NbVar==2)
      f.Val.setSize(f.Np[0],f.Np[1]);
   else if (f.NbVar==3)
      f.Val.setSize(f.Np[0],f.Np[1],f.Np[2]);
   else if (f.NbVar==4)
      f.Val.setSize(f.Np[0],f.Np[1],f.Np[2],f.Np[3]);
}


real_t Tabulation::getValue(string funct,
                            real_t x)
{
   fct &f = Funct[_funct_id[funct]-1];
   real_t h = (f.Max[0]-f.Min[0])/(f.Np[0]-1);
   size_t i = size_t((x-f.Min[0])/h);
   real_t xs = (x-f.Min[0])/h - i;
   if (x<f.Min[0])
      i=0, xs=0;
   else if (x>f.Max[0])
      i=f.Np[0]-1, xs=1;
   return (1-xs)*f.Val(i+1) + xs*f.Val(i+2);
}


real_t Tabulation::getValue(string funct,
                            real_t x,
                            real_t y)
{
   fct &f = Funct[_funct_id[funct]-1];
   size_t N1=f.Np[0], N2=f.Np[1];
   real_t h1 = (f.Max[0]-f.Min[0])/(N1-1),
          h2 = (f.Max[1]-f.Min[1])/(N2-1);
   size_t i=size_t((x-f.Min[0])/h1), j=size_t((y-f.Min[1])/h2);
   real_t xs = (x-f.Min[0])/h1 - i,
          ys = (y-f.Min[1])/h2 - j;
   if (x<f.Min[0])
      i=0, xs=0;
   if (x>f.Max[0])
      i=N1-1, xs=1;
   if (y<f.Min[1])
      j=0, ys=0;
   if (y>f.Max[1])
      j=N2-1, ys=1;
   return   (1-xs)*(1-ys)*f.Val(i+1,j+1) + xs*(1-ys)*f.Val(i+2,j+1)
          + xs*ys*f.Val(i+2,j+2) + (1-xs)*ys*f.Val(i+1,j+2);
}


real_t Tabulation::getValue(string funct,
                            real_t x,
                            real_t y,
                            real_t z)
{
   fct &f = Funct[_funct_id[funct]-1];
   size_t N1=f.Np[0], N2=f.Np[1], N3=f.Np[2];
   real_t h1 = (f.Max[0]-f.Min[0])/(N1-1),
          h2 = (f.Max[1]-f.Min[1])/(N2-1),
          h3 = (f.Max[2]-f.Min[2])/(N3-1);
   size_t i=size_t((x-f.Min[0])/h1),
          j=size_t((y-f.Min[1])/h2),
          k=size_t((z-f.Min[2])/h3);
   real_t xs = (x-f.Min[0])/h1 - i,
          ys = (y-f.Min[1])/h2 - j,
          zs = (z-f.Min[2])/h3 - k;
   if (x<f.Min[0])
      i=0, xs=0;
   else if (x>f.Max[0])
      i=N1-1, xs=1;
   if (y<f.Min[1])
      j=0, ys=0;
   else if (y>f.Max[1])
      j=N2-1, ys=1;
   if (z<f.Min[2])
      k=0, zs=0;
   else if (z>f.Max[2])
      k=N3-1, zs=1;
   real_t w1 = (1-xs)*(1-ys)*(1-zs)*f.Val(i+1,j+1,k+1)
             + (1-xs)*ys*(1-zs)*f.Val(i+1,j+2,k+1)
             + (1-xs)*(1-ys)*zs*f.Val(i+1,j+1,k+2)
             + (1-xs)*ys*zs*f.Val(i+1,j+2,k+2);
   real_t w2 = xs*(1-ys)*(1-zs)*f.Val(i+2,j+1,k+1)
             + xs*ys*(1-zs)*f.Val(i+2,j+2,k+1)
             + xs*(1-ys)*zs*f.Val(i+2,j+1,k+2)
             + xs*ys*zs*f.Val(i+2,j+2,k+2);
   return w1+w2;
}


real_t Tabulation::getValue(string funct,
                            real_t x,
                            real_t y,
                            real_t z,
                            real_t t)
{
   fct &f = Funct[_funct_id[funct]-1];
   size_t N1=f.Np[0], N2=f.Np[1], N3=f.Np[2], N4=f.Np[3];
   real_t h1 = (f.Max[0]-f.Min[0])/(N1-1),
          h2 = (f.Max[1]-f.Min[1])/(N2-1),
          h3 = (f.Max[2]-f.Min[2])/(N3-1),
          h4 = (f.Max[3]-f.Min[3])/(N4-1);
   size_t i=size_t((x-f.Min[0])/h1),
          j=size_t((y-f.Min[1])/h2),
          k=size_t((z-f.Min[2])/h3),
          l=size_t((t-f.Min[3])/h4);
   real_t xs = (x-f.Min[0])/h1 - i,
          ys = (y-f.Min[1])/h2 - j,
          zs = (z-f.Min[2])/h3 - k,
          ts = (t-f.Min[3])/h4 - l;
   if (x<f.Min[0])
      i=0, xs=0;
   else if (x>f.Max[0])
      i=N1-1, xs=1;
   if (y<f.Min[1])
      j=0, ys=0;
   else if (y>f.Max[1])
      j=N2-1, ys=1;
   if (z<f.Min[2])
      k=0, zs=0;
   else if (z>f.Max[2])
      k=N3-1, zs=1;
   if (t<f.Min[3])
      l=0, ts=0;
   else if (t>f.Max[3])
      l=N4-1, ts=1;
   real_t w1 = (1-xs)*(1-ys)*(1-zs)*(1-ts)*f.Val(i+1,j+1,k+1,l+1)
             + (1-xs)*(1-ys)*zs*(1-ts)*f.Val(i+1,j+1,k+2,l+1)
             + (1-xs)*(1-ys)*(1-zs)*ts*f.Val(i+1,j+1,k+1,l+2)
             + (1-xs)*(1-ys)*zs*ts*f.Val(i+1,j+1,k+2,l+2);
   real_t w2 = (1-xs)*ys*(1-zs)*(1-ts)*f.Val(i+1,j+2,k+1,l+1)
             + (1-xs)*ys*zs*(1-ts)*f.Val(i+1,j+2,k+2,l+1)
             + (1-xs)*ys*(1-zs)*ts*f.Val(i+1,j+2,k+1,l+2)
             + (1-xs)*ys*zs*ts*f.Val(i+1,j+2,k+2,l+2);
   real_t w3 = xs*(1-ys)*(1-zs)*(1-ts)*f.Val(i+2,j+1,k+1,l+1)
             + xs*(1-ys)*zs*(1-ts)*f.Val(i+2,j+1,k+2,l+1)
             + xs*(1-ys)*(1-zs)*ts*f.Val(i+2,j+1,k+1,l+2)
             + xs*(1-ys)*zs*ts*f.Val(i+2,j+1,k+2,l+2);
   real_t w4 = xs*ys*(1-zs)*(1-ts)*f.Val(i+2,j+2,k+1,l+1)
             + xs*ys*zs*(1-ts)*f.Val(i+2,j+2,k+2,l+1)
             + xs*ys*(1-zs)*ts*f.Val(i+2,j+2,k+1,l+2)
             + xs*ys*zs*ts*f.Val(i+2,j+2,k+2,l+2);
   return w1+w2+w3+w4;
}


real_t Tabulation::getDerivative(string funct,
                                 real_t x)
{
   fct &f = Funct[_funct_id[funct]-1];
   size_t N = f.Np[0];
   real_t h = (f.Max[0]-f.Min[0])/(N-1);
   size_t i = size_t((x-f.Min[0])/h);
   if (x<f.Min[0])
      i = 0;
   else if (x>f.Max[0])
      i = N - 1;
   return f.Val(i+2) - f.Val(i+1);
}


ostream& operator<<(ostream&          s,
                    const Tabulation& p)
{
   size_t n = p._nb_funct;
   s << "\n===============================================================================\n";
   s << "Summary of Tabulated function data:\n";
   s << "-------------------------------------------------------------------------------\n";
   for (size_t i=1; i<=n; ++i) {
      const fct &f = p.Funct[i-1];
      s << " Function: " << f.Label << ", Number of variables: " << f.NbVar << endl;
      for (size_t j=0; j<f.NbVar; ++j) {
         s << "    Variable: " << f.VarLabel[j] << ", Min. = " << f.Min[j];
         s << ", Max = " << f.Max[j] << ", Number of points = " << f.Np[j] << endl;
      }
      s << "-------------------------------------------------------------------------------" << endl;
   }
   return s;
}

} /* namespace OFELI */
