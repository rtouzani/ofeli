/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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
   _nb_funct++;
   _funct(_nb_funct).Label = label;
   _funct(_nb_funct).NbVar = 0;
   _funct_id[label] = _nb_funct;
}


void Tabulation::setVariable(string label)
{
   size_t n = _funct(_nb_funct).NbVar;
   _funct(_nb_funct).VarLabel[n] = label;
   _funct(_nb_funct).Id[label] = n + 1;
   _funct(_nb_funct).NbVar++;
}


void Tabulation::setSizes()
{
   if (_funct(_nb_funct).NbVar==1)
      _funct(_nb_funct).Val.setSize(_funct(_nb_funct).Np[0]);
   else if (_funct(_nb_funct).NbVar==2)
      _funct(_nb_funct).Val.setSize(_funct(_nb_funct).Np[0],_funct(_nb_funct).Np[1]);
   else if (_funct(_nb_funct).NbVar==3)
      _funct(_nb_funct).Val.setSize(_funct(_nb_funct).Np[0],_funct(_nb_funct).Np[1],_funct(_nb_funct).Np[2]);
}


real_t Tabulation::getValue(string funct,
                            real_t v)
{
   size_t n=_funct_id[funct], N=_funct(n).Np[0];
   real_t h = (_funct(n).Max[0]-_funct(n).Min[0])/(N-1);
   size_t i = size_t((v-_funct(n).Min[0])/h);
   real_t s = (v-_funct(n).Min[0])/h - i;
   if (v<_funct(n).Min[0])
      i=0, s=0;
   else if (v>_funct(n).Max[0])
      i=N-1, s=1;
   return (1-s)*_funct(n).Val(i+1) + s*_funct(n).Val(i+2);
}


real_t Tabulation::getValue(string funct,
                            real_t v1,
                            real_t v2)
{
   size_t n=_funct_id[funct];
   size_t N1=_funct(n).Np[0], N2=_funct(n).Np[1];
   real_t h1 = (_funct(n).Max[0]-_funct(n).Min[0])/(N1-1),
          h2 = (_funct(n).Max[1]-_funct(n).Min[1])/(N2-1);
   size_t i=size_t((v1-_funct(n).Min[0])/h1), j=size_t((v2-_funct(n).Min[1])/h2);
   real_t s = (v1-_funct(n).Min[0])/h1 - i,
          t = (v2-_funct(n).Min[1])/h2 - j;
   if (v1<_funct(n).Min[0])
      i=0, s=0;
   if (v1>_funct(n).Max[0])
      i=N1-1, s=1;
   if (v2<_funct(n).Min[1])
      j=0, t=0;
   if (v2>_funct(n).Max[1])
      j=N2-1, t=1;
   return   (1-s)*(1-t)*_funct(n).Val(i+1,j+1) + s*(1-t)*_funct(n).Val(i+2,j+1)
          + s*t*_funct(n).Val(i+2,j+2) + (1-s)*t*_funct(n).Val(i+1,j+2);
}


real_t Tabulation::getValue(string funct,
                            real_t v1,
                            real_t v2,
                            real_t v3)
{
   size_t n=_funct_id[funct];
   size_t N1=_funct(n).Np[0], N2=_funct(n).Np[1], N3=_funct(n).Np[2];
   real_t h1 = (_funct(n).Max[0]-_funct(n).Min[0])/(N1-1),
          h2 = (_funct(n).Max[1]-_funct(n).Min[1])/(N2-1),
          h3 = (_funct(n).Max[2]-_funct(n).Min[2])/(N3-1);
   size_t i=size_t((v1-_funct(n).Min[0])/h1), j=size_t((v2-_funct(n).Min[1])/h2), k=size_t((v3-_funct(n).Min[2])/h3);
   real_t r = (v1-_funct(n).Min[0])/h1 - i,
          s = (v2-_funct(n).Min[1])/h2 - j,
          t = (v3-_funct(n).Min[2])/h3 - k;
   if (v1<_funct(n).Min[0])
      i=0, r=0;
   if (v1>_funct(n).Max[0])
      i=N1-1, r=1;
   if (v2<_funct(n).Min[1])
      j=0, s=0;
   if (v2>_funct(n).Max[1])
      j=N2-1, s=1;
   if (v3<_funct(n).Min[2])
      k=0, t=0;
   if (v3>_funct(n).Max[2])
      k=N3-1, t=1;
   return   (1-r)*(1-s)*(1-t)*_funct(n).Val(i+1,j+1,k+1)
          + (1-r)*s*(1-t)*_funct(n).Val(i+1,j+2,k+1)
          + (1-r)*(1-s)*t*_funct(n).Val(i+1,j+1,k+2)
          + (1-r)*s*t*_funct(n).Val(i+1,j+2,k+2)
          + r*(1-s)*(1-t)*_funct(n).Val(i+2,j+1,k+1)
          + r*s*(1-t)*_funct(n).Val(i+2,j+2,k+1)
          + r*(1-s)*t*_funct(n).Val(i+2,j+1,k+2)
          + r*s*t*_funct(n).Val(i+2,j+2,k+2);
}


real_t Tabulation::getDerivative(string funct,
                                 real_t v)
{
   size_t n = _funct_id[funct];
   size_t N = _funct(n).Np[0];
   real_t h = (_funct(n).Max[0]-_funct(n).Min[0])/(N-1);
   size_t i = size_t((v-_funct(n).Min[0])/h);
   if (v<_funct(n).Min[0])
      i=0;
   else if (v > _funct(n).Max[0])
      i=N-1;
   return _funct(n).Val(i+2) - _funct(n).Val(i+1);
}


ostream& operator<<(      ostream&    s,
                    const Tabulation& p)
{
   size_t n = p._nb_funct;
   s << "\n===============================================================================\n";
   s << "Summary of Tabulated function data:\n";
   s << "-------------------------------------------------------------------------------\n";
   for (size_t i=1; i<=n; i++) {
      s << " Function: " << p._funct(i).Label << ", Number of variables: " << p._funct(i).NbVar << endl;
      for (size_t j=0; j<p._funct(i).NbVar; j++) {
         s << "    Variable: " << p._funct(i).VarLabel[j] << ", Min. = " << p._funct(i).Min[j];
         s << ", Max = " << p._funct(i).Max[j] << ", Number of points = " << p._funct(i).Np[j] << endl;
      }
   s << "-------------------------------------------------------------------------------" << endl;
   }
   return s;
}

} /* namespace OFELI */
