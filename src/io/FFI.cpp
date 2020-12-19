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

                      A class to read Free Format Files

  ==============================================================================*/

#include "io/FFI.h"
#include "util/util.h"
#include "OFELIException.h"

using std::ifstream;
using std::istringstream;

namespace OFELI {

FFI::FFI()
    : _msg(1), _string_nb(0), _ident(" "), _in(true), _non_fatal(false), _eol(true),
      _is(nullptr), _iss(nullptr)
{
}


FFI::FFI(const string& file)
    : _msg(2), _string_nb(0), _input_file(file), _ident(" "), _in(false), _non_fatal(false),
      _eol(true), _iss(nullptr)
{
   _in = false;
   _is = new ifstream;
   _is->open(file.c_str());
   if (_is->fail())
      throw OFELIException("In FFI::FFI(string): Trying to open file " + file + ". File not found.");
   *_is >> _ident;
}


FFI::FFI(const string& file,
         const string& ident)
    : _msg(1), _string_nb(1), _non_fatal(false), _eol(true)
{
   _in = false;
   _is = new ifstream;
   _is->open(file.c_str());
   if (_is->fail())
      throw OFELIException("FFI::FFI(string,string): Trying to open file " + file + ". File not found.");
   _input_file = file;
   _is->getline(_buffer,120);
   if (_is->fail())
      throw OFELIException("In FFI::FFI(string,string): End of file " + file + " reached.");
   _iss = new istringstream(string(_buffer),istringstream::in);
   *_iss >> _ident;
   if (_ident != ident)
      throw OFELIException("In FFI::FFI(string,string): Trying to read file " + file +
                           ". File must start with string: " + ident);
   delete _iss;
   _iss = nullptr;
}


FFI::FFI(const FFI& ff)
    : _is(ff._is), _iss(ff._iss)
{
}


FFI::~FFI()
{
   if (_is)
      delete _is;
   if (_iss)
      delete _iss;
}


void FFI::open(const string& file,
               const string& ident)
{
   _in = false;
   _is = new ifstream;
   _string_nb = 1;
   _is->open(file.c_str());
   if (_is->fail())
      throw OFELIException("In FFI::open(string,string): Trying to open file " +
                           file + ". File not found.");
   _input_file = file;
   _is->getline(_buffer,120);
   if (_iss)
      delete _iss;
   _iss = new istringstream(string(_buffer),istringstream::in);
   *_iss >> _ident;
   if (_ident != ident)
      throw OFELIException("In FFI::FFI(string,string): Trying to read file " + file +
                           ". File must start with string: " + ident);
   _eol = true;
   _msg = 1;
   _non_fatal = false;
   delete _iss;
   _iss = nullptr;
}


void FFI::open(const string& file,
               const string& ident,
               int           c)
{
   _in = false;
   _is = new ifstream;
   _string_nb = 1;
   _is->open(file.c_str());
   if (_is->fail())
      throw OFELIException("In FFI::open(string,string,int): Trying to open file "
                           + file + ". File not found.");
   _input_file = file;
   _is->getline(_buffer,120);
   if (_iss)
      delete _iss;
   _iss = new istringstream(string(_buffer),istringstream::in);
   *_iss >> _ident;
   if (_ident != ident)
      throw OFELIException("In FFI::open(string,string,int): Trying to read file " + file +
                           ". File must start with string: " + ident);
   _eol = true;
   _msg = c;
   _non_fatal = false;
   delete _iss;
   _iss = nullptr;
}


void FFI::open(const string& file)
{
   _in = false;
   _is = new ifstream;
   _string_nb = 1;
   _is->open(file.c_str());
   if (_is->fail())
      throw OFELIException("In FFI::open(string): Trying to open file " + file + ". File not found.");
   _input_file = file;
   *_is >> _ident;
   _msg = 1;
   _iss = nullptr;
}


void FFI::setKeywords(const vector<string>& s)
{
   _kw.resize(s.size());
   _kw = s;
}


int FFI::get_token()
{
   do {
      if (_eol) {
         if (_in)
            cin.getline(_buffer,120);
         else
            _is->getline(_buffer,120);
         Trim(_buffer);
         if (!_in)
            if (_is->eof())
               return 1;
         if (_iss)
            delete _iss;
         _iss = new istringstream(string(_buffer),istringstream::in);
      }
      _comment = false;
      *_iss >> _token;
      _eol = false;
      if (_iss->fail() == true)
         _eol = true;
      if (_token.c_str()[0]=='#') {
         _comment = true;
         _eol = true;
      }
   }
   while (_comment || _eol);
   return 0;
}


int FFI::getKW(const string& msg)
{
   if (msg != "0")
      cout << msg;
   string s = getS();
   for (size_t i=0; i<_kw.size(); i++) {
      int n = int(_kw[i].find("$"));
      if (_kw[i].substr(0,n) == s.substr(0,n))
         return int(i);
   }
   return -1;
}


int FFI::getI(const string& msg)
{
   if (msg != "0")
      cout << msg;
   if (get_token())
      throw OFELIException("In FFI::getI(string): End of file " +
                           _input_file + " reached.");
   int i = stringTo<int>(_token);
   _string_nb++;
   return i;
}


real_t FFI::getD(const string& msg)
{
   if (msg != "0")
      cout << msg;
   if (get_token())
      throw OFELIException("In FFI::getD(string): End of file "
                           + _input_file + " reached.");
   real_t d = stringTo<real_t>(_token);
   _string_nb++;
   return d;
}


string FFI::getS(const string& msg)
{
   if (msg != "0")
      cout << msg;
   if (get_token())
      throw OFELIException("In FFI::getS(string): End of file " + _input_file + " reached.");
   _string_nb++;
   return _token;
}


string FFI::getE(const string& msg)
{
   if (msg != "0")
      cout << msg;
   string s = getS();
   if (s[0]=='"' && s[s.size()-1]=='"')
      s[0] = s[s.size()-1] = ' ';
   return s;
}


FFI & FFI::operator=(const FFI& ff)
{
   _is = ff._is;
   _iss = ff._iss;
   return *this;
}

} /* namespace OFELI */
